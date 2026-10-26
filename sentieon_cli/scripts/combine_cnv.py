#!/usr/bin/env python3
"""
Combine CNVscope calls with SV calls converted to CNV by indel2cnv.py into one CNV VCF.

The rule values are DEFAULTS, those for paired-end CNVscope models, updated by the JSON file given with --preset
(presets/SE.json for single-end models), a file or a member of a Sentieon model bundle (MODEL.bundle/combine_cnv.json).
- Converted calls are written with SOURCE=SV_CNV, CN = CN_NEUTRAL + CNDIFF (floored at 0; on haploid sequence a
  homozygous call counts one copy) and FILTER PASS. Dropped: converted losses longer than max_converted_loss and,
  on haploid sequence, non-homozygous losses and gains longer than max_converted_gain_haploid.
- A converted gain whose array is at least nodepth_size long, whose array-average CN shift reaches
  nodepth_min_shift, and over which CNVscope's raw segmentation shows depth but almost no gain, gets FILTER NoDepth.
  The segments are the records of --raw (CNVscope's output) or, without --raw, every record of --cnv, which
  CNVModelApply then wrote with --emit_conf 0. The gain evidence is the segments' HMM state (HMMCN) where they carry
  it, else the PASS calls of --cnv.
- A CNVscope PASS call that repeats a kept converted call gets FILTER SVdup (INFO DEDUP names the array).
- With loss_support_size set (SE), a shorter CNVscope PASS loss that no converted loss overlaps gets FILTER NoSV.
- Other CNVscope PASS calls are written with SOURCE=CNV. Non-PASS CNVscope records (FILTER MLrejected: the segments
  the model rejected) are written as CNVModelApply wrote them, with SOURCE=CNV appended to INFO.
The final call set is the PASS records. CN_NEUTRAL follows CNVscope's --sex/--par or --blocks ploidy.
"""
import argparse
import bisect
import datetime
import io
import json
import logging
import math
import os
import re
import sys
from collections import defaultdict

import vcflib

VERSION = '2.3.0'  # combine_sv_cnv release this file belongs to

# The rule values for paired-end CNVscope models; a --preset JSON file overrides any of them (presets/SE.json)
DEFAULTS = {
    'max_converted_loss': 1500,          # converted losses longer than this are dropped (0 = keep all)
    'loss_support_size': 0,              # NoSV: CNVscope losses shorter than this need a converted loss (0 = off)
    'max_converted_gain_haploid': 3000,  # converted gains on haploid sequence longer than this are dropped
    'loss_min_overlap': 0.5,             # SVdup: share of a CNVscope loss inside one converted loss array
    'gain_min_overlap': 0.8,             # SVdup: share of a CNVscope gain inside one converted gain array,
    'gain_slack': 600,                   # or at most this many bp of it outside the array,
    'gain_min_shift': 0.4,               # and the CN shift the converted gains must explain over it
    'nodepth_size': 10000,               # NoDepth: converted gains with arrays at least this long are checked,
    'nodepth_min_shift': 0.5,            # if their array-average CN shift is visible to depth,
    'nodepth_min_segmented': 0.5,        # and raw segments of any state cover at least this share of the array;
    'nodepth_min_cov': 0.2,              # they fail when raw gain segments cover less than this share
}
FRACTIONS = ('loss_min_overlap', 'gain_min_overlap', 'nodepth_min_segmented', 'nodepth_min_cov')  # rules in [0, 1]
AR_MAGIC = b'!<arch>\n'  # a Sentieon model bundle is an ar archive

AUTOSOME_CN = 2
SEX_CHROMS = {'chrX': 'X', 'X': 'X', 'chrY': 'Y', 'Y': 'Y'}
log = logging.getLogger(__name__)

USAGE = ('combine_cnv.py --cnv VCF --converted VCF -o VCF [--raw VCF] [--preset FILE]\n'
         '                      [--sex {M,F} --par BED | --blocks BED] [-v]')
EXAMPLE = '''Production use, one sample:
  sentieon driver -r REF.fa --algo CNVModelApply --model MODEL --vcf S.cnv.raw.vcf.gz --emit_conf 0 S.cnv.vcf.gz
  python3 indel2cnv.py REF.fa S.sv.vcf.gz S.sv.cnv.vcf.gz -t 8
  python3 combine_cnv.py --preset presets/SE.json --cnv S.cnv.vcf.gz --converted S.sv.cnv.vcf.gz \\
      --sex M --par hs38.PAR.bed -o S.combined.cnv.vcf.gz
The final call set is the PASS records: bcftools view -f PASS S.combined.cnv.vcf.gz'''


def tool_version():
    """COMBINE_SV_CNV_VERSION when a driver set it, else combine_sv_cnv-v<VERSION>, followed by
    ' (sentieon-cli-<version>)' when this file sits in the sentieon_cli package (sentieon_cli/scripts/)."""
    if os.environ.get('COMBINE_SV_CNV_VERSION'):
        return os.environ['COMBINE_SV_CNV_VERSION']
    ver = 'combine_sv_cnv-v' + VERSION
    if os.path.basename(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))) == 'sentieon_cli':
        try:
            from sentieon_cli import __version__ as scli_version
        except Exception:
            try:
                from importlib.metadata import version
                scli_version = version('sentieon_cli')
            except Exception:
                scli_version = 'unknown'
        ver += ' (sentieon-cli-%s)' % scli_version
    return ver


def provenance_line(tool='combine_cnv'):
    """The ##CommandLine header line, in the form of Sentieon's ##SentieonCommandLine."""
    date = datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
    cmd = ' '.join(sys.argv).replace('\\', '\\\\').replace('"', '\\"')
    return '##CommandLine.%s=<ID=%s,Version="%s",Date="%s",CommandLine="%s">' % (tool, tool, tool_version(), date, cmd)


def fmt_value(v):
    return '%g' % v if isinstance(v, float) else str(v)


def describe(rules):
    return ', '.join(f'{k}={fmt_value(v)}' for k, v in rules.items())


def ar_member(archive, name):
    """The bytes of the member called name in an ar archive (a model bundle), reading only the member headers and
    that member. Member names as sentieon_cli.archive.ar_load reads them: GNU/SysV "name/" and "/<offset>" into the
    "//" table, BSD "#1/<length>" and plain names; the symbol tables are skipped. Raises ValueError naming the
    archive."""
    try:
        f = open(archive, 'rb')
    except OSError as e:
        raise ValueError(f'{archive}: cannot read ({e.strerror})') from None
    with f:
        if f.read(len(AR_MAGIC)) != AR_MAGIC:
            raise ValueError(f'{archive}: not an ar archive')
        file_size, want = os.fstat(f.fileno()).st_size, os.fsencode(name)
        names, table = [], None
        while True:
            hdr = f.read(60)
            if not hdr:
                break
            if len(hdr) < 60:
                raise ValueError(f'{archive}: truncated ar archive')
            field, size = hdr[:16].rstrip(b' '), hdr[48:58].strip()
            if hdr[58:] != b'`\n' or not size.isdigit():
                raise ValueError(f'{archive}: corrupted ar archive')
            start, size = f.tell(), int(size)
            if start + size > file_size:
                raise ValueError(f'{archive}: truncated ar archive')
            member = None
            if field in (b'/', b'/SYM64/'):  # symbol tables
                pass
            elif field == b'//':  # GNU long-name table
                table = f.read(size)
            elif field[:1] == b'/':  # GNU long name: /<offset> into the table
                off = field[1:]
                stop = table.find(b'/\n', int(off)) if table is not None and off.isdigit() else -1
                if stop < 0:
                    raise ValueError(f'{archive}: corrupted ar archive')
                member = table[int(off):stop]
            elif field[-1:] == b'/':  # GNU/SysV short name
                member = field[:-1]
            elif field[:3] == b'#1/':  # BSD long name: the first <length> bytes of the data
                n = field[3:]
                if not n.isdigit() or int(n) > size:
                    raise ValueError(f'{archive}: corrupted ar archive')
                member = f.read(int(n)).rstrip(b'\0')
            else:  # BSD short name
                member = field
            if member == want:
                return f.read(start + size - f.tell())
            if member is not None:
                names.append(member.decode('utf-8', 'replace'))
            f.seek(start + size + size % 2)
    raise ValueError(f'{archive}: no member "{name}"; members: {", ".join(names) or "none"}')


def read_preset(path):
    """The bytes of a preset: the file path, or member M of the ar archive (model bundle) A when path is A/M."""
    archive, member = os.path.split(path)
    if os.path.isfile(archive):
        return ar_member(archive, member)
    try:
        with open(path, 'rb') as f:
            head = f.read(len(AR_MAGIC))
            if head == AR_MAGIC:
                raise ValueError(f'{path}: an ar archive (model bundle); name a member, e.g. {path}/combine_cnv.json')
            return head + f.read()
    except (FileNotFoundError, NotADirectoryError):  # NotADirectoryError: the path runs through a file
        raise ValueError(f'{path}: no such file') from None
    except OSError as e:
        raise ValueError(f'{path}: cannot read ({e.strerror})') from None


def load_preset(path):
    """DEFAULTS updated by a preset: one JSON object of DEFAULTS keys and an optional "description" string, in a file
    or a member of a model bundle (MODEL.bundle/combine_cnv.json). Raises ValueError naming the path and the key."""
    def unique(pairs):
        keys = [k for k, _ in pairs]
        for k in keys:
            if keys.count(k) > 1:
                raise ValueError(f'{path}: duplicate key "{k}"')
        return dict(pairs)

    text = io.TextIOWrapper(io.BytesIO(read_preset(path)), encoding='utf-8')  # decoded as open(path) would
    try:
        data = json.load(text, object_pairs_hook=unique)
    except UnicodeDecodeError:
        raise ValueError(f'{path}: invalid JSON (not UTF-8 text)') from None
    except json.JSONDecodeError as e:
        raise ValueError(f'{path}: invalid JSON ({e})') from None
    if not isinstance(data, dict):
        raise ValueError(f'{path}: the top-level value is not a JSON object')
    rules = dict(DEFAULTS)
    for key, value in data.items():
        shown = json.dumps(value)
        if key == 'description':
            if not isinstance(value, str):
                raise ValueError(f'{path}: "description" is {shown}, not a string')
            continue
        if key not in DEFAULTS:
            raise ValueError(f'{path}: unknown key "{key}"; valid keys: description, {", ".join(DEFAULTS)}')
        is_int = isinstance(DEFAULTS[key], int)
        if isinstance(value, bool):
            raise ValueError(f'{path}: "{key}" is {shown}, a boolean, not a number')
        if is_int and not isinstance(value, int):
            raise ValueError(f'{path}: "{key}" is {shown}, not an integer')
        if not isinstance(value, (int, float)) or isinstance(value, float) and not math.isfinite(value):
            raise ValueError(f'{path}: "{key}" is {shown}, not a number')
        if key in FRACTIONS and not 0 <= value <= 1:
            raise ValueError(f'{path}: "{key}" is {shown}, outside [0, 1]')
        if value < 0:
            raise ValueError(f'{path}: "{key}" is {shown}, negative')
        rules[key] = value if is_int else float(value)
    return rules


def preset_rules(parser, path):
    """The rule values in effect and where they come from: DEFAULTS, updated by the --preset file when given."""
    if path is None:
        return dict(DEFAULTS), "the code's defaults (the PE values)"
    try:
        return load_preset(path), f'--preset {path}'
    except ValueError as e:
        parser.error(f'--preset {e}')


def parse_par_bed(path):
    """chrom -> [(start, end)] 0-based half-open, same as cnv_eval.py."""
    intervals = {}
    with open(path) as f:
        for line in f:
            flds = line.rstrip().split('\t')
            if not line.strip() or line.startswith('#') or len(flds) < 3:
                continue
            intervals.setdefault(flds[0], []).append((int(flds[1]), int(flds[2])))
    return intervals


def parse_blocks_bed(path):
    """CNVscope --blocks BED: chrom -> [(start, end, cn_neutral)] plus the sex= header hint."""
    blocks, sex = {}, None
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                m = re.search(r'\bsex=([MF])\b', line)
                if m:
                    sex = m.group(1)
                continue
            flds = line.rstrip().split('\t')
            if not line.strip() or len(flds) < 4:
                continue
            blocks.setdefault(flds[0], []).append((int(flds[1]), int(flds[2]), int(flds[3])))
    return blocks, sex


class Ploidy:
    """Neutral CN lookup mirroring CNVscope: the --blocks block covering pos, else the --sex/--par preset."""
    def __init__(self, sex, par_intervals, blocks=None):
        self.sex = sex
        self.par = par_intervals or {}
        self.blocks = {}
        for chrom, ivs in (blocks or {}).items():
            ivs = sorted(ivs)
            self.blocks[chrom] = ([s for s, _, _ in ivs], ivs)

    def cn_neutral(self, chrom, pos):
        """-1 = uncovered (female chrY), else the local neutral CN."""
        if chrom in self.blocks:
            starts, ivs = self.blocks[chrom]
            i = bisect.bisect_right(starts, pos) - 1
            if i >= 0 and pos < ivs[i][1]:
                return ivs[i][2]
        sc = SEX_CHROMS.get(chrom)
        if not self.sex or sc is None:
            return AUTOSOME_CN
        if self.sex == 'F':
            return -1 if sc == 'Y' else AUTOSOME_CN
        if sc == 'X' and any(s <= pos < e for s, e in self.par.get(chrom, ())):
            return AUTOSOME_CN
        return 1


def header_sex(vcf):
    for h in vcf.headers:
        if h.startswith('##SampleSex='):
            return h.split('=', 1)[1].strip()
    return None


def gt_alleles(v):
    """The sample's GT alleles as strings ('1|1' -> ['1', '1'], '1' -> ['1'])."""
    gt = v.samples[0].get('GT', '.') if v.samples else '.'
    return re.split(r'[/|]', str(gt))


def is_hom_alt(v):
    """Every allele non-reference and called: 1/1, 1|1, 1 (haploid), 1|2."""
    return all(a not in ('0', '.', '') for a in gt_alleles(v))


def over_loss_cap(v, max_loss):
    return v.info.get('CNDIFF', 0) < 0 and bool(max_loss) and v.end - v.pos > max_loss


def haploid_drop(v, cn0, max_gain_haploid):
    """Why a converted call on haploid sequence (CN_NEUTRAL 1) is dropped, or None."""
    if cn0 != 1:
        return None
    cndiff = v.info.get('CNDIFF', 0)
    # the SV caller genotypes male X/Y as diploid: a het deletion on single-copy sequence cannot be real
    if cndiff < 0 and not is_hom_alt(v):
        return 'haploid_het_loss'
    if cndiff > 0 and max_gain_haploid and v.end - v.pos > max_gain_haploid:
        return 'haploid_gain_oversize'
    return None


def converted_cn(v, cn0):
    """CN_NEUTRAL + CNDIFF floored at 0, and whether it was floored. On haploid sequence a hom call's CNDIFF counts
    the one copy twice (the converter sums both diploid haplotypes), so it is halved there."""
    cndiff = v.info.get('CNDIFF', 0)
    if cn0 == 1 and is_hom_alt(v) and len(gt_alleles(v)) > 1:
        half = (abs(cndiff) + 1) // 2      # an odd count (merged units) rounds away from zero
        cndiff = half if cndiff > 0 else -half
    return max(0, cn0 + cndiff), cn0 + cndiff < 0


def load_converted_intervals(converted_vcf, max_loss=0, ploidy=None, max_gain_haploid=0):
    """(chrom, direction) -> sorted (refstart, refstop, cndiff, svlen) of the converted calls; with ploidy, only the
    calls the size cap and the haploid rules keep."""
    intervals = defaultdict(list)
    for chrom in converted_vcf.contigs:
        for v in converted_vcf.range(chrom):
            cndiff = v.info.get('CNDIFF', 0)
            refstart = v.info.get('REFSTART', 0)
            refstop = v.info.get('REFSTOP', 0)
            if cndiff == 0 or refstart >= refstop or over_loss_cap(v, max_loss):
                continue
            if ploidy is not None and haploid_drop(v, ploidy.cn_neutral(chrom, v.pos), max_gain_haploid):
                continue
            intervals[(chrom, -1 if cndiff < 0 else 1)].append((refstart, refstop, cndiff, v.info.get('SVLEN', 0)))
    for key in intervals:
        intervals[key].sort()
    return intervals


def overlaps_any(intervals, pos, end, min_frac, slack=0):
    """The first interval holding >= min_frac of [pos, end), or leaving at most `slack` bp of it outside."""
    length = end - pos
    if length <= 0:
        return None
    for rs, re_, _, _ in intervals:
        if rs >= end:
            break
        if re_ <= pos:
            continue
        ovl = min(end, re_) - max(pos, rs)
        if ovl / length >= min_frac or length - ovl <= slack:
            return (rs, re_)
    return None


def explained_shift(intervals, pos, end):
    """Expected CN shift over [pos, end) from converted calls (CNDIFF*SVLEN spread over each array), and the
    largest contribution."""
    length = end - pos
    shift = 0.0
    best = None
    for rs, re_, cndiff, svlen in intervals:
        if rs >= end:
            break
        if re_ <= pos or re_ <= rs:
            continue
        part = cndiff * svlen * (min(end, re_) - max(pos, rs)) / ((re_ - rs) * length)
        shift += part
        if best is None or abs(part) > abs(best[2]):
            best = (rs, re_, part)
    return shift, best


def loss_duplicate(intervals, pos, end, rules):
    """The converted loss array holding loss_min_overlap of the CNVscope loss, or None."""
    return overlaps_any(intervals, pos, end, rules['loss_min_overlap'])


def gain_duplicate(intervals, pos, end, rules):
    """The converted gain array holding the CNVscope gain (gain_min_overlap, or gain_slack bp outside at most),
    when the converted gains also explain gain_min_shift of CN over it; else None."""
    hit = overlaps_any(intervals, pos, end, rules['gain_min_overlap'], rules['gain_slack'])
    if hit is None:
        return None
    shift, best = explained_shift(intervals, pos, end)
    return hit if best is not None and abs(shift) >= rules['gain_min_shift'] else None


def loss_supported(intervals, pos, end):
    """Does any converted loss array overlap [pos, end)?"""
    if end <= pos:
        return False
    for rs, re_, _, _ in intervals:
        if rs >= end:
            break
        if re_ > pos:
            return True
    return False


def load_raw_segments(records, chrom, ploidy):
    """CNVscope's raw segments among the records on chrom: -1 losses, 1 gains, 0 every segment (where CNVscope had
    depth to judge)."""
    segs = {-1: [], 0: [], 1: []}
    for v in records:
        cn0 = v.info.get('CN_NEUTRAL', ploidy.cn_neutral(chrom, v.pos))
        cn = v.info.get('HMMCN', cn0)
        if cn0 < 0:
            continue
        segs[0].append((v.pos, v.end, cn))
        if cn != cn0:
            segs[-1 if cn < cn0 else 1].append((v.pos, v.end, cn))
    for d in segs:
        segs[d].sort()
    return segs


def load_pass_calls(records, chrom, ploidy):
    """CNVscope's PASS calls among the records on chrom: -1 losses, 1 gains (the gain evidence for raw segments without
    HMM states)."""
    calls = {-1: [], 1: []}
    for v in records:
        if v.filter and set(v.filter).difference(('PASS',)):
            continue
        cn0 = v.info.get('CN_NEUTRAL', ploidy.cn_neutral(chrom, v.pos))
        cn = v.info.get('CN', cn0)
        if cn0 >= 0 and cn != cn0:
            calls[-1 if cn < cn0 else 1].append((v.pos, v.end, cn))
    for d in calls:
        calls[d].sort()
    return calls


def with_source(line):
    """A record line as it was, with SOURCE=CNV appended to its INFO."""
    f = line.split('\t')
    f[7] = 'SOURCE=CNV' if f[7] == '.' else f[7] + ';SOURCE=CNV'
    return '\t'.join(f)


def has_line(vcf, prefix):
    return any(h.startswith(prefix) for h in vcf.headers)


def emit_conf(vcf):
    """The --emit_conf value on the last CNVModelApply command line in the header, or None."""
    lines = [h for h in vcf.headers if h.startswith('##SentieonCommandLine.CNVModelApply=')]
    values = re.findall(r'--emit_conf[=\s]+([^\s">]+)', lines[-1]) if lines else []
    try:
        return float(values[-1]) if values else None
    except ValueError:
        return None


def has_hmm_states(raw_vcf):
    """Whether the raw segments carry their HMM state (HMMCN)."""
    return has_line(raw_vcf, '##INFO=<ID=HMMCN,')


def covered_fraction(segs, start, end):
    """Fraction of [start, end) covered by the (sorted, possibly overlapping) segments."""
    length = end - start
    if length <= 0:
        return 0.0
    covered = 0
    reach = start
    for cs, ce, _ in segs:
        if cs >= end:
            break
        if ce <= start:
            continue
        lo, hi = max(cs, reach), min(ce, end)
        if hi > lo:
            covered += hi - lo
            reach = hi
    return covered / length


def no_depth(cndiff, svlen, rs, re_, raw, rules):
    """A long converted gain over which CNVscope had depth but saw (almost) no gain."""
    if raw is None or cndiff <= 0 or re_ - rs < rules['nodepth_size']:
        return False
    if abs(cndiff * svlen) / (re_ - rs) < rules['nodepth_min_shift']:
        return False
    return (covered_fraction(raw[0], rs, re_) >= rules['nodepth_min_segmented']
            and covered_fraction(raw[1], rs, re_) < rules['nodepth_min_cov'])


def build_parser():
    parser = argparse.ArgumentParser(
        prog='combine_cnv.py', usage=USAGE, epilog=EXAMPLE, formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Combine CNVscope calls with SV calls converted to CNV (indel2cnv.py) into one CNV VCF.')
    req = parser.add_argument_group('required')
    req.add_argument('--cnv', required=True, metavar='VCF', help='CNVscope calls (CNVModelApply output)')
    req.add_argument('--converted', required=True, metavar='VCF', help='SV calls converted by indel2cnv.py')
    req.add_argument('-o', '--output', required=True, metavar='VCF', help='combined CNV VCF (.vcf.gz)')
    parser.add_argument('--raw', metavar='VCF',
                        help="CNVscope's own output (the input of CNVModelApply); default: the segments in --cnv, "
                             'which CNVModelApply must have written with --emit_conf 0')
    parser.add_argument('--preset', metavar='FILE',
                        help='JSON file of rule values, or a member of a model bundle, e.g. '
                             'MODEL.bundle/combine_cnv.json; without it the PE values apply (paired-end models). '
                             'presets/SE.json is for single-end models (cnv.se.model: Roche SBX, Ultima)')
    parser.add_argument('--sex', choices=['M', 'F'], help='sample sex (default: ##SampleSex in --cnv)')
    parser.add_argument('--par', metavar='BED', help='PAR BED, required for male samples, as for CNVscope')
    parser.add_argument('--blocks', metavar='BED', help='instead of --sex/--par: the CNVscope --blocks file')
    parser.add_argument('-v', '--verbose', action='store_true', help='log progress and counts')
    return parser


def setup(parser, args):
    """Check the inputs and resolve the ploidy as CNVscope does; returns (ploidy, sex)."""
    cnv_vcf = vcflib.VCF(args.cnv, 'r')
    raw_vcf = vcflib.VCF(args.raw, 'r') if args.raw is not None else None
    try:
        if not has_line(cnv_vcf, '##SentieonCommandLine.CNVModelApply='):
            parser.error(f'--cnv {args.cnv} is not CNVModelApply output')
        if raw_vcf is None:
            conf = emit_conf(cnv_vcf)
            if conf is None or conf > 0:
                parser.error(f'--cnv {args.cnv} was not written with CNVModelApply --emit_conf 0 and so lacks most '
                             "rejected segments: pass --raw (CNVscope's output) or rerun CNVModelApply with "
                             '--emit_conf 0')
        elif (not has_line(raw_vcf, '##SentieonCommandLine.CNVscope=')
                or has_line(raw_vcf, '##SentieonCommandLine.CNVModelApply=')):
            parser.error(f'--raw {args.raw} is not CNVscope\'s own output (the input of CNVModelApply)')
        if args.blocks and (args.sex or args.par):
            parser.error('--sex/--par and --blocks are mutually exclusive (as for CNVscope)')
        blocks, blocks_sex = parse_blocks_bed(args.blocks) if args.blocks else ({}, None)
        hsex = header_sex(cnv_vcf)
        if args.blocks and hsex:
            log.warning('--blocks given; ignoring ##SampleSex=%s in %s', hsex, args.cnv)
        sex = blocks_sex if args.blocks else (args.sex or hsex)
        if args.sex and hsex and args.sex != hsex:
            log.warning('--sex %s differs from ##SampleSex=%s in %s', args.sex, hsex, args.cnv)
        if not sex and not args.blocks:
            # read as diploid, a male's chrX/chrY would turn CNVscope's single-copy calls into losses
            parser.error(f'--cnv {args.cnv} names no sex (no ##SampleSex): give --sex (with --par for a male sample) '
                         'or --blocks')
        if sex == 'M' and not args.par and not args.blocks:
            parser.error('--sex M requires --par <par.bed>')
        par_intervals = parse_par_bed(args.par) if args.par else {}
        for chrom in par_intervals:
            if chrom not in cnv_vcf.contigs:
                parser.error(f'PAR BED contig {chrom!r} not in {args.cnv}')
        for chrom, ivs in blocks.items():
            if chrom not in cnv_vcf.contigs:
                parser.error(f'blocks BED contig {chrom!r} not in {args.cnv}')
            if any(cn < 1 or cn > 4 for _, _, cn in ivs):
                parser.error(f'--blocks cn_neutral must be in [1,4] ({chrom})')
    finally:
        cnv_vcf.close()
        if raw_vcf is not None:
            raw_vcf.close()
    if args.blocks:
        log.info('Ploidy: blocks=%s (%d blocks; uncovered sex chroms use preset %s)', args.blocks,
                 sum(len(b) for b in blocks.values()), sex or 'diploid')
    else:
        log.info('Ploidy: sex=%s par=%s', sex or 'unset (diploid)', args.par or 'none')
    return Ploidy(sex, par_intervals, blocks), sex


def output_headers(sex, tool='combine_cnv'):
    hdrs = [
        '##INFO=<ID=CN_NEUTRAL,Number=1,Type=Integer,Description="Copy number neutral state for this region">',
        '##INFO=<ID=SOURCE,Number=1,Type=String,Description="Call source: CNV or SV_CNV">',
        '##INFO=<ID=DEDUP,Number=1,Type=String,Description="Matching converted interval (chrom:refstart-refstop)">',
        '##FILTER=<ID=SVdup,Description="CNVscope call repeating a converted SV call in the same repeat array">',
        '##FILTER=<ID=NoDepth,Description="Long converted gain over which CNVscope saw depth but no gain">',
        '##FILTER=<ID=NoSV,Description="Short CNVscope loss that no converted SV loss overlaps">',
        provenance_line(tool),
    ]
    if sex:
        hdrs.append(f'##SampleSex={sex}')
    return hdrs


def combine(cnv_path, conv_path, raw_path, out_path, rules, ploidy, headers, on_converted=None):
    """Write the combined VCF and return the counts. Without raw_path, cnv_path (an --emit_conf 0 file) gives the
    segments. on_converted(chrom, v, rs, re_, d, raw) is called for every kept converted call before its filter is set
    (the dev tool annotates through it)."""
    conv_vcf = vcflib.VCF(conv_path, 'r')
    conv_intervals = load_converted_intervals(conv_vcf, rules['max_converted_loss'], ploidy,
                                              rules['max_converted_gain_haploid'])
    support = None
    if rules['loss_support_size']:
        # every converted loss counts, before the size cap and the haploid rules: did the SV caller see a loss here?
        support = load_converted_intervals(conv_vcf)
    conv_vcf.close()
    conv_vcf = vcflib.VCF(conv_path, 'r')
    cnv_vcf = vcflib.VCF(cnv_path, 'r')
    raw_vcf = vcflib.VCF(raw_path, 'r') if raw_path is not None else None
    seg_vcf = cnv_vcf if raw_vcf is None else raw_vcf  # without --raw the records of --cnv are the segments
    raw_states = has_hmm_states(seg_vcf)

    # header: CNVscope's, plus the converter's own lines (not its contigs)
    conv_extra = [h for h in conv_vcf.headers if h.startswith('##') and not h.startswith(('##contig=', '##fileformat='))]
    out_vcf = vcflib.VCF(out_path, 'wb')
    out_vcf.copy_header(cnv_vcf, update=tuple(conv_extra) + tuple(headers))
    out_vcf.emit_header()

    all_contigs = list(cnv_vcf.contigs.keys()) + [c for c in conv_vcf.contigs if c not in cnv_vcf.contigs]
    conv_n = dict.fromkeys(('kept', 'floored', 'uncovered', 'over_cap', 'nodepth', 'haploid_het_loss',
                            'haploid_gain_oversize'), 0)
    cnv_n = dict.fromkeys(('kept', 'svdup', 'svdup_gain', 'nosv', 'rejected', 'neutral', 'uncovered'), 0)
    seq = 0
    for chrom in all_contigs:
        records = []
        cnv_records = list(cnv_vcf.range(chrom))  # read once: segments, PASS-call evidence and output
        raw = None
        if chrom in seg_vcf.contigs:
            raw = load_raw_segments(cnv_records if raw_vcf is None else raw_vcf.range(chrom), chrom, ploidy)
            if not raw_states:
                raw.update(load_pass_calls(cnv_records, chrom, ploidy))

        for v in conv_vcf.range(chrom):
            if over_loss_cap(v, rules['max_converted_loss']):
                conv_n['over_cap'] += 1
                continue
            cn0 = ploidy.cn_neutral(chrom, v.pos)
            if cn0 < 0:
                conv_n['uncovered'] += 1
                continue
            why = haploid_drop(v, cn0, rules['max_converted_gain_haploid'])
            if why:
                conv_n[why] += 1
                continue
            cndiff = v.info.get('CNDIFF', 0)
            cn, floored = converted_cn(v, cn0)
            if floored:
                conv_n['floored'] += 1
                log.info('%s:%d CNDIFF=%d with CN_NEUTRAL=%d: CN floored to 0', chrom, v.pos + 1, cndiff, cn0)
            v.info['CN'] = cn
            v.info['CN_NEUTRAL'] = cn0
            v.info['SOURCE'] = 'SV_CNV'
            d = -1 if cndiff < 0 else 1
            rs, re_ = v.info.get('REFSTART', v.pos), v.info.get('REFSTOP', v.end)
            if on_converted:
                on_converted(chrom, v, rs, re_, d, raw)
            if no_depth(cndiff, v.info.get('SVLEN', 0), rs, re_, raw, rules):
                v.filter = ['NoDepth']
                conv_n['nodepth'] += 1
                # an unsupported gain must not suppress the CNVscope gain
                try:
                    conv_intervals[(chrom, d)].remove((rs, re_, cndiff, v.info.get('SVLEN', 0)))
                except (KeyError, ValueError):
                    pass
            elif not v.filter:
                v.filter = ['PASS']
            v.line = None
            conv_n['kept'] += 1
            seq += 1
            records.append((v.pos, 0, seq, v))

        for v in cnv_records:
            cn0 = v.info.get('CN_NEUTRAL', ploidy.cn_neutral(chrom, v.pos))
            if cn0 < 0:
                cnv_n['uncovered'] += 1
                continue
            cn = v.info.get('CN', cn0)
            if v.filter and set(v.filter).difference(('PASS',)):
                cnv_n['rejected'] += 1
                v.line = with_source(v.line)  # a segment the model rejected: CNVModelApply's own line
            elif cn == cn0:
                cnv_n['neutral'] += 1
                continue
            else:
                d = -1 if cn < cn0 else 1
                ivs = conv_intervals.get((chrom, d), [])
                hit = (loss_duplicate if d < 0 else gain_duplicate)(ivs, v.pos, v.end, rules)
                if hit:
                    v.filter = ['SVdup']
                    v.info['DEDUP'] = f'{chrom}:{hit[0]}-{hit[1]}'
                    cnv_n['svdup'] += 1
                    cnv_n['svdup_gain'] += d > 0
                elif (d < 0 and support is not None and v.end - v.pos < rules['loss_support_size']
                        and not loss_supported(support.get((chrom, -1), []), v.pos, v.end)):
                    v.filter = ['NoSV']
                    cnv_n['nosv'] += 1
                else:
                    cnv_n['kept'] += 1
                v.info['SOURCE'] = 'CNV'
                v.line = None
            seq += 1
            records.append((v.pos, 1, seq, v))

        records.sort()
        for _, _, _, v in records:
            out_vcf.emit(v)
        del records, cnv_records  # free this chromosome's records before the next is read

    out_vcf.close()
    conv_vcf.close()
    cnv_vcf.close()
    if raw_vcf is not None:
        raw_vcf.close()
    log.info('Converted: %d written (%d NoDepth gains, %d floored to CN 0); dropped %d over the loss cap, %d haploid '
             'het losses, %d haploid gains over the cap, %d on uncovered sequence', conv_n['kept'], conv_n['nodepth'],
             conv_n['floored'], conv_n['over_cap'], conv_n['haploid_het_loss'], conv_n['haploid_gain_oversize'],
             conv_n['uncovered'])
    log.info('CNVscope: %d PASS calls kept, %d SVdup (%d gains), %d NoSV, %d rejected segments written; left out %d '
             'neutral PASS records, %d on uncovered sequence', cnv_n['kept'], cnv_n['svdup'], cnv_n['svdup_gain'],
             cnv_n['nosv'], cnv_n['rejected'], cnv_n['neutral'], cnv_n['uncovered'])
    return conv_n, cnv_n


def main():
    parser = build_parser()
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING,
                        format='%(asctime)s %(levelname)s %(message)s')
    rules, source = preset_rules(parser, args.preset)
    log.info('%s, rules from %s: %s', tool_version(), source, describe(rules))
    ploidy, sex = setup(parser, args)
    combine(args.cnv, args.converted, args.raw, args.output, rules, ploidy, output_headers(sex))
    log.info('Written to %s', args.output)


if __name__ == '__main__':
    main()
