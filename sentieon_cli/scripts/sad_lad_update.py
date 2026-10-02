#!/usr/bin/env python3
"""sad_lad_update.py -- alignment-free pangenome AD/DP update.

The dnascope-pangenome joint pileup calls a bwa CRAM (all reads) together with
a minimap2 pangenome-realignment CRAM (a subset of the same reads, tagged
LR:1).  Reads that appear in both CRAMs are counted twice, so FORMAT/AD =
SAD + LAD, FORMAT/DP and INFO/DP are inflated.  Plain VCF and gVCF input are
both accepted.

This script picks, per record and per sample, EITHER the short-read depth
vector (FORMAT/SAD) OR the pangenome one (FORMAT/LAD) as the updated
FORMAT/AD -- whichever likelihood vector (SPL / LPL) better supports the
called genotype -- and derives FORMAT/DP and INFO/DP from that choice.  No
alignment file is read; only the fields already present in the VCF are used.
gVCF hom-ref blocks are the one exception: nothing there is chosen, read or
written (see "gVCF input" below).

Decision rule, per record and sample.  Let g be the PL index of the called
genotype (VCF ordering) and, for a likelihood vector X,

    conf(X) = min(X[j] for j != g) - X[g],   capped at --conf_cap

so conf is GQ-like: large positive means X strongly favours the called GT, a
negative value means X prefers a different genotype.  conf is None
("uninformative") when X is missing, has fewer than two entries, or has no
value at index g.

  branch 1  GT missing (./., .) or SAD or LAD missing
            -> sample left untouched, ADSRC=NONE, ADRULE=0
  branch 2  sum(LAD) == 0                        -> SAD   (nothing to gain;
            sum(SAD) == 0 and sum(LAD) > 0       -> LAD    this also covers
            every LPL=. record)
  branch 3  exactly one conf is None             -> the informative one
            both None -> SAD if sum(SAD) >= sum(LAD) else LAD
  branch 4  conf_S - conf_L >  --similar_threshold -> SAD
            conf_L - conf_S >  --similar_threshold -> LAD
  branch 5  "similar": LAD, UNLESS sum(SAD) > (1 + --lad_ratio) * sum(LAD),
            in which case SAD.

ADRULE is written as 0 for branch 1 (untouched) as declared in the header, and
as 2..5 for the branches that actually choose a source; the value 1 is never
emitted.

Output.  FORMAT/AD becomes the chosen vector, FORMAT/DP its sum, INFO/DP the
sum of the new FORMAT/DP over the samples of the record.  The originals are
kept as FORMAT/AD_ORIG, FORMAT/DP_ORIG and INFO/DP_ORIG, and the choice is
recorded as FORMAT/ADSRC (SAD/LAD/NONE) and FORMAT/ADRULE.  For an untouched
sample AD and DP keep their original values and AD_ORIG / DP_ORIG are written
as '.' -- the '.' is the signal that nothing was updated there.

gVCF input.  A record whose ALT is exactly <NON_REF> (or <*>) is a
homozygous-reference block and is copied to the output byte for byte: its
original line is emitted, nothing is parsed back out, and it is counted only
as records_ref_block -- never in samples_seen, the by-filter tables or
records_updated / records_untouched.  No update is possible there
(a block carries no SAD/LAD, so no source can be chosen) and nothing may be
added either: INFO/DP or extra FORMAT keys on a block change what GVCFtyper
reads back, so the safe action and the right one are the same.  gVCF
*variant* records are handled exactly like plain-VCF ones -- their trailing
<NON_REF> allele is present in AD/SAD/LAD (Number=R) and in PL/SPL
(Number=G) alike, so the vectors stay parallel and gt_index() / conf_of() /
select() need no gVCF special case.  Output naming needs none either:
'out.g.vcf.gz' satisfies the --output_vcf '.vcf.gz' suffix check.

Implementation.  Built on Sentieon's vcflib: records are read as parsed
Variant objects, v.info / v.samples are mutated in place, v.line is cleared
and vcflib's VCF.format() re-emits the record.  Re-formatting is therefore
expected and intended -- format() sorts the INFO keys, writes FORMAT as GT
followed by the remaining keys in sorted order, prints QUAL as '%4.2f' and
re-encodes floats through str() (QUAL 0 -> 0.00, FS=0 -> FS=0.0).  Values are
unchanged; only their spelling is.  The header declarations are added through
vcflib's copy_header(update=...), so the new lines land at the end of their
##INFO / ##FORMAT group.  vcflib also writes the output .tbi during the write
pass; it answers region queries like a tabix index but carries no htslib count
pseudo-bin, so `bcftools index -n` reports 0 on it.

Usage:
  sad_lad_update.py --input_vcf in.vcf.gz --output_vcf out.vcf.gz \\
                 [--stats stats.tsv] [--similar_threshold 10] \\
                 [--lad_ratio 0.10] [--conf_cap 99] [--threads N]

With --threads > 1 vcflib buffers each shard in a temporary file under
$SENTIEON_TMPDIR (or the system temp dir); point it at a writable filesystem
with room for the uncompressed output (~4 GB for a WGS VCF).
"""

import argparse
import json
import resource
import sys
import time

try:
    import vcflib
except ImportError:
    sys.stderr.write(
        "sad_lad_update: cannot import vcflib.\n"
        "  Install it with `pip install vcflib`, or point PYTHONPATH at a\n"
        "  vcflib checkout, e.g.\n"
        "    PYTHONPATH=/path/to/vcflib python3 sad_lad_update.py ...\n")
    sys.exit(1)

INFO_DP_ORIG_HDR = ('##INFO=<ID=DP_ORIG,Number=1,Type=Integer,Description='
                    '"Original combined depth before the pangenome AD/DP '
                    'update">')
FMT_AD_ORIG_HDR = ('##FORMAT=<ID=AD_ORIG,Number=R,Type=Integer,Description='
                   '"Original allelic depths before the pangenome AD/DP '
                   'update">')
FMT_DP_ORIG_HDR = ('##FORMAT=<ID=DP_ORIG,Number=1,Type=Integer,Description='
                   '"Original read depth before the pangenome AD/DP '
                   'update">')
FMT_ADSRC_HDR = ('##FORMAT=<ID=ADSRC,Number=1,Type=String,Description='
                 '"Source of the updated AD: SAD, LAD or NONE">')
FMT_ADRULE_HDR = ('##FORMAT=<ID=ADRULE,Number=1,Type=Integer,Description='
                  '"Rule branch that selected ADSRC (1-5 as in '
                  'sad_lad_update.py; 0 = untouched)">')

SRC_NONE, SRC_SAD, SRC_LAD = 0, 1, 2
SRC_NAME = ('NONE', 'SAD', 'LAD')
FILTER_CLASSES = ('PASS', 'MLrejected', 'LowQual', 'other')

COUNTERS = (
    'records',
    'records_ref_block',
    'records_updated',
    'records_untouched',
    'records_no_samples',
    'samples_seen',
    'samples_updated',
    'samples_untouched',
    'src_SAD', 'src_LAD', 'src_NONE',
    'rule_0', 'rule_1', 'rule_2', 'rule_3', 'rule_4', 'rule_5',
    'info_DP_added',
    'fmt_AD_added',
    'fmt_DP_added',
    'AD_changed',
    'FMT_DP_changed',
    'INFO_DP_changed',
)
BYFILT_KEYS = ('records', 'src_SAD', 'src_LAD', 'src_NONE',
               'rule_0', 'rule_1', 'rule_2', 'rule_3', 'rule_4', 'rule_5')


# ------------------------------------------------------------ decision rule --
def gt_index(gt):
    """PL index of the called genotype, or None if GT is missing.

    VCF ordering: for alleles sorted ascending a_0..a_{P-1},
    index = sum_k C(a_{k-1} + k - 1, k).  Diploid a<=b gives b(b+1)/2 + a and
    haploid a gives a.
    """
    if not gt or gt == '.' or gt == './.' or gt == '.|.':
        return None
    parts = gt.split('|') if '|' in gt else gt.split('/')
    try:
        alleles = sorted(int(p) for p in parts)
    except ValueError:          # any missing allele -> treat GT as missing
        return None
    if not alleles:
        return None
    if len(alleles) == 1:
        return alleles[0]
    if len(alleles) == 2:
        a, b = alleles
        return (b * (b + 1)) // 2 + a
    idx = 0                     # general polyploid fallback
    for k, a in enumerate(alleles, start=1):
        n, num, den = a + k - 1, 1, 1
        for i in range(k):
            num *= (n - i)
            den *= (i + 1)
        idx += num // den
    return idx


def conf_of(pl, g, cap):
    """min(pl[j], j != g) - pl[g], capped; None when uninformative."""
    if pl is None:
        return None
    n = len(pl)
    if n < 2 or g >= n:
        return None
    x = pl[g]
    if x is None:
        return None
    best = None
    for j in range(n):
        if j == g:
            continue
        y = pl[j]
        if y is not None and (best is None or y < best):
            best = y
    if best is None:
        return None
    c = best - x
    return c if c < cap else cap


def sum_ad(v):
    s = 0
    for x in v:
        if x is not None:
            s += x
    return s


def select(g, sad, lad, spl, lpl, cap, thr, ratio):
    """-> (source, branch).  See the module docstring."""
    if g is None or sad is None or lad is None:
        return SRC_NONE, 0
    ssum = sum_ad(sad)
    lsum = sum_ad(lad)
    if lsum == 0:
        return SRC_SAD, 2
    if ssum == 0:
        return SRC_LAD, 2
    cs = conf_of(spl, g, cap)
    cl = conf_of(lpl, g, cap)
    if cs is None or cl is None:
        if cs is None and cl is None:
            return (SRC_SAD if ssum >= lsum else SRC_LAD), 3
        return (SRC_LAD if cs is None else SRC_SAD), 3
    if cs - cl > thr:
        return SRC_SAD, 4
    if cl - cs > thr:
        return SRC_LAD, 4
    if ssum > (1.0 + ratio) * lsum:
        return SRC_SAD, 5
    return SRC_LAD, 5


def filt_class(filters):
    """v.filter is a list; [] means '.'"""
    if not filters or filters == ['PASS']:
        return 'PASS'
    if filters == ['MLrejected']:
        return 'MLrejected'
    if filters == ['LowQual']:
        return 'LowQual'
    return 'other'


# ----------------------------------------------------------------- the map --
def fix_shard(vcfi, vcfo, cap=99.0, thr=10.0, ratio=0.10):
    """Rewrite one shard (or, at --threads 1, the whole file).

    vcfi/vcfo are VCF or VCFReader/VCFWriter.  Returns (counters, by_filter)
    which the parent reduces in shard order.
    """
    cnt = dict.fromkeys(COUNTERS, 0)
    byfilt = dict((c, dict.fromkeys(BYFILT_KEYS, 0)) for c in FILTER_CLASSES)
    # VCFWriter.emit drops v.pos >= end but NOT v.pos < start, so a record
    # that straddles the shard start belongs to the previous shard and would
    # otherwise be written twice (and make tabix.add assert).
    lo = getattr(vcfo, 'start', 0)

    for v in vcfi:
        if v.pos < lo:
            continue
        cnt['records'] += 1
        if v.alt == ['<NON_REF>'] or v.alt == ['<*>']:
            # gVCF hom-ref block: no SAD/LAD to choose between, and adding
            # INFO/DP or FORMAT keys would change what GVCFtyper reads back.
            cnt['records_ref_block'] += 1
            vcfo.emit(v)                        # v.line is intact
            continue
        if not v.samples:
            cnt['records_no_samples'] += 1
            vcfo.emit(v)                        # v.line is intact
            continue

        bf = byfilt[filt_class(v.filter)]
        bf['records'] += 1

        info = v.info
        has_info_dp = 'DP' in info
        old_info_dp = info.get('DP')
        total = 0
        any_upd = False

        for s in v.samples:
            cnt['samples_seen'] += 1
            sad = s.get('SAD')
            lad = s.get('LAD')
            src, rule = select(gt_index(s.get('GT')), sad, lad,
                               s.get('SPL'), s.get('LPL'), cap, thr, ratio)
            has_ad = 'AD' in s
            has_dp = 'DP' in s
            old_ad = s.get('AD')
            old_dp = s.get('DP')
            if not has_ad:
                cnt['fmt_AD_added'] += 1
            if not has_dp:
                cnt['fmt_DP_added'] += 1

            if src == SRC_NONE:
                cnt['samples_untouched'] += 1
                cnt['src_NONE'] += 1
                bf['src_NONE'] += 1
                s['AD_ORIG'] = None
                s['DP_ORIG'] = None
                if old_dp is not None:
                    total += old_dp
            else:
                cnt['samples_updated'] += 1
                any_upd = True
                new_ad = list(sad if src == SRC_SAD else lad)
                new_dp = sum_ad(new_ad)
                if src == SRC_SAD:
                    cnt['src_SAD'] += 1
                    bf['src_SAD'] += 1
                else:
                    cnt['src_LAD'] += 1
                    bf['src_LAD'] += 1
                s['AD_ORIG'] = old_ad if has_ad else None
                s['DP_ORIG'] = old_dp if has_dp else None
                s['AD'] = new_ad
                s['DP'] = new_dp
                total += new_dp
                if new_ad != old_ad:
                    cnt['AD_changed'] += 1
                if new_dp != old_dp:
                    cnt['FMT_DP_changed'] += 1
            s['ADSRC'] = SRC_NAME[src]
            s['ADRULE'] = rule
            cnt['rule_%d' % rule] += 1          # rule_0 .. rule_5
            bf['rule_%d' % rule] += 1

        if any_upd:
            cnt['records_updated'] += 1
        else:
            cnt['records_untouched'] += 1

        if has_info_dp:
            info['DP_ORIG'] = old_info_dp
        else:
            cnt['info_DP_added'] += 1
        if total != old_info_dp:
            cnt['INFO_DP_changed'] += 1
        info['DP'] = total

        v.line = None                           # force VCF.format() to re-emit
        vcfo.emit(v)

    return cnt, byfilt


# -------------------------------------------------------------- the driver --
def sharded_run(nthr, contigs, func, step, *args, **kwargs):
    """vcflib.Sharder driver, after sentieon_cli.vcf_mod.sharded_run()."""
    sharder = vcflib.Sharder(nthr)
    intervals = ((c, 0, int(t['length'])) for c, t in contigs.items())
    shards = sharder.cut(intervals, step)
    return sharder.run(shards, func, [], *args, **kwargs)


def open_vcfs(input_vcf, output_vcf):
    """after sentieon_cli.vcf_mod.open_vcfs(); the header is copied from the
    input with the five new declarations added through vcflib."""
    try:
        vcfi = vcflib.VCF(input_vcf, 'r')
    except EnvironmentError as err:
        sys.exit('sad_lad_update: cannot read %s or its index: %s' %
                 (input_vcf, err))
    if not vcfi.contigs:
        sys.exit('sad_lad_update: the input VCF header has no ##contig lines; '
                 'vcflib cannot index the output without them')
    for k in ('SAD', 'LAD'):
        if k not in vcfi.formats:
            sys.stderr.write('sad_lad_update: WARNING the input does not declare '
                             'FORMAT/%s; every sample will be left '
                             'untouched\n' % k)

    if ('DP_ORIG' in vcfi.infos or 'AD_ORIG' in vcfi.formats or
            'DP_ORIG' in vcfi.formats):
        sys.stderr.write(
            'sad_lad_update: WARNING the input already declares DP_ORIG and/or '
            'AD_ORIG (has this VCF already been updated?); the *current* '
            'DP/AD are what will be stored as _ORIG\n')
    update = [INFO_DP_ORIG_HDR, FMT_AD_ORIG_HDR, FMT_DP_ORIG_HDR,
              FMT_ADSRC_HDR, FMT_ADRULE_HDR]

    try:
        vcfo = vcflib.VCF(output_vcf, 'w')
        vcfo.copy_header(vcfi, update=update)
        vcfo.emit_header()
    except EnvironmentError as err:
        sys.exit('sad_lad_update: cannot write %s: %s' % (output_vcf, err))
    return vcfi, vcfo


def write_stats(path, cnt, byfilt, args, elapsed, rss_mb, nshards):
    meta = {
        'input_vcf': args.input_vcf,
        'output_vcf': args.output_vcf,
        'similar_threshold': args.similar_threshold,
        'lad_ratio': args.lad_ratio,
        'conf_cap': args.conf_cap,
        'threads': args.threads,
        'shards': nshards,
        'elapsed_s': round(elapsed, 1),
        'peak_rss_mb': round(rss_mb, 1),
    }
    if path.endswith('.json'):
        with open(path, 'w') as fh:
            json.dump({'params': meta, 'counters': cnt, 'by_filter': byfilt},
                      fh, indent=1, sort_keys=True)
            fh.write('\n')
        return
    with open(path, 'w') as fh:
        fh.write('#class\tkey\tvalue\n')
        for k in sorted(meta):
            fh.write('param\t%s\t%s\n' % (k, meta[k]))
        for k in sorted(cnt):
            fh.write('total\t%s\t%d\n' % (k, cnt[k]))
        for c in FILTER_CLASSES:
            for k in sorted(byfilt[c]):
                fh.write('filter:%s\t%s\t%d\n' % (c, k, byfilt[c][k]))


def main(args):
    if not args.output_vcf.endswith('.vcf.gz'):
        sys.exit('sad_lad_update: --output_vcf must end in .vcf.gz')
    t0 = time.time()
    vcfi, vcfo = open_vcfs(args.input_vcf, args.output_vcf)

    kw = dict(cap=args.conf_cap, thr=args.similar_threshold,
              ratio=args.lad_ratio)
    if args.threads > 1:
        results = sharded_run(args.threads, vcfi.contigs, fix_shard,
                              args.step_size, vcfi, vcfo, **kw)
    else:
        results = [fix_shard(vcfi, vcfo, **kw)]

    total = dict.fromkeys(COUNTERS, 0)
    byfilt = dict((c, dict.fromkeys(BYFILT_KEYS, 0)) for c in FILTER_CLASSES)
    for cnt, bf in results:
        for k, n in cnt.items():
            total[k] += n
        for c in FILTER_CLASSES:
            for k, n in bf[c].items():
                byfilt[c][k] += n

    for f in (vcfi, vcfo):
        f.close()

    elapsed = time.time() - t0
    rss = (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss +
           resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss) / 1024.0
    for k in COUNTERS:
        sys.stderr.write('sad_lad_update: %-22s= %d\n' % (k, total[k]))
    sys.stderr.write('sad_lad_update: %-22s= %d\n' % ('shards', len(results)))
    sys.stderr.write('sad_lad_update: %-22s= %.1f\n' % ('elapsed_s', elapsed))
    sys.stderr.write('sad_lad_update: %-22s= %.1f\n' % ('peak_rss_mb', rss))
    if args.stats:
        write_stats(args.stats, total, byfilt, args, elapsed, rss,
                    len(results))
    return 0


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__.split('\n\n')[1],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--input_vcf', required=True,
                   help='input VCF or gVCF (bgzf + index); gVCF hom-ref '
                        'blocks are passed through unchanged')
    p.add_argument('--output_vcf', required=True,
                   help='output VCF; must end in .vcf.gz (so .g.vcf.gz is '
                        'accepted), .tbi is written too')
    p.add_argument('--stats', default=None,
                   help='write the counters here (.json -> JSON, else TSV)')
    p.add_argument('--similar_threshold', type=float, default=10.0,
                   help='|conf_S - conf_L| at or below which the two sources '
                        'are called "similar" (branch 4)')
    p.add_argument('--lad_ratio', type=float, default=0.10,
                   help='branch 5 uses LAD unless sum(SAD) exceeds sum(LAD) '
                        'by more than this fraction')
    p.add_argument('--conf_cap', type=float, default=99.0,
                   help='cap applied to conf(SPL) and conf(LPL)')
    p.add_argument('--threads', type=int, default=1,
                   help='worker processes; 1 runs in-process, >1 shards')
    p.add_argument('--step_size', type=int, default=25 * 1000 * 1000,
                   help='shard size in bp (only used when --threads > 1)')
    return p.parse_args(argv)


if __name__ == '__main__':
    sys.exit(main(parse_args()))

# vim: ts=4 sw=4 expandtab
