#!/usr/bin/env python3
"""
Combine CNVscope calls with indel2cnv-converted SV calls (both in CNV space).

The converted CNV calls have REFSTART/REFSTOP defining the repeat array extent.
CNVscope loss calls overlapping these ranges are marked FILTER=SVdup (converted
call is preferred — it has precise breakpoints from the SV caller).
Gain calls are never marked.

Both call sets use CN for evaluation: CNVscope has CN natively, converted calls
get CN = neutral_cn + CNDIFF.

Usage:
    python3 combine_cnv.py \\
        --cnv test/pred.HG002_30x_hs38.vcf.gz \\
        --converted test/hg002-svcaller.cnv.vcf.gz \\
        -o test/combined_cnv.vcf.gz
"""
import argparse
import sys
from collections import defaultdict

import vcflib


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def load_converted_intervals(converted_vcf):
    """Load REFSTART/REFSTOP intervals from converted calls, keyed by (chrom, direction)."""
    intervals = defaultdict(list)
    n = 0
    for chrom in converted_vcf.contigs:
        for v in converted_vcf.range(chrom):
            cndiff = v.info.get('CNDIFF', 0)
            refstart = v.info.get('REFSTART', 0)
            refstop = v.info.get('REFSTOP', 0)
            if cndiff == 0 or refstart >= refstop:
                continue
            d = -1 if cndiff < 0 else 1
            intervals[(chrom, d)].append((refstart, refstop))
            n += 1
    for key in intervals:
        intervals[key].sort()
    return intervals, n


def overlaps_any(intervals, pos, end, min_frac):
    """Check if [pos, end) overlaps any interval by >= min_frac of the call."""
    length = end - pos
    if length <= 0:
        return None
    for rs, re_ in intervals:
        if rs >= end:
            break
        if re_ <= pos:
            continue
        ovl = min(end, re_) - max(pos, rs)
        if ovl / length >= min_frac:
            return (rs, re_)
    return None


def main():
    parser = argparse.ArgumentParser(
        description='Combine CNVscope and converted SV calls (both CNV space)')
    parser.add_argument('--cnv', required=True, help='CNVscope call VCF')
    parser.add_argument('--converted', required=True,
        help='indel2cnv.py output VCF (converted SV calls)')
    parser.add_argument('-o', '--output', required=True, help='Output VCF (.vcf.gz)')
    parser.add_argument('--min_overlap', type=float, default=0.5,
        help='Min overlap fraction to mark SVdup (default: 0.5)')
    parser.add_argument('--neutral_cn', type=int, default=2,
        help='Neutral copy number (default: 2)')
    args = parser.parse_args()

    conv_vcf = vcflib.VCF(args.converted, 'r')
    cnv_vcf = vcflib.VCF(args.cnv, 'r')

    # Build dedup intervals from converted calls
    eprint('Loading converted intervals...')
    conv_intervals, n_conv = load_converted_intervals(conv_vcf)
    eprint(f'  {n_conv} intervals')

    # Reopen for full iteration
    conv_vcf.close()
    conv_vcf = vcflib.VCF(args.converted, 'r')

    # Merged header: CNV as base, add converted-specific headers
    conv_extra = [h for h in conv_vcf.headers if h.startswith('##')]
    new_hdrs = (
        '##INFO=<ID=SOURCE,Number=1,Type=String,Description="Call source: CNV or SV_CNV">',
        '##INFO=<ID=DEDUP,Number=1,Type=String,Description="Matching converted interval (chrom:refstart-refstop)">',
        '##FILTER=<ID=SVdup,Description="CNV call redundant with converted SV call in same repeat array">',
    )
    out_vcf = vcflib.VCF(args.output, 'wb')
    out_vcf.copy_header(cnv_vcf, update=tuple(conv_extra) + new_hdrs)
    out_vcf.emit_header()

    all_contigs = list(cnv_vcf.contigs.keys())
    for c in conv_vcf.contigs:
        if c not in cnv_vcf.contigs:
            all_contigs.append(c)

    conv_stats = {'total': 0}
    cnv_stats = {'kept': 0, 'marked': 0, 'neutral': 0, 'filtered': 0}

    seq = 0
    for chrom in all_contigs:
        records = []

        # Converted SV calls: emit all, add CN from CNDIFF
        for v in conv_vcf.range(chrom):
            cndiff = v.info.get('CNDIFF', 0)
            v.info['CN'] = args.neutral_cn + cndiff
            v.info['SOURCE'] = 'SV_CNV'
            v.line = None
            conv_stats['total'] += 1
            seq += 1
            records.append((v.pos, 0, seq, v))

        # CNVscope calls: mark losses overlapping converted REFSTART/REFSTOP
        for v in cnv_vcf.range(chrom):
            if v.filter and set(v.filter).difference(('PASS',)):
                cnv_stats['filtered'] += 1
                continue
            cn = v.info.get('CN', args.neutral_cn)
            if cn == args.neutral_cn:
                cnv_stats['neutral'] += 1
                continue
            v.info['SOURCE'] = 'CNV'
            if cn < args.neutral_cn:
                key = (chrom, -1)
                interval = overlaps_any(conv_intervals.get(key, []),
                                        v.pos, v.end, args.min_overlap)
                if interval:
                    v.filter = ['SVdup']
                    v.info['DEDUP'] = f'{chrom}:{interval[0]}-{interval[1]}'
                    cnv_stats['marked'] += 1
                    v.line = None
                    seq += 1
                    records.append((v.pos, 1, seq, v))
                    continue
            cnv_stats['kept'] += 1
            v.line = None
            seq += 1
            records.append((v.pos, 1, seq, v))

        records.sort()
        for _, _, _, v in records:
            out_vcf.emit(v)

    out_vcf.close()
    conv_vcf.close()
    cnv_vcf.close()

    eprint(f'Converted: {conv_stats["total"]} records')
    eprint(f'CNV: {cnv_stats["kept"]} kept, {cnv_stats["marked"]} marked SVdup, '
           f'{cnv_stats["neutral"]} neutral, {cnv_stats["filtered"]} non-PASS')
    eprint(f'Written to {args.output}')


if __name__ == '__main__':
    main()
