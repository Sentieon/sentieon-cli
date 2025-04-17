#!/usr/bin/env python
from __future__ import print_function
import argparse
import math
import multiprocessing as mp
import sys
import vcflib
import os
import bisect
from vcflib.compat import *


extra_headers = (
    '##INFO=<ID=LHC,Number=1,Type=Integer,Description="Longread hap count">',
)

remove_headers = (
)

def annotate_lhc(invcf, outvcf, ctg2s_e_c):
    start = getattr(invcf, 'start', -1)
    bed_ctg = ''
    bed_idx = 0
    for v in invcf:
        ctg = v.chrom
        pos = v.pos
        hap_cnt = -1
        if (bed_ctg != ctg) :
            if ctg in ctg2s_e_c:
                bed_ctg = ctg
                bed_idx = bisect.bisect_left(ctg2s_e_c[ctg], (pos, pos+1, 0))
                if (bed_idx > 0):
                    bed_idx -= 1
        if (bed_ctg == ctg) :
            while (bed_idx < len(ctg2s_e_c[ctg])):
                if (ctg2s_e_c[ctg][bed_idx][1] <= pos):
                    bed_idx += 1
                else :
                    break
            if (bed_idx < len(ctg2s_e_c[ctg])):
                if (ctg2s_e_c[ctg][bed_idx][0] <= pos):
                    hap_cnt = ctg2s_e_c[ctg][bed_idx][2]

        if hap_cnt != -1:
            cols = v.line.split('\t')
            cols[7] += ';LHC='+str(hap_cnt)
            v.line='\t'.join(cols)
        if v.pos >= start:
            outvcf.emit(v)

def main(args):
    if not os.path.exists(args.bed):
        print('Error: input bed file %s does not exist' % args.bed)
        return -1

    ctg2s_e_c=dict()

    with open(args.bed, 'r') as bedf:
        for line in bedf:
            cols = line.rstrip().split('\t')
            if len(cols) < 4:
                print('Error: wrong format in line %s' % line)
                return -1
            ctg = cols[0]
            if ctg not in ctg2s_e_c:
                ctg2s_e_c[ctg] = list()
            ctg2s_e_c[ctg].append((int(cols[1]), int(cols[2]), int(cols[3])))

    if not os.path.exists(args.vcf):
        print('Error: input file %s does not exist' %  args.vcf)
        return -1

    invcf = vcflib.VCF(args.vcf, 'r')

    outvcf = vcflib.VCF(args.output, 'w')
    outvcf.copy_header(invcf, extra_headers, remove_headers)
    outvcf.emit_header()

    if args.threads < 2:
        annotate_lhc(invcf, outvcf, ctg2s_e_c)
    else:
        sharder = vcflib.Sharder(args.threads)
        try:
            contig_lengths = [
                (contig, 0, int(d["length"]))
                for contig, d in invcf.contigs.items()
            ]
        except KeyError:
            return (
                "This script requires a VCF with contig lengths in"
                " the header when using multiple threads"
            )
        shards = sharder.cut(contig_lengths, args.step_size)
        sharder.run(shards, annotate_lhc, [], invcf, outvcf, ctg2s_e_c)

    outvcf.close()
    invcf.close()
    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='sentieon pyexec hybrid_anno.py', usage='%(prog)s [options] -v VCF -b BED output')
    parser.add_argument('output', help='Output vcf file name')
    parser.add_argument('-v','--vcf',required=True, help='Input vcf file name')
    parser.add_argument('-b','--bed',required=True, help='region haplotype count from long reads')
    parser.add_argument('-t','--threads',type=int,default=mp.cpu_count(),help='number of threads')
    parser.add_argument('--step-size',type=int,default=10*1000*1000,help=argparse.SUPPRESS)
    sys.exit(main(parser.parse_args()))

# vim: ts=4 sw=4 expandtab
