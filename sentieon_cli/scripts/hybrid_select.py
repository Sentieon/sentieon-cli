#!/usr/bin/env python
from __future__ import print_function
import argparse
import math
import multiprocessing as mp
import sys
import vcflib
import os
from vcflib.compat import *

class HybridFilter:

    params = {
        'min_conf_longread' : [ 30., "Minimum call confidence in long read", ""],
        'min_depth_longread' : [ 2, "Minimum depth in long read", ""],
        'min_conf_shortread' : [ 30., "Minimum call confidence in short read", ""],
    }

    @classmethod
    def add_arguments(cls, parser):
        for k,v in cls.params.items():
            parser.add_argument('--'+k, default=v[0], type=type(v[0]), help=v[1] + ' (default: '+str(v[0])+')', metavar=v[2])

    def __init__(self, args):
        self.args = args

    def applyFilters(self, in_vcf, out_vcf):
        start = getattr(in_vcf, 'start', -1)
        for v in in_vcf:
            filters = []
            lad = v.samples[0].get('LAD')
            sad = v.samples[0].get('SAD')
            lpl = v.samples[0].get('LPL')
            spl = v.samples[0].get('SPL')

            if not lad:
                lad=[0,0,0]
            if not lpl:
                lpl=[0,0,0]
            if not spl:
                spl=[0,0,0]

            lminpos = lpl.index(0) #same as GT
            sminpos = spl.index(0)
            lpl_ref = lpl[0]
            lpl.remove(0)
            lpl_conf = min(lpl) #2nd min
            spl.remove(0)
            spl_conf = min(spl) #2nd min

            if not filters:
                if sum(lad) < self.args.min_depth_longread:
                    filters.append('longread_lowdepth')

            if not filters:
                if ((v.info.get('STR') is None) and (lminpos != 0)):
                    lpl_conf = lpl_ref
                if ( lpl_conf < self.args.min_conf_longread):
                    filters.append('longread_lowconf')

            if not filters:
                if ( spl_conf >= self.args.min_conf_shortread):
                    if (lminpos == sminpos):
                        filters.append('same_gt')
                    elif (lminpos == 0) and (v.info.get('STR') is not None):
                        filters.append('longread_str')

            flds = v.line.split('\t')
            flds[6] = filters and ';'.join(sorted(set(filters))) or 'PASS'
            v.line = '\t'.join(flds)
            if v.pos >= start:
                out_vcf.emit(v)
        return

extra_headers = (
    '##FILTER=<ID=PASS,Description="Accept as a confident somatic mutation">',
    '##FILTER=<ID=longread_lowdepth,Description="low depth in long read sample">',
    '##FILTER=<ID=longread_lowconf,Description="low confidence in long read sample">',
    '##FILTER=<ID=shortread_lowconf,Description="low confidence in short read sample">',
    '##FILTER=<ID=same_gt,Description="same gt in long/sort read sample">',
    '##FILTER=<ID=longread_str,Description="str region in long read sample">',
)

remove_headers = (
    '##source=*',
)

expect_types = {
    'GT': {'Number': '1', 'Type': 'String' },
}

def check_header_types(vcf):
    vcf_types = {}
    for k, v in iteritems(vcf.infos):
        if k in expect_types:
            dd = {}
            for kk, vv in iteritems(v):
                if kk in ('Number', 'Type'):
                    dd[kk] = vv
                vcf_types[k] = dd
    for k, v in iteritems(vcf.formats):
        if k in expect_types:
            dd = {}
            for kk, vv in iteritems(v):
                if kk in ('Number', 'Type'):
                    dd[kk] = vv
                vcf_types[k] = dd
    return vcf_types == expect_types

def main(args):
    if not os.path.exists(args.vcf):
        print('Error: input file %s does not exist' %  args.vcf)
        return -1

    invcf = vcflib.VCF(args.vcf, 'r')

    if not check_header_types(invcf):
        print('Error: vcf format is not expected')
        return -1

    filter = HybridFilter(args)

    outvcf = vcflib.VCF(args.output, 'w')
    outvcf.copy_header(invcf, extra_headers, remove_headers)
    outvcf.emit_header()

    if args.threads < 2:
        filter.applyFilters(invcf, outvcf)
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
        sharder.run(shards, filter.applyFilters, [], invcf, outvcf)

    outvcf.close()
    invcf.close()
    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='sentieon pyexec hybrid_select.py', usage='%(prog)s [options] -v VCF output.vcf.gz')
    parser.add_argument('output', help='Output vcf file name')
    parser.add_argument('-v','--vcf',required=True, help='Input vcf file name')
    parser.add_argument('-t','--threads',type=int,default=mp.cpu_count(),help='number of threads')
    parser.add_argument('--step-size',type=int,default=10*1000*1000,help=argparse.SUPPRESS)
    HybridFilter.add_arguments(parser)
    sys.exit(main(parser.parse_args()))

# vim: ts=4 sw=4 expandtab
