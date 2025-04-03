#!/usr/bin/env python

"""
Functionality for manipulating DNAscope-LR VCFs
"""

# Copyright (c) 2023-2025 Sentieon Inc. All rights reserved

from __future__ import print_function
import argparse
import bisect
import copy
import functools
import heapq
import io
import itertools
import multiprocessing
import operator
import os
import re
import resource
import sys
import time
import vcflib

from vcflib.compat import *

class IntervalList(object):
    def __init__(self, path):
        self.regions = self.load_bed(path)

    def load_bed(self, path):
        regions = {}
        if path.endswith('.gz'):
            fp = bgzf.open(path, 'rb')
        else:
            fp = io.open(path, 'rb')
        for line in fp:
            if line.startswith(b'#') or line.startswith(b'track'):
                continue
            cols = line.rstrip().decode().split('\t')
            if len(cols) < 3:
                continue
            chrom,start,end = cols[:3]
            regions.setdefault(chrom,[]).append([int(start),int(end)])
        fp.close()
        for chrom in regions:
            v = []
            s0 = e0 = None
            for (s,e) in sorted(regions[chrom]):
                if e0 is None:
                    s0 = s
                    e0 = e
                elif s > e0:
                    v.extend((s0,e0))
                    s0 = s
                    e0 = e
                else:
                    e0 = e
            if e0 is not None:
                v.extend((s0,e0))
            regions[chrom] = v
        return regions

    def get(self, c, s, e):
        r = self.regions.get(c)
        v = []
        if r is not None:
            i = bisect.bisect_right(r, s)
            if i % 2:
                v.append((r[i-1], r[i]))
                i += 1
            while i < len(r) and r[i] < e:
                v.append((r[i], r[i+1]))
                i += 2
        return v

    def __contains__(self, cs):
        if isinstance(cs, basestring):
            return cs in self.regions
        c, s = cs
        r = self.regions.get(c)
        return r and bisect.bisect_right(r, s) % 2 != 0

def grouper(*vcfs):
    start = -1
    q = []
    for k, vcf in enumerate(vcfs):
        if vcf is None:
            continue
        i = iter(vcf)
        v = next(i, None)
        if v:
            e = (v.pos+1 + len(v.info['RU'])*v.info['RPA'][0]
                if v.info.get('STR') else v.end)
            heapq.heappush(q, (v.pos, e, k, v, i))
        if start < 0:
            start = getattr(vcf, 'start', start)
    while q:
        grp = [None] * len(vcfs)
        ovl = []
        pos, end, k, v, i = heapq.heappop(q)
        grp[k] = v
        v = i and next(i, None) or None
        if v:
            e = (v.pos+1 + len(v.info['RU'])*v.info['RPA'][0]
                if v.info.get('STR') else v.end)
            heapq.heappush(q, (v.pos, e, k, v, i))
        while q and q[0][0] >= pos and q[0][0] < end:
            _, e, k, v, i = heapq.heappop(q)
            end = max(end, e)
            if v.pos > pos:
                ovl.append((k,v))
            elif grp[k] is None:
                grp[k] = v
            else:
                ovl.append((k,v))
            v = i and next(i, None) or None
            if v:
                e = (v.pos+1 + len(v.info['RU'])*v.info['RPA'][0]
                    if v.info.get('STR') else v.end)
                heapq.heappush(q, (v.pos, e, k, v, i))
        if pos >= start:
            yield (pos, grp, ovl)
        for k,v in ovl:
            if v and v.pos > pos:
                e = (v.pos+1 + len(v.info['RU'])*v.info['RPA'][0]
                    if v.info.get('STR') else v.end)
                heapq.heappush(q, (v.pos, e, k, v, None))

def trim(r, a):
    i, n = 1, min(len(r), len(a))
    while i < n and r[-i] == a[-i]:
        i += 1
    return (r[:len(r)-i+1], a[:len(a)-i+1])

def combine(r1, a1, r2, a2):
    if len(r1) == len(r2):
        return r1, [a1, a2]
    if len(r1) < len(r2):
        return r2, [a1 + r2[len(r1):], a2]
    else:
        return r1, [a1, a2 + r1[len(r2):]]

def compatible(v, d):
    for a in v.alt:
        if all(trim(d.ref, b) != trim(v.ref, a) for b in d.alt):
            return False
    return True


def sharded_run(
    nthr, vcf_contigs, func, step=100 * 1000 * 1000, bed=None, *args
):
    sharder = vcflib.Sharder(nthr)
    contigs = ((c, 0, int(t["length"])) for c, t in iteritems(vcf_contigs))
    if bed:
        contigs = ((c, s, e) for c, s, e in contigs if c in bed)
    shards = sharder.cut(contigs, step)
    return sharder.run(shards, func, None, *args)

def open_vcfs(input_vcf_fns, output_vcf_fns, update=None, remove=None, hdr_idx=0):
    input_vcfs = []
    for vcf_fn in input_vcf_fns:
        try:
            input_vcfs.append(vcflib.VCF(vcf_fn, "r"))
        except EnvironmentError as err:
            print(str(err))
            print(
                "Error: failed to read the input VCF or VCF index file - '{}'".format(
                    vcf_fn
                )
            )
            sys.exit(1)

    output_vcfs = []
    for i, vcf_fn in enumerate(output_vcf_fns):
        try:
            out_vcf = vcflib.VCF(vcf_fn, "w")
            _hdr_idx = i if hdr_idx is None else hdr_idx
            out_vcf.copy_header(input_vcfs[_hdr_idx], update=update, remove=remove)
            out_vcf.emit_header()
            output_vcfs.append(out_vcf)
        except EnvironmentError as err:
            print(str(err))
            print(
                "Error: failed to open the output VCF or VCF index file - '{}'".format(
                    vcf_fn
                )
            )
            sys.exit(1)

    return (input_vcfs, output_vcfs)

def haploid_patch_main(args):
    input_vcfs, output_vcfs = open_vcfs(
        (args.vcf, args.vcf_hp),
        (args.patch_vcf,),
        update=(
            '##INFO=<ID=DELTA,Number=0,Type=Flag,Description="Delta flag">',
        ),
        hdr_idx=None,
    )
    vcfs = input_vcfs + output_vcfs
    result = sharded_run(
        args.thread_count,
        input_vcfs[0].contigs,
        haploid_patch,
        100 * 1000 * 1000,
        None,
        *vcfs
    )

    for f in vcfs:
        if f:
            f.close()
    return 0

def haploid_patch(vcfi1, vcfd1, vcfo1, **kwargs):
    for pos, grp, ovl in grouper(vcfi1, vcfd1):
        v1, d1 = grp

        if d1 and d1.samples[0].get('GT') is None:
            d1 = None
        if not v1 or d1 and compatible(v1, d1):
            v1 = d1
        if not v1:
            continue

        v1 = trim1(v1 == d1 and vcfd1 or vcfi1, v1)
        i1, p1 = getpl(v1, None)

        # skip high conf ref calls
        if (i1 == 0 and p1 == 0 and 
            (not v1 or v1.samples[0].get('GQ') >= 20)):
            continue

        if v1:
            if v1 == d1 or v1.info.get('STR') and v1.info.get('RPA')[0]>=4:
                v1.info['DELTA'] = True
            v1.info.pop('ML_PROB', None)
            v1.filter = []
            v1.line = None
            vcfo1.emit(v1)

def merge_main(args):
    input_vcfs, output_vcfs = open_vcfs(
        (args.hap1, args.hap2, args.unphased, args.phased),
        (args.output_vcf,),
        update=(
            '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">',
        ),
        remove=(
            '##FORMAT=<ID=PID>',
            '##FORMAT=<ID=PGT>',
        ),
        hdr_idx=3,
    )
    bed = IntervalList(args.bed)
    merge_args = input_vcfs + output_vcfs + [bed]
    result = sharded_run(
        args.thread_count,
        input_vcfs[0].contigs,
        merge2,
        100 * 1000 * 1000,
        bed,
        *merge_args
    )

    for f in input_vcfs + output_vcfs:
        if f:
            f.close()
    return 0

def sub1(f, v, *args):
    v.alt = [v.alt[i] for i in args]
    n = len(os.path.commonprefix([v.ref[:0:-1]] + [a[:0:-1] for a in v.alt]))
    if n > 0:
        v.ref, v.alt = v.ref[:-n], [a[:-n] for a in v.alt]
        v.end = v.pos + len(v.ref)
    n = {}
    for k,u in iteritems(v.info):
        d = f.infos.get(k)
        t = d and d['Number'] or None
        if t == 'A':
            n[k] = [u[i] for i in args]
        elif t == 'R' or t == 'G':
            n[k] = [u[0]] + [u[i+1] for i in args]
    if n:
        v.info.update(n)
    n = {}
    for k,u in iteritems(v.samples[0]):
        d = f.formats.get(k)
        t = d and d['Number'] or None
        if t == 'A':
            n[k] = [u[i] for i in args]
        elif t == 'R' or t == 'G':
            n[k] = [u[0]] + [u[i+1] for i in args]
        if k == 'GT':
            gt = u and int(u) or 0
            if gt > 0:
                try:
                    gt = args.index(gt-1) + 1
                except:
                    gt = 0
                n[k] = str(gt)
    if n:
        v.samples[0].update(n)
    v.line = None
    return v

def trim1(f, v):
    if v is None or len(v.alt) <= 2:
        return v
    pl = v.samples[0].get('PL')
    if pl:
        pl = sorted(enumerate(pl[1:]), key=operator.itemgetter(1))
        pl = (i for i,_ in pl[:2])
    else:
        pl = (0, 1)
    return sub1(f, v, *pl)

def getpl(v, g):
    if v is None:
        if g is None:
            return (0, 0)
    else:
        pl = v.samples[0].get('PL')
        if pl and len(v.alt) > 0:
            if g is None:
                return (0, pl[0])
            for i,a in enumerate(v.alt):
                if trim(v.ref, a) == g:
                    return (i+1, pl[i+1])
    return (-1, 999999999)

def split1(f, v):
    gt = re.split(r'([/|])', v.samples[0].get('GT'))
    assert len(gt) == 3
    if gt[1] != '|' and gt[0] != gt[2] or gt[0] == '.' or gt[2] == '.':
        # workaround phaser not marking homvar phased
        gt[0] = gt[2] = 0

    for i in (int(gt[0]), int(gt[2])):
        if i == 0:
            u = None
        else:
            ia, ig = (0, i), (0, i*(i+1)//2)
            u = copy.deepcopy(v)
            sub2(f, u, ia, ig)
            u.samples[0]['GT'] = '1'
        yield u

def patch1_main(args):
    if (not args.phased and (not args.hap1 or not args.hap2)) or (
        args.phased and (args.hap1 or args.hap2)
    ):
        print(
            "Either the --phased VCF or both the --hap1 and --hap2 VCFs are"
            " required"
        )
        sys.exit(1)

    input_vcfs, output_vcfs = [], []
    if args.phased:
        input_vcfs, output_vcfs = open_vcfs(
            (args.phased, args.hap1_hp, args.hap2_hp),
            (args.patch1, args.patch2),
            update=(
                '##INFO=<ID=DELTA,Number=0,Type=Flag,Description="Delta flag">',
            ),
        )
        input_vcfs.insert(1, None)
    else:
        input_vcfs, output_vcfs = open_vcfs(
            (args.hap1, args.hap2, args.hap1_hp, args.hap2_hp),
            (args.patch1, args.patch2),
            update=(
                '##INFO=<ID=DELTA,Number=0,Type=Flag,Description="Delta flag">',
            ),
            hdr_idx=None,
        )

    vcfs = input_vcfs + output_vcfs
    result = sharded_run(
        args.thread_count,
        input_vcfs[0].contigs,
        patch1,
        100 * 1000 * 1000,
        None,
        *vcfs
    )

    for f in vcfs:
        if f:
            f.close()
    return 0

def patch1(vcfi1, vcfi2, vcfd1, vcfd2, vcfo1, vcfo2, **kwargs):
    for pos, grp, ovl in grouper(vcfi1, vcfi2, vcfd1, vcfd2):
        v1, v2, d1, d2 = grp

        if vcfi2 is None and v1:
            v1, v2 = split1(vcfi1, v1)

        if d1 and d1.samples[0].get('GT') is None:
            d1 = None
        if d2 and d2.samples[0].get('GT') is None:
            d2 = None
        if not v1 or d1 and compatible(v1, d1):
            v1 = d1
        if not v2 or d2 and compatible(v2, d2):
            v2 = d2
        if not v1 and not v2:
            continue

        v1 = trim1(v1 == d1 and vcfd2 or vcfi1, v1)
        v2 = trim1(v2 == d2 and vcfd2 or vcfi2 or vcfi1, v2)

        i1, p1 = getpl(v1, None)
        i2, p2 = getpl(v2, None)

        # skip high conf ref calls
        if (i1 == 0 and i2 == 0 and p1 == 0 and p2 == 0 and
            (not v1 or v1.samples[0].get('GQ') >= 20) and
            (not v2 or v2.samples[0].get('GQ') >= 20)):
            continue

        if v1:
            if v1 == d1 or vcfi2 and v1.info.get('STR') and v1.info.get('RPA')[0]>=4:
                v1.info['DELTA'] = True
            v1.info.pop('ML_PROB', None)
            v1.filter = []
            v1.line = None
            vcfo1.emit(v1)
        if v2:
            if v2 == d2 or vcfi2 and v2.info.get('STR') and v2.info.get('RPA')[0]>=4:
                v2.info['DELTA'] = True
            v2.info.pop('ML_PROB', None)
            v2.filter = []
            v2.line = None
            vcfo2.emit(v2)

def sub2(f, v, ia, ig):
    n = min(len(v.ref), min(len(v.alt[i-1]) for i in ia[1:]))
    v.ref = v.ref[:len(v.ref)-n+1]
    v.alt = [v.alt[i-1][:len(v.alt[i-1])-n+1] for i in ia[1:]]
    v.end = v.pos + len(v.ref)
    v.line = None

    n = {}
    for k,u in iteritems(v.info):
        d = f.infos.get(k)
        t = d and d['Number'] or None
        if t == 'A':
            n[k] = [u[i-1] for i in ia[1:]]
        elif t == 'R':
            n[k] = [u[i] for i in ia]
    if n:
        v.info.update(n)
    n = {}
    for k,u in iteritems(v.samples[0]):
        d = f.formats.get(k)
        t = d and d['Number'] or None
        if t == 'A':
            n[k] = [u[i-1] for i in ia[1:]]
        elif t == 'R':
            n[k] = [u[i] for i in ia]
        elif t == 'G':
            n[k] = [u[i] for i in ig]
    if n:
        v.samples[0].update(n)

def trim2(f, v, j1, j2):
    if j2 < j1:
        j1, j2 = j2, j1
    if v is None or len(v.alt) == 1:
        return (v, j1, j2)

    gt = [(i1,i2) for i2 in range(0,len(v.alt)+1) for i1 in range(0,i2+1)]
    pl = v.samples[0].get('PL')
    if pl is None or len(pl) != len(gt):
        return (None, 0, 0)

    p,i,i1,i2 = min((pl[k],k,k1,k2) for k,(k1,k2) in enumerate(gt) if k > 0)
    if i1 == 0 or i1 == i2:
        p,i,i1,i3 = min((pl[k],k,k1,k2) for k,(k1,k2) in enumerate(gt)
            if k > 0 and k != i and (k1 in (0,i2) or k2 in (0,i2)))
    else:
        i3 = 0
    ia = [i for i,_ in itertools.groupby(sorted((0, i1, i2, i3)))]
    ig = [i for i,(i1,i2) in enumerate(gt) if i1 in ia and i2 in ia]
    try:
        j1, j2 = ia.index(j1), ia.index(j2)
    except ValueError:
        j1, j2 = -1, -1
    v.samples[0]['GT'] = len(ia) == 3 and '1/2' or '1/1'
    sub2(f, v, ia, ig)
    return (v, j1, j2)

def patch2_main(args):
    input_vcfs, output_vcfs = open_vcfs(
        (args.vcf, args.vcf_hp),
        (args.output_vcf,),
        update=(
            '##INFO=<ID=DELTA,Number=0,Type=Flag,Description="Delta flag">',
        ),
    )
    vcfs = input_vcfs + output_vcfs
    result = sharded_run(
        args.thread_count,
        input_vcfs[0].contigs,
        patch2,
        100 * 1000 * 1000,
        None,
        *vcfs
    )

    for f in vcfs:
        if f:
            f.close()
    return 0

def patch2(vcfi, vcfd, vcfo, **kwargs):
    for pos, grp, ovl in grouper(vcfi, vcfd):
        v, d = grp

        if d and d.samples[0].get('GT') == './.':
            d = None
        if not v or d and compatible(v, d):
            v = d
        if not v:
            continue

        if vcfd:
            i1, p1 = getpl(v, None)
            i2, p2 = getpl(v, None)
            v, i1, i2 = trim2(v == d and vcfd or vcfi, v, i1, i2)
        else:
            i1, p1 = -1, 0
            i2, p2 = -1, 0

        # skip high conf ref calls
        if (i1 == 0 and i2 == 0 and p1 == 0 and p2 == 0 and
            v == d and d and d.samples[0].get('GQ') >= 20):
            continue

        if v:
            if v == d:
                v.info['DELTA'] = True
            v.info.pop('ML_PROB', None)
            v.filter = []
            v.line = None
            vcfo.emit(v)

def join2(f0, f1, f2, v0, v1, v2, pos, bed):
    if v1:
        gt = None if v1.filter else v1.samples[0].get('GT')
        i1 = int(gt) if gt else 0
        r1, a1 = trim(v1.ref, v1.alt[i1-1] if i1 else '')
    else:
        i1 = -1
    if v2:
        gt = None if v2.filter else v2.samples[0].get('GT')
        i2 = int(gt) if gt else 0
        r2, a2 = trim(v2.ref, v2.alt[i2-1] if i2 else '')
    else:
        i2 = -1

    if i1 <= 0 and i2 <= 0:
        return None

    v = vcflib.Variant(f0.chrom, pos, '.', '.', [], 0, [], {}, [{}])
    if v0:
        v.info.update(v0.info)

    v.id = v.id == '.' and v0 and v0.id or v.id
    v.id = v.id == '.' and v1 and v1.id or v.id
    v.id = v.id == '.' and v2 and v2.id or v.id

    ps = bed and bed.get(v.chrom, pos, pos) or None
    if ps:
        v.samples[0]['PS'] = ps[0][0]+1

    if i1 > 0 and i2 > 0:
        if r1 == r2 and a1 == a2:
            v.ref, v.alt = r1, [a1]
            v.samples[0]['GT'] = '1|1' if ps else '1/1'
            v.info['AC'] = [2]
            v.info['AF'] = [1.0]
        else:
            v.ref, v.alt = combine(r1, a1, r2, a2)
            v.samples[0]['GT'] = '1|2' if ps else '1/2'
            v.info['AC'] = [1,1]
            v.info['AF'] = [0.5,0.5]
        v.qual = v1.qual + v2.qual
        v.info['AN'] = 2
    else:
        if i1 > 0:
            v.ref, v.alt = r1, [a1]
            v.samples[0]['GT'] = '1|0' if ps else '0/1'
            v.qual = v1.qual
        else:
            v.ref, v.alt = r2, [a2]
            v.samples[0]['GT'] = '0|1' if ps else '0/1'
            v.qual = v2.qual
        v.info['AC'] = [1]
        v.info['AF'] = [0.5]
        v.info['AN'] = 2

    # allele map
    m0 = v0 and [0] + [next((i+1 for i,b in enumerate(v0.alt)
        if trim(v0.ref, b) == trim(v.ref, a)), -1) for a in v.alt] or None
    m1 = v1 and [0] + [next((i+1 for i,b in enumerate(v1.alt)
        if trim(v1.ref, b) == trim(v.ref, a)), -1) for a in v.alt] or None
    m2 = v2 and [0] + [next((i+1 for i,b in enumerate(v2.alt)
        if trim(v2.ref, b) == trim(v.ref, a)), -1) for a in v.alt] or None

    for k,d in iteritems(f0.infos):
        t = d and d['Number'] or None
        if t == 'A':
            if k in ('AC','AF','AN'):
                continue
        elif t == 'R':
            if k == 'RPA':
                d0 = d1 = d2 = [-1] * (len(v.alt)+1)
                if v0 and v0.info.get(k):
                    s = v0.info.get(k)
                    d0 = [i >= 0 and s[i] or -1 for i in m0]
                if v1 and v1.info.get(k):
                    s = v1.info.get(k)
                    d1 = [i >= 0 and s[i] or -1 for i in m1]
                if v2 and v2.info.get(k):
                    s = v2.info.get(k)
                    d2 = [i >= 0 and s[i] or -1 for i in m2]
                d = [max(e) for e in zip(d0,d1,d2)]
                if min(d) > 0:
                    v.info[k] = d
                continue
        else:
            if k in ('STR', 'RU'):
                d = None
                d = d and d or v1 and v1.info.get(k)
                d = d and d or v2 and v2.info.get(k)
                d = d and d or v0 and v0.info.get(k)
                if d:
                    v.info[k] = d
                continue
            if k == 'DP':
                d = -1
                if d < 0 and v0:
                    d = v0.info.get(k, -1)
                if d < 0:
                    d = 0
                    d += v1 and v1.info.get(k, 0) or 0
                    d += v2 and v2.info.get(k, 0) or 0
                if d >= 0:
                    v.info[k] = d
                continue
            if k == 'ML_PROB':
                d = (u.info.get(k) for u in (v1,v2) if u)
                d = [e for e in d if e is not None]
                if d:
                    d = functools.reduce(operator.mul, d, 1)
                    v.info[k] = round(d, 2)
                continue
            if k not in ('QD',):
                continue
        v.info.pop(k, None)

    for k,d in iteritems(f0.formats):
        t = d and d['Number'] or None
        if t == 'R':
            if k == 'AD':
                d0 = d1 = d2 = None
                if v0 and v0.samples[0].get(k):
                    s = v0.samples[0].get(k)
                    d0 = [i >= 0 and s[i] or 0 for i in m0]
                if v1 and v1.samples[0].get(k):
                    s = v1.samples[0].get(k)
                    d1 = [i >= 0 and s[i] or 0 for i in m1]
                if v2 and v2.samples[0].get(k):
                    s = v2.samples[0].get(k)
                    d2 = [i >= 0 and s[i] or 0 for i in m2]
                if d1 and d2:
                    d = [sum(e) for e in zip(d1,d2)]
                    d0 = None
                else:
                    d = d1 or d2 or d0
                if d and d0 and d != d0:
                    d[0] = max(d[0], sum(d0) - sum(d[1:]))
                if d:
                    v.samples[0][k] = d
                continue
        if t == 'G':
            if k == 'PL':
                infty = 999999999
                d1 = d2 = [0] + [infty] * len(v.alt)
                if v1 and v1.samples[0].get(k):
                    s = v1.samples[0].get(k)
                    d1 = [s[i] if i >= 0 else infty for i in m1]
                if v2 and v2.samples[0].get(k):
                    s = v2.samples[0].get(k)
                    d2 = [s[i] if i >= 0 else infty for i in m2]
                gt = [(i,j) for j in range(len(v.alt)+1) for i in range(j+1)]
                pl = [min(d1[i]+d2[j],d1[j]+d2[i]) for i,j in gt]
                if v0 and v0.samples[0].get(k):
                    s = v0.samples[0].get(k)
                    d0 = [infty if i < 0 or j < 0 else s[j*(j+1)/2+i]
                        for i,j in (sorted((m0[i],m0[j])) for i,j in gt)]
                    pl = [e if e < infty else f for e,f in zip(pl,d0)]
                mx = sorted(e for e in pl if e < infty)
                gq = mx[1] - mx[0]
                mn, mx = mx[0], max(mx[-1] - mx[0] + gq, gq * 10)
                pl = [min(e - mn, mx) for e in pl]
                v.samples[0][k] = pl
                v.samples[0]['GQ'] = gq
                continue
        else:
            if k == 'DP':
                d = -1
                if d < 0 and v0:
                    d = v0.samples[0].get(k, -1)
                if d < 0:
                    d = 0
                    d += v1 and v1.samples[0].get(k, 0) or 0
                    d += v2 and v2.samples[0].get(k, 0) or 0
                if d >= 0:
                    v.samples[0][k] = d
                continue

    v.end = v.pos + len(v.ref)
    v.line = None
    return v

def merge2(vcfi1, vcfi2, vcfi3, vcfi0, vcfo, bed=None, **kwargs):
    for pos, grp, ovl in grouper(vcfi1, vcfi2, vcfi3, vcfi0):
        v1, v2, v3, v0 = grp

        if v0 and (v0.filter or v0.samples[0].get('GT') in (None, '0/0')):
            v = None
        else:
            v = v0
        if v and any(len(v.ref) == len(a) for a in v.alt):
            v1 = v2 = v3 = None

        if v1 and v1.info.get('DELTA') or v2 and v2.info.get('DELTA'):
            v = join2(vcfi0, vcfi1, vcfi2, v0, v1, v2, pos, bed)
        elif v3 and v3.info.get('DELTA'):
            if v3.filter or v3.samples[0].get('GT') in (None, '0/0'):
                v = None
            else:
                v = v3
                v.info.pop('DELTA', None)
                v.line = None
        elif v:
            ps = (bed and v and v.samples[0].get('PS') and
                  bed.get(v.chrom, pos, pos) or None)
            if ps:
                v.samples[0]['PS'] = ps[0][0]+1
                v.line = None

        if v:
            pgt = v.samples[0].get('PGT')
            pid = v.samples[0].get('PID')
            if pgt or pid:
                v.samples[0].pop('PGT', None)
                v.samples[0].pop('PID', None)
                v.line = None
            vcfo.emit(v)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-t",
        "--thread_count",
        type=int,
        default=multiprocessing.cpu_count(),
        help="Number of threads",
    )
    subparsers = parser.add_subparsers()

    merge_parser = subparsers.add_parser("merge", help="Merge multiple VCFs")
    merge_parser.add_argument("--hap1", required=True, help="The hap1 VCF")
    merge_parser.add_argument("--hap2", required=True, help="The hap2 VCF")
    merge_parser.add_argument(
        "--unphased", required=True, help="The unphased VCF"
    )
    merge_parser.add_argument(
        "--phased", required=True, help="The phased SNV VCF"
    )
    merge_parser.add_argument(
        "--bed", required=True, help="The BED file of phased regions"
    )
    merge_parser.add_argument("output_vcf", help="The merged output VCF")
    merge_parser.set_defaults(func=merge_main)

    haploid_patch_parser = subparsers.add_parser(
        "haploid_patch", help="Patch haploid DNAscope and DNAscopeHP VCF files"
    )
    haploid_patch_parser.add_argument(
        "--patch1", required=True, help="The output patch1 VCF"
    )
    haploid_patch_parser.add_argument(
        "--patch2", required=True, help="The output patch2 VCF"
    )
    haploid_patch_parser.add_argument(
        "--hap1_hp", required=True, help="The hap1 DNAscopeHP VCF"
    )
    haploid_patch_parser.add_argument(
        "--hap2_hp", required=True, help="The hap2 DNAscopeHP VCF"
    )
    haploid_patch_parser.add_argument("--phased", help="The phased VCF")
    haploid_patch_parser.add_argument("--hap1", help="The hap1 VCF")
    haploid_patch_parser.add_argument("--hap2", help="The hap2 VCF")
    haploid_patch_parser.set_defaults(func=patch1_main)

    patch_parser = subparsers.add_parser(
        "patch", help="Patch DNAscope and DNAscopeHP VCF files"
    )
    patch_parser.add_argument("--vcf", required=True, help="The DNAscope VCF")
    patch_parser.add_argument(
        "--vcf_hp", required=True, help="The DNAscopeHP VCF"
    )
    patch_parser.add_argument("output_vcf", help="The patched output VCF")
    patch_parser.set_defaults(func=patch2_main)

    haploid_patch2_parser = subparsers.add_parser(
        "haploid_patch2",
        help="Patch a single pair of haploid DNAscope and DNAscopeHP VCFs",
    )
    haploid_patch2_parser.add_argument(
        "--vcf", required=True, help="The DNAscope vcf"
    )
    haploid_patch2_parser.add_argument(
        "--vcf_hp", required=True, help="The DNAscopeHP vcf"
    )
    haploid_patch2_parser.add_argument(
        "--patch_vcf", required=True, help="The output patch vcf"
    )
    haploid_patch2_parser.set_defaults(func=haploid_patch_main)

    return (parser.parse_args(argv), parser, subparsers)


def main(args, parser, subparsers):
    if not hasattr(args, "func"):
        subcommands = list(subparsers.choices.keys())
        msg = (
            "{}: error: a subcommand is required. Possible subcommands are: {}"
        )
        print(msg.format(parser.prog, ", ".join(subcommands)), file=sys.stderr)
        sys.exit(1)
    return args.func(args)


if __name__ == "__main__":
    args, parser, subparsers = parse_args()
    sys.exit(main(args, parser, subparsers))

# vim: ts=4 sw=4 expandtab
