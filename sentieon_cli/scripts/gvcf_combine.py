from __future__ import print_function
import heapq
import sys
import argparse
import vcflib
import copy
import time
import resource
import collections
import io

from vcflib.compat import *

Contig = collections.namedtuple('Contig', 'length offset width skip')

class Reference(object):
    def __init__(self, path):
        self.path = path
        self.index = collections.OrderedDict()
        with io.open(self.path+'.fai', 'rb') as fp:
            for line in fp:
                flds = line.rstrip().decode().split('\t')
                if len(flds) != 5:
                    raise RuntimeError('Corrupted index')
                self.index[flds[0]] = Contig._make(map(int,flds[1:]))
        self.fp = io.open(path, 'rb')

    def __getstate__(self):
        odict = self.__dict__.copy()
        odict['fp'] = self.fp.tell()
        return odict

    def __setstate__(self, ndict):
        path = ndict['path']
        fp = io.open(path, 'rb')
        fp.seek(ndict.pop('fp'))
        ndict['fp'] = fp
        self.__dict__.update(ndict)

    def get(self, c, s, e):
        ci = self.index[c]
        e = min(e, ci.length)
        if s > e:
            raise ValueError
        seq = b''
        while s < e:
            o = ci.offset + s // ci.width * ci.skip + s % ci.width
            n = ci.width - s % ci.width
            n = min(n, e-s)
            self.fp.seek(o)
            seq += self.fp.read(n)
            s += n
        return seq.decode()

    def __iter__(self):
        return iter(self.index.items())

class Combiner(vcflib.Shardable, vcflib.ShardResult):

    info_vals = {"AC": 0, "MLEAC": 0, "AF": 0, "MLEAF": 0, "RPA": 0}
    sample_vals = {"AD": lambda s: 0 if 'DP' not in s else s['DP'] - sum(s['AD'])}

    def __init__(self, ref, gvcf_input, vcf_input, gvcf_output):
        self.ref = Reference(ref)
        self.vcf_in = vcflib.VCF(vcf_input, 'r')
        self.gvcf_in = vcflib.VCF(gvcf_input, 'r')
        self.contigs = self.gvcf_in.contigs
        extras = self.extra_headers(self.gvcf_in, self.vcf_in)
        self.vcftmp = None
        self.gvcf_out = vcflib.VCF(gvcf_output, 'w')
        self.gvcf_out.copy_header(self.gvcf_in, extras)
        self.gvcf_out.emit_header()
        self.remove_info_keys = [k for k, v in self.vcf_in.infos.items() if v['Number'] in ('A', 'R') and k not in self.info_vals]
        self.remove_sample_keys = [k for k, v in self.vcf_in.formats.items() if v['Number'] in ('A', 'R') and k not in self.sample_vals]
    
    @staticmethod
    def extra_headers(vcf1, vcf2):
        filter_diff = ['##FILTER=<ID=%s,' %f for f in vcf2.filters.keys() if f not in vcf1.filters]
        info_diff = ['##INFO=<ID=%s,' %f for f in vcf2.infos.keys() if f not in vcf1.infos]
        fmt_diff = ['##FORMAT=<ID=%s,' %f for f in vcf2.formats.keys() if f not in vcf1.formats]
        diffs = filter_diff + info_diff + fmt_diff
        return [h for h in vcf2.headers if any([h.startswith(s) for s in diffs])]

    def __del__(self):
        self.vcf_in.close()
        self.gvcf_in.close()
        if self.gvcf_out != '-':
            self.gvcf_out.close()

    def grouper(self, gvcfi, vcfi):
        q = []
        iters = [iter(i) for i in (gvcfi, vcfi)]
        for k, i in enumerate(iters):
            v = next(i, None)
            if v:
                heapq.heappush(q, (v.pos, v.end, k, v, i))
        while q:
            grp = [None, None]
            pos, end, k, v, i = heapq.heappop(q)
            grp[k] = v
            vv = next(i, None)
            if vv:
                heapq.heappush(q, (vv.pos, vv.end, k, vv, i))
            if q and q[0][2] != k:
                pos2, end2, k2, v2, i2 = heapq.heappop(q)
                if k < k2:
                    g, v = v, v2
                else:
                    g = v2
                if len(g.alt) == 1:
                    if Combiner.ovl(g, v):
                        for gg in self.split_g(g, v):
                            heapq.heappush(q, (gg.pos, gg.end, 0, gg, iters[0]))
                        heapq.heappush(q, (v.pos, v.end, 1, v, iters[1]))
                        continue
                    heapq.heappush(q, (pos2, end2, k2, v2, i2))
                elif pos2 == pos:
                    if g.end > v.end:
                        g.ref = g.ref[v.end - v.pos]
                        g.pos = v.end
                        g.alt = ['<NON_REF>']
                        vlist = []
                        while q:
                            vv = heapq.heappop(q)
                            vlist.append(vv)
                            if vv[2] == 0:
                                g.end = vv[3].pos
                                break
                        for vv in vlist:
                            heapq.heappush(q, vv)
                        if not vlist or vv[2] != 0:
                            vv = next(iters[0], None)
                            if vv:
                                g.end = vv.pos
                                heapq.heappush(q, (vv.pos, vv.end, 0, vv, iters[0]))
                        g.info = {'END': g.end}
                        g.qual = None
                        g.samples[0] = {'GT': '0/0', 'DP': v2.samples[0].get('DP'), 'GQ': 15, 'PL': [0, 15, 99]}
                        g.line = None
                        heapq.heappush(q, (g.pos, g.end, 0, g, iters[0]))
                    else:
                        grp[k2] = v2
                        nv = next(i2, None)
                        if nv:
                            heapq.heappush(q, (nv.pos, nv.end, k2, nv, i2))
                else:
                    heapq.heappush(q, (pos2, end2, k2, v2, i2))
            yield (pos, grp)

    @staticmethod
    def ovl(g, v):
        if not len(g.alt) == 1:
            return False
        return v.pos >= g.pos and v.pos < g.end

    def split_g(self, g, v):
        gs = []
        if g.pos < v.pos:
            g1 = copy.deepcopy(g)
            g1.end = v.pos
            g1.info['END'] = v.pos
            g1.line = None
            gs.append(g1)
        if g.end > v.end:
            g.ref = self.ref.get(g.chrom, v.end, v.end+1)
            g.pos = v.end
            g.line = None
            gs.append(g)
        return gs
    
    def __shard__(self, cse):
        self.shard = cse
        self.vcftmp = self.gvcf_out.__shard__(cse)
        return self    

    def __getstate__(self):
        odict = self.__dict__.copy()
        return odict
        
    def __getdata__(self):
        return self.vcftmp.__getdata__()
    
    def __accum__(self, data):
        self.gvcf_out.__accum__(data)

    def combine(self, shard=None):
        shard = shard or self.shard
        chrom, start, end = shard
        vcfs = [vcf.__shard__(shard) for vcf in (self.gvcf_in, self.vcf_in)]
        for pos, grp in self.grouper(*vcfs):
            g, v = grp
            if v is None:
                if len(g.alt) == 1:
                    if pos < start:
                        g.pos = start
                        g.line = None
                    if g.end > end:
                        g.info['END'] = end
                        g.line = None
                    self.vcftmp.emit(g)
                elif pos >= start:
                    g.samples[0]['GT'] = '0/0'
                    g.line = None
                    self.vcftmp.emit(g)
            if v and pos >= start:
                v.alt.append('<NON_REF>')
                if "PL" in v.samples[0]:
                    n_alt = len(v.alt)
                    idx = -1
                    for i in range(n_alt+1):
                        idx += n_alt + 1 - i
                        v.samples[0]['PL'].insert(idx, 100)
                for k in self.remove_info_keys:
                    if k in v.info:
                        del v.info[k]
                for k in self.remove_sample_keys:
                    if k in v.samples[0]:
                        del v.samples[0][k]
                for k, u in self.info_vals.items():
                    if k in v.info:
                        v.info[k].append(u if not callable(u) else u(v.info))
                for k, u in self.sample_vals.items():
                    if k in v.samples[0]:
                        v.samples[0][k].append(u if not callable(u) else u(v.samples[0]))
                v.line = None
                self.vcftmp.emit(v)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ref", help="Reference file")
    parser.add_argument("gvcf", help="Input GVCF file")
    parser.add_argument("vcf", help="Input VCF file")
    parser.add_argument("gvcf_out", help="Output GVCF file")
    parser.add_argument("-t", "--threads", help="Number of threads", type=int)
    
    args = parser.parse_args()
    combiner = Combiner(args.ref, args.gvcf, args.vcf, args.gvcf_out)

    t0 = time.time()
    nthr = args.threads
    step = 10*1000*1000

    sharder = vcflib.Sharder(nthr)
    contigs = ((c,0,int(t['length'])) for c,t in iteritems(combiner.contigs))
    shards = sharder.cut(contigs, step)
    _ = sharder.run(shards, Combiner.combine, [], combiner)

    t1 = time.time()
    mm, ut, st = 0, 0, 0
    for who in (resource.RUSAGE_SELF, resource.RUSAGE_CHILDREN):
        ru = resource.getrusage(who)
        mm += ru.ru_maxrss * 1024
        ut += ru.ru_utime
        st += ru.ru_stime
    print('overall: %d mem %.3f user %.3f sys %.3f real' %
        (mm, ut, st, t1-t0), file=sys.stderr)

    return 0

if __name__ == '__main__':
    sys.exit(main())

# vim: ts=4 sw=4 expandtab
