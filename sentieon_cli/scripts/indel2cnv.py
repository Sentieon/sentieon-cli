import collections
import Levenshtein
import math
import io
import vcflib
from scipy.signal import find_peaks
from multiprocessing import Pool
import argparse
import logging
import time
import resource
import sys
import platform
import mappy
import os
import re

# Convert long INDELs (from assembly-based SV calls) into tandem CNV truth events.
# Pipeline: find_repeats (detect period) -> match_ref (locate repeat on reference)
#   -> extend array boundaries -> emit DEL/DUP with REFSTART/REFSTOP for cnv_eval.
# Handles segmental duplications: insertion may be distant from (POSSHIFT2 > 0) or
# a truncated copy of (SVLEN < PERIOD) an existing tandem array — both are valid CNVs
# because reads from the new copy map back to the reference array, increasing depth.

MIN_SEQ_LEN = 400
REF_SEARCH_RNG = 50000
REF_SEARCH_CNT = 6
MAX_EXTEND_ITERS = 200
MAX_EXTEND_SEQ = 20000
MAX_SEQ_LEN = 200000
MIN_PERIOD_RATIO = 0.1  # skip if matched portion < 10% of period (coincidental similarity)

log = logging.getLogger(__name__)

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
        s = max(0, s)
        e = min(e, ci.length)
        if s >= e:
            return ''
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

def copy_variant(v, **overrides):
    args = [getattr(v, s) for s in v.__slots__]
    new_v = vcflib.Variant(*args)
    new_v.samples = [dict(s) for s in v.samples]
    new_v.alt = list(v.alt)
    for k, val in overrides.items():
        setattr(new_v, k, val)
    return new_v

# Detect repeat period(s) in an INDEL sequence via autocorrelation-like scanning.
# Compares var_seq[0:n] vs var_seq[n:2n] at increasing n, finds peaks, refines with
# binary search, and explores harmonic multipliers (2x, 3x, ...) as candidate periods.
def find_repeats(var_seq, thresh):
    if len(var_seq) < MIN_SEQ_LEN:
        return []
    length = len(var_seq)
    step = 200 if length > 2000 else 100
    start = step * 2
    scan_limit = min(length//2+step*5, length-step*1)
    dist = [Levenshtein.ratio(var_seq[:n], var_seq[n:n*2]) for n in range(start, scan_limit, step)]
    if not dist or max(dist) < thresh:
        return [len(var_seq)]
    peaks, _ = find_peaks(dist, distance=start//step, height=0.95)
    if peaks.size == 0:
        mid = length//2
        if Levenshtein.ratio(var_seq[:mid], var_seq[mid:mid*2]) > thresh:
            peaks = [mid/step]
    else:
        peaks = [p+start//step for p in peaks]
    if len(peaks) > 0:
        est_period = min([peaks[i+1] - peaks[i] for i in range(len(peaks)-1)] + [peaks[0]])
        pidx = [int(round(p/est_period)) for p in peaks]
        period = int((sum([pidx[i] * p for i, p in enumerate(peaks)])/sum([pidx[i]**2 for i, p in enumerate(peaks)]) if len(peaks) > 1 else peaks[0]) * step)
        cnt = int(round(length/period))
        if cnt < 2:
            return [len(var_seq)]
        p_min, p_max = int(length/(cnt+0.5)), int(length/(cnt-0.5))
        p_step = max(1, (p_max - p_min)//10)
        while True:
            ratios = [max([Levenshtein.ratio(var_seq[i*p:(i+1)*p], var_seq[(i+1)*p:(i+2)*p]) for i in range(cnt-1)]) for p in range(p_min, p_max+1, p_step)]
            best_ratio = max(ratios)
            best_period = p_min + (ratios.index(best_ratio))*p_step
            if p_step == 1:
                break
            p_min, p_max = best_period - 2 * p_step, best_period + 2 * p_step
            p_step = max(1, (p_max - p_min)//10)
        ratio_orig = max([Levenshtein.ratio(var_seq[i*period:(i+1)*period], var_seq[(i+1)*period:(i+2)*period]) for i in range(cnt-1)])
        if best_ratio > ratio_orig:
            period = best_period
        periods = {period: max(best_ratio, ratio_orig)}
        i = 2
        while True:
            cur_period = period * i
            if cur_period * 2 > length:
                break
            cur_cnt = int(round(length/cur_period))
            if cur_cnt >= 2:
                periods[cur_period] = max([Levenshtein.ratio(var_seq[j*cur_period:(j+1)*cur_period], var_seq[(j+1)*cur_period:(j+2)*cur_period]) for j in range(cur_cnt-1)])
            i += 1
        periods = {k: v for k, v in periods.items() if k >= MIN_SEQ_LEN}
        if len(periods) > 1:
            min_ratio = max(thresh, max(periods.values()) - 0.03)
            periods = {k: v for k, v in periods.items() if v > min_ratio}
        return list(periods.keys())
    return [len(var_seq)]

# Align one repeat unit (seq) to a reference window around the INDEL position using mappy.
# If only part of seq aligned (partial match), extend the alignment boundary via binary
# search on Levenshtein ratio to handle split period boundaries.
# Then walk outward in steps of matched_ref_length to find the full tandem repeat array.
# Returns: (match_start, match_stop), (array_start, array_stop), (seq_start, seq_stop), ratio
def mm2_match(seq, ref, chrom, pos, half_win, thresh):
    ref_search_start = max(0, pos - half_win)
    ref_search_stop = pos + half_win
    max_nm_ratio = 1 - thresh
    min_match_length = MIN_SEQ_LEN
    ref_seq = ref.get(chrom, ref_search_start, ref_search_stop)
    if not ref_seq:
        return None
    aligner = mappy.Aligner(seq=ref_seq, preset='map-ont')
    hits = list(aligner.map(seq))
    if not hits:
        return None
    # filter hits
    filtered = []
    for h in hits:
        if h.strand != 1 or h.is_primary == False:
            continue
        nm = h.NM if hasattr(h, 'NM') else (h.blen - h.mlen)
        nm_rate = nm / h.blen if h.blen > 0 else 1.0
        if nm_rate >= max_nm_ratio or h.blen <= min_match_length:
            continue
        filtered.append((h, nm_rate))
    if not filtered:
        return None
    min_nm_rate = min(f[1] for f in filtered)
    filtered = [(h, nr) for h, nr in filtered if nr <= min(min_nm_rate + 0.01, max_nm_ratio)]
    best_hit, best_nm_rate = max(filtered, key=lambda x: x[0].blen)
    seq_start, seq_stop = best_hit.q_st, best_hit.q_en
    ref_start = best_hit.r_st + ref_search_start
    ref_stop = best_hit.r_en + ref_search_start
    matched_seq_length = seq_stop - seq_start
    matched_ref_length = ref_stop - ref_start
    best_matched_ref = ref.get(chrom, ref_start, ref_stop)
    matched_seq = seq[seq_start:seq_stop]
    best_ratio = Levenshtein.ratio(best_matched_ref, matched_seq)
    if best_ratio < thresh:
        return None
    # Partial match extension: when mappy aligned only part of seq, the period boundary
    # may be split (e.g. seq = [tail|head] of repeat). Extend the match into adjacent
    # reference to recover the full period, using binary search on Levenshtein ratio.
    if abs(len(seq) - matched_seq_length) > 50 and (len(seq)/matched_seq_length-1) < thresh and len(seq) <= MAX_EXTEND_SEQ:
        if len(seq) - seq_stop < 10:
            ref_after = ref.get(chrom, ref_stop, ref_stop + len(seq) - matched_seq_length)
            step = min(100, len(ref_after)//5)
            step = max(1, step)
            start_i = step
            end_i = len(ref_after)
            seq_len_diff = len(seq) - seq_stop
            last_valid_ratio = -1
            last_valid = 0
            last_invalid = -1
            best_ext = 0
            while True:
                for l in range(start_i, end_i, step):
                    ratio = Levenshtein.ratio(best_matched_ref + ref_after[:l], seq[seq_start:]+seq[:l-seq_len_diff])
                    if ratio > thresh:
                        last_valid = l
                        last_valid_ratio = ratio
                    else:
                        last_invalid = l
                        break
                if last_invalid == -1 or step == 1:
                    best_ext = last_valid
                    break
                step //= 2
                start_i = last_valid
                end_i = last_invalid
            if last_valid_ratio > 0:
                best_ratio = last_valid_ratio
                ref_stop += best_ext
                best_matched_ref += ref_after[:best_ext]
                matched_seq = seq[seq_start:] + seq[:best_ext - seq_len_diff]
                matched_seq_length = len(matched_seq)
                seq_stop = seq_start + matched_seq_length
                matched_ref_length += best_ext
        elif seq_start < 10:
            ref_before = ref.get(chrom, ref_start - len(seq) + matched_seq_length, ref_start)
            step = min(100, len(ref_before)//5)
            step = max(1, step)
            start_i = step
            end_i = len(ref_before)
            last_valid_ratio = -1
            last_valid = 0
            last_invalid = -1
            best_ext = 0
            while True:
                for l in range(start_i, end_i, step):
                    ratio = Levenshtein.ratio(ref_before[-l:] + best_matched_ref, seq[-l+seq_start:] + seq[:seq_stop])
                    if ratio > thresh:
                        last_valid = l
                        last_valid_ratio = ratio
                    else:
                        last_invalid = l
                        break
                if last_invalid == -1 or step == 1:
                    best_ext = last_valid
                    break
                step //= 2
                start_i = last_valid
                end_i = last_invalid
            if last_valid_ratio > 0:
                best_ratio = last_valid_ratio
                ref_start -= best_ext
                best_matched_ref = ref_before[-best_ext:] + best_matched_ref
                matched_seq = seq[-best_ext+seq_start:] + seq[:seq_stop]
                matched_seq_length = len(matched_seq)
                seq_start = seq_stop - len(matched_seq)
                matched_ref_length += best_ext
    # Walk outward from the match to find the full tandem repeat array on reference.
    # REFSTART/REFSTOP define where read-depth change is expected (used by cnv_eval).
    ref_start_ext = ref_start
    ref_stop_ext = ref_stop
    for _ in range(MAX_EXTEND_ITERS):
        ref_start_ext -= matched_ref_length
        if ref_start_ext < 0 or Levenshtein.ratio(ref.get(chrom, ref_start_ext, ref_start_ext+matched_ref_length), best_matched_ref) < thresh:
            break
    for _ in range(MAX_EXTEND_ITERS):
        ref_stop_ext += matched_ref_length
        if Levenshtein.ratio(ref.get(chrom, ref_stop_ext-matched_ref_length, ref_stop_ext), best_matched_ref) < thresh:
            break
    return (ref_start, ref_stop), (ref_start_ext, ref_stop_ext), (seq_start, seq_stop), best_ratio

# Try to match one repeat unit in the reference immediately adjacent to the INDEL.
# For insertions: align alt_seq[:period] against ref_after[:p] + ref_before[-(period-p):]
# at varying split positions p (binary search). Handles wrap-around period boundaries.
# For deletions: the deleted ref sequence itself is the repeat unit.
def local_match(v, ref, period, thresh):
    if len(v.ref) < len(v.alt[0]):
        alt_seq = v.alt[0]
        ref_before = ref.get(v.chrom, v.pos-period, v.pos)
        ref_after = ref.get(v.chrom, v.pos, v.pos+period)
        step = max(period//40, 100)
        cnt = int(round(len(alt_seq)/period))
        offset = len(alt_seq) - cnt * period
        if offset > 0:
            ref_after = alt_seq[-offset:] + ref_after
        elif offset < 0:
            ref_after = ref_after[-offset:]
        alt_period_seq = alt_seq[:period]
        ratios = {p: Levenshtein.ratio(ref_after[:p] + ref_before[-(period-p):], alt_period_seq) for p in range(0, period, step)}
        ratios[period] = Levenshtein.ratio(ref_after, alt_period_seq)
        while True:
            best_pos = max(ratios, key=ratios.get)
            max_ratio = ratios[best_pos]
            if max_ratio < 0.8:
                return None
            if step == 1:
                break
            step = step//2
            p_lo = max(0, best_pos - step)
            p_hi = min(best_pos + step, period)
            for p in (p_lo, p_hi):
                if p not in ratios:
                    ratios[p] = Levenshtein.ratio(ref_after[:p] + ref_before[-(period-p):], alt_period_seq)
        if max_ratio < thresh:
            return None
        ref_start = v.pos - (period - best_pos)
        ref_end = v.pos + best_pos
        ref_seq = ref_before[-(period-best_pos):] + alt_seq[:best_pos] if ref_start != v.pos else alt_period_seq
    else:
        ref_start = v.pos
        ref_end = v.pos + period
        ref_seq = v.ref[:period]
        max_ratio = 1.
    ref_start_ext = ref_start
    ref_end_ext = ref_end
    for _ in range(MAX_EXTEND_ITERS):
        ref_start_ext -= period
        if ref_start_ext < 0 or Levenshtein.ratio(ref.get(v.chrom, ref_start_ext, ref_start_ext+period), ref_seq) < thresh:
            break
    for _ in range(MAX_EXTEND_ITERS):
        ref_end_ext += period
        if Levenshtein.ratio(ref.get(v.chrom, ref_end_ext-period, ref_end_ext), ref_seq) < thresh:
            break
    return (ref_start, ref_end), (ref_start_ext, ref_end_ext), (0, period), max_ratio

# Try local_match first (fast, exact adjacent match). Fall back to mm2_match
# (mappy alignment in a wider window) for segdups where the repeat may be distant.
def match_ref(v, ref, period, thresh):
    local_result = local_match(v, ref, period, thresh)
    if local_result is not None:
        return local_result
    return mm2_match(v.alt[0][:period], ref, v.chrom, v.pos, min(max(period, len(v.alt[0]), len(v.ref))*REF_SEARCH_CNT, REF_SEARCH_RNG), thresh)

# Process one input INDEL variant into CNV call(s).
# Multi-allelic variants are split; nearby het INS with opposite phase are merged.
# For each allele: find_repeats -> match_ref -> pick best period -> emit DUP/DEL.
def proc_variant(v, ref, thresh):
    stime = time.time()
    result_vs = []
    # Split multi-allelic into per-allele variants, trim shared prefix/suffix
    split_var = len(v.alt) > 1
    vs = []
    is_insert = []
    if split_var:
        for i, a in enumerate(v.alt):
            if a == '*':
                continue
            if len(a) < len(v.ref):
                if len(v.ref) > MIN_SEQ_LEN and len(v.ref) - len(a) > MIN_SEQ_LEN:
                    v1 = copy_variant(v, alt=[a])
                    if len(a) > 1:
                        if a == v1.ref[:len(a)] or len(a) < 10 or Levenshtein.ratio(a, v1.ref[:len(a)]) > thresh:
                            v1.pos += len(a) - 1
                            v1.ref = v.ref[len(a):]
                            v1.alt = [a[-1]]
                        elif Levenshtein.ratio(a, v.ref[-len(a):]) > thresh:
                            v1.ref = v.ref[:-len(a)]
                            v1.alt = [v.ref[0]]
                    v1.samples[0]['GT'] = '0/1'
                    vs.append(v1)
                    is_insert.append(False)
            elif len(a) - len(v.ref) > MIN_SEQ_LEN:
                v1 = copy_variant(v)
                v1.samples[0]['GT'] = '0/1'
                if len(v.ref) < 6 or Levenshtein.ratio(v.alt[i][:len(v.ref)], v.ref) > min(thresh, 1-5/len(v.ref)):
                    v1.alt = [v1.alt[i][len(v.ref)-1:]]
                    v1.pos += len(v.ref) - 1
                    v1.ref = v.ref[-1]
                elif len(v.ref) < 6 or Levenshtein.ratio(v.alt[i][-len(v.ref):], v.ref) > min(thresh, 1-5/len(v.ref)):
                    v1.alt = [v1.alt[i][:-len(v.ref)]]
                    v1.ref = v.ref[0]
                else:
                    v1.alt = [v1.alt[i]]
                vs.append(v1)
                is_insert.append(True)
    else:
        vs = [v]
        is_insert = [len(v.alt[0]) > len(v.ref)]
    periods = []
    for si, vv in enumerate(vs):
        if not is_insert[si]:
            del_len = len(vv.ref) - len(vv.alt[0])
            if del_len >= MIN_SEQ_LEN:
                if del_len > MAX_SEQ_LEN:
                    gt = re.split(r"\||/", vv.samples[0]['GT'].replace('.', '0'))
                    gt_sum = sum(int(g) for g in gt)
                    end_pos = vv.pos + len(vv.ref)
                    vv.info = {'SVTYPE': 'DEL', 'CNDIFF': -gt_sum, 'POSSHIFT': 0, 'POSSHIFT2': 0,
                               'DIST': 1.0, 'SVLEN': del_len, 'PERIOD': del_len,
                               'REFSTART': vv.pos + 1, 'REFSTOP': end_pos, 'END': end_pos,
                               'INDEL': '%s:%d' % (vv.chrom, vv.pos + 1)}
                    vv.end = end_pos - 1
                    vv.ref = ref.get(vv.chrom, vv.pos, vv.pos + 1)
                    vv.alt = ['<DEL>']
                    vv.line = None
                    result_vs.append(vv)
                    periods.append([])
                else:
                    periods.append(find_repeats(vv.ref[1:], thresh))
            else:
                periods.append([])
        else:
            if len(vv.alt[0]) - len(vv.ref) >= MIN_SEQ_LEN and len(vv.alt[0]) <= MAX_SEQ_LEN:
                periods.append(find_repeats(vv.alt[0], thresh))
            else:
                periods.append([])
    if len(sum(periods, [])) == 0:
        return [], time.time() - stime
    # If both alleles are same direction with matching periods, merge as homozygous
    if (all(is_insert) or not any(is_insert)) and len(periods) > 1 and min(len(p) for p in periods) > 0:
        updated_periods = []
        for p0 in periods[0]:
            p1s = [p1 for p1 in periods[1] if abs(p0/p1-1) <= 0.05]
            if p1s and Levenshtein.ratio(*(a[:min(p0, p1s[0])] for a in v.alt)) > 1 - (1-thresh)/2:
                updated_periods.append(min(p0, p1s[0]))
        if updated_periods:
            periods = [updated_periods]
            vs = [vs[1]] if len(vs[1].alt[0]) > len(vs[0].alt[0]) else [vs[0]]
            vs[0].samples[0]['GT'] = '1/1'

    for si, sv in enumerate(vs):
        period = periods[si]
        if not period:
            continue
        results = []
        for j, p in enumerate(period):
            try:
                result = match_ref(sv, ref, p, thresh)
            except Exception as e:
                log.warning("match_ref failed for %s:%d period=%d: %s", sv.chrom, sv.pos, p, e)
                result = None
            if result:
                results.append((j, *result))
        if not results:
            continue
        # Pick the best period: highest Levenshtein ratio, break ties by matched length
        if len(results) == 1:
            result = results[0]
        else:
            best_ratios = [r[4] for r in results]
            best_ratio = max(best_ratios)
            best_i = best_ratios.index(best_ratio)
            result = results[best_i]
            best_seq_length = result[3][1] - result[3][0]
            for ri, r in enumerate(results):
                if r and best_i != ri and best_ratio - r[4] < 0.01:
                    seq_length = r[3][1] - r[3][0]
                    if seq_length > best_seq_length:
                        result = r
                        best_i = ri
                        best_seq_length = seq_length

        ri, best_match, ref_range, seq_range, best_ratio = result
        best_period = period[ri]
        if is_insert[si]:
            period_cnt = len(sv.alt[0]) / best_period
        else:
            overlap = max(0, min(sv.pos + len(sv.ref), ref_range[1]) - max(sv.pos, ref_range[0]))
            period_cnt = overlap / best_period
        period_cnt = int(round(period_cnt))
        seq_length = seq_range[1] - seq_range[0]
        # Filter: matched portion must be large enough and a meaningful fraction of period
        if seq_length < MIN_SEQ_LEN or seq_length / best_period < MIN_PERIOD_RATIO:
            continue
        gt = re.split(r"\||/", sv.samples[0]['GT'].replace('.', '0'))
        gt_sum = sum(int(g) for g in gt)
        is_dup = is_insert[si]
        # POSSHIFT: distance from INDEL to matched repeat on reference
        # POSSHIFT2: distance from INDEL to nearest edge of tandem array (0 = inside array)
        # REFSTART/REFSTOP: tandem array boundaries (where depth change is expected)
        indel_pos = '%s:%d' % (sv.chrom, sv.pos + 1)
        sv.info = {'SVTYPE': 'DUP' if is_dup else 'DEL',
                'CNDIFF': period_cnt * (gt_sum if is_dup else -gt_sum),
                'POSSHIFT': best_match[0] - sv.pos,
                'POSSHIFT2': sv.pos - ref_range[1] if sv.pos > ref_range[1] else (ref_range[0] - sv.pos if ref_range[0] > sv.pos else 0),
                'DIST': best_ratio,
                'SVLEN': seq_length,
                'PERIOD': best_period,
                'REFSTART': ref_range[0] + 1,
                'REFSTOP': ref_range[1] + 1,
                'END': best_match[1] + 1,
                'INDEL': indel_pos}
        sv.pos = best_match[0]
        sv.end = best_match[1]
        sv.ref = ref.get(sv.chrom, sv.pos, sv.pos+1)
        sv.alt = ['<DUP>'] if is_dup else ['<DEL>']
        sv.line = None
        result_vs.append(sv)
    return result_vs, time.time() - stime

def proc_batch(ref, variants, thresh):
    results = []
    for v in variants:
        try:
            results.append(proc_variant(v, ref, thresh))
        except Exception as e:
            log.error("Failed processing %s:%d: %s", v.chrom, v.pos, e)
            results.append(([], 0.0))
    return results

# Merge two adjacent het INS variants with opposite phase (1|0 + 0|1) into one hom variant.
def merge(ref, v1, v2):
    if v1.chrom != v2.chrom or v2.pos - v1.pos > 10:
        return None
    if set([v1.samples[0]['GT'], v2.samples[0]['GT']]) != set(['1|0', '0|1']):
        return None
    if len(v1.alt) != len(v2.alt) or len(v2.alt) != 1:
        return None
    if len(v1.ref) > len(v1.alt[0]) or len(v1.ref) > len(v2.alt[0]):
        return None
    ref_gap = ref.get(v1.chrom, v1.end, v2.pos+1)
    v1.ref += ref_gap
    v1.alt[0] += ref_gap
    v1.alt.append(v1.ref[:-1] + v2.alt[0])
    v1.samples[0]['GT'] = '1|1'
    return v1

def gt_set(gt):
    return [int(g) for g in re.split(r"\||/", gt) if g != '.']

# Adjust CNDIFF/SVLEN so that the repeat count fits the array size.
# E.g. CNDIFF=2 with SVLEN=2*period in a 4-copy array -> CNDIFF=1 with SVLEN=4*period.
def extend_v(v):
    cndiff = v.info['CNDIFF']
    period = v.info['PERIOD']
    if not cndiff or not period:
        return v
    num_periods = round(v.info['SVLEN'] / period)
    gt_sum = sum(gt_set(v.samples[0]['GT'])) or 1
    diff = cndiff // gt_sum
    max_num_period = (v.info['REFSTOP'] - v.info['REFSTART']) / period
    if max_num_period - round(max_num_period) < 0.1:
        max_num_period = round(max_num_period)
    else:
        max_num_period = math.ceil(max_num_period) - 1
    best_d = -1
    for d in range(abs(diff), 0, -1):
        sp = num_periods * abs(diff) / d
        if abs(sp - round(sp)) > 0.1 or sp > max_num_period:
            continue
        best_d = d
    if abs(diff) == best_d or best_d == -1:
        return v
    svlen_periods = round(num_periods * abs(diff) / best_d)
    v.info['CNDIFF'] = best_d if cndiff > 0 else -best_d
    v.info['SVLEN'] = svlen_periods * period
    v.pos = v.info['REFSTART']
    v.end = v.pos + v.info['SVLEN']
    v.line = None
    return v

# Merge overlapping output CNV calls from the same tandem array.
# Handles: (1) opposite-phase same-type pairs -> hom, (2) opposite-type pairs -> net CNDIFF,
# (3) same-type DEL/DUP pairs with compatible periods -> combined CNDIFF.
def merge_output_vars(variants):
    if not variants:
        return []
    result = []
    last = variants[0]
    for v in variants[1:]:
        if last is None:
            last = v
            continue
        merged = _try_merge_output(last, v)
        if merged == 'drop':
            last = None
        elif merged:
            last = merged
        else:
            result.append(last)
            last = v
    if last is not None:
        result.append(last)
    return result

def _try_merge_output(v1, v2):
    if v1.chrom != v2.chrom:
        return None
    rs1, re1 = v1.info['REFSTART'], v1.info['REFSTOP']
    rs2, re2 = v2.info['REFSTART'], v2.info['REFSTOP']
    overlap = min(re1, re2) - max(rs1, rs2)
    if overlap <= 0:
        return None
    len1, len2 = re1 - rs1, re2 - rs2
    if len1 <= 0 or len2 <= 0:
        return None
    r1, r2 = overlap / len1, overlap / len2
    if max(r1, r2) < 0.9:
        return None
    # Opposite types (DUP+DEL or DEL+DUP): net the CNDIFF
    if v1.alt[0] != v2.alt[0] and min(r1, r2) > 0.9:
        net = v1.info['CNDIFF'] + v2.info['CNDIFF']
        if net == 0:
            return 'drop'  # cancel out
        keep = v1 if abs(v1.info['CNDIFF']) > abs(v2.info['CNDIFF']) else v2
        keep.info['CNDIFF'] = net
        keep.line = None
        return keep
    # Same type DEL
    if v1.alt[0] == '<DEL>' and v2.alt[0] == '<DEL>':
        keep = v1 if len1 >= len2 else v2
        if keep.info['CNDIFF'] == -2 or v1.samples[0]['GT'] == v2.samples[0]['GT']:
            return keep
        if min(r1, r2) > 0.9:
            gt1, gt2 = gt_set(v1.samples[0]['GT']), gt_set(v2.samples[0]['GT'])
            dot = '.' in v1.samples[0]['GT'] or '.' in v2.samples[0]['GT']
            if not dot:
                gt = [gt1[i] | gt2[i] for i in range(min(len(gt1), len(gt2)))]
                keep.samples[0]['GT'] = '|'.join(str(g) for g in gt)
            else:
                keep.samples[0]['GT'] = str(gt1[0] | gt2[0]) + '/.'
            keep.info['CNDIFF'] = min(-2, v1.info['CNDIFF'] + v2.info['CNDIFF'])
            keep.line = None
            return keep
    # Same type DUP with compatible periods
    if v1.alt[0] == '<DUP>' and v2.alt[0] == '<DUP>':
        p1, p2 = v1.info['PERIOD'], v2.info['PERIOD']
        r = max(p1, p2) / min(p1, p2) if min(p1, p2) > 0 else 999
        if abs(r - round(r)) > 0.1:
            log.warning("Incompatible periods %d vs %d at %s:%d", p1, p2, v1.chrom, v1.pos)
            return None
        period = min(p1, p2)
        keep, other = (v1, v2) if len1 >= len2 else (v2, v1)
        # Normalize CNDIFF to common period; use round to avoid truncating partial matches
        n1 = max(1, round(keep.info['SVLEN'] / period))
        n2 = max(1, round(other.info['SVLEN'] / period))
        cd1 = keep.info['CNDIFF'] * n1
        cd2 = other.info['CNDIFF'] * n2
        gt1, gt2 = gt_set(v1.samples[0]['GT']), gt_set(v2.samples[0]['GT'])
        if gt1 == gt2 or (cd1 == cd2 and [1, 1] not in (gt1, gt2)):
            dot = '.' in v1.samples[0]['GT'] or '.' in v2.samples[0]['GT']
            if not dot:
                gt = [gt1[i] | gt2[i] for i in range(min(len(gt1), len(gt2)))]
                keep.samples[0]['GT'] = '|'.join(str(g) for g in gt)
            else:
                keep.samples[0]['GT'] = str(gt1[0] | gt2[0]) + '/.'
            keep.info['PERIOD'] = period
            keep.info['SVLEN'] = period
            keep.info['CNDIFF'] = cd1 + cd2
            keep = extend_v(keep)
            keep.line = None
            return keep
    return None

def main(args):
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=log_level, format='%(asctime)s %(levelname)s %(message)s',
                        datefmt='%H:%M:%S')
    ref = Reference(args.ref)
    input_vcf = vcflib.VCF(args.input_vcf, 'r')
    fout = vcflib.VCF(args.out_vcf, 'w')
    update = ('##INFO=<ID=POSSHIFT,Number=1,Type=Integer,Description="Best matching CNV sequence position shifted from INDEL position">',
              '##INFO=<ID=POSSHIFT2,Number=1,Type=Integer,Description="Distance from original INDEL position to CNV sequence reference boundary">',
              '##INFO=<ID=CNDIFF,Number=1,Type=Integer,Description="Copy number state difference">',
              '##INFO=<ID=DIST,Number=1,Type=Float,Description="Levenshtein distance">',
              '##INFO=<ID=PERIOD,Number=1,Type=Integer,Description="Length of the repeating unit">',
              '##INFO=<ID=REFSTART,Number=1,Type=Integer,Description="Reference start locus of the tandem repeat array">',
              '##INFO=<ID=REFSTOP,Number=1,Type=Integer,Description="Reference stop locus of the tandem repeat array">',
              '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SVLEN">',
              '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SVTYPE">',
              '##ALT=<ID=DEL,Description="Deletion relative to the reference">',
              '##ALT=<ID=DUP,Description="Duplication relative to the reference">',
              '##INFO=<ID=END,Number=1,Type=Integer,Description="END location">',
              '##INFO=<ID=INDEL,Number=1,Type=String,Description="Source INDEL locus (chrom:pos)">',
            )
    fout.copy_header(input_vcf, update=update)
    fout.emit_header()
    all_contigs = [c for c, t in input_vcf.contigs.items()]
    if args.contig:
        input_contigs = args.contig.split(',')
        contigs = [c for c in all_contigs if c in input_contigs]
        if not contigs:
            log.error("Contig %s not found in the VCF.", args.contig)
            sys.exit(1)
    else:
        contigs = all_contigs
    t0 = time.time()
    variants = []
    for c in contigs:
        vcf_c = input_vcf.range(c)
        for v in vcf_c:
            if (len(v.filter) == 0 or 'PASS' in v.filter) and max(len(a) for a in [v.ref] + v.alt) >= MIN_SEQ_LEN:
                if variants:
                    v_new = merge(ref, variants[-1], v)
                    if v_new:
                        variants[-1] = v_new
                        continue
                variants.append(v)
    log.info("Loaded %d qualifying variants from %d contigs", len(variants), len(contigs))
    if args.threads == 1:
        all_vars = []
        for v in variants:
            vs, t = proc_variant(v, ref, args.thresh)
            all_vars += vs
            if vs:
                for sv in vs:
                    log.info('%s:%d %s SVLEN=%d shift=%d,%d (%.3fs)',
                             sv.chrom, sv.pos, sv.alt[0], sv.info['SVLEN'],
                             sv.info['POSSHIFT'], sv.info['POSSHIFT2'], t)
            else:
                log.debug('%s:%d No CNV (%.3fs)', v.chrom, v.pos, t)
    else:
        n_threads = args.threads or os.cpu_count()
        # Sort by expected difficulty (sequence length) descending so that
        # expensive variants are dispatched first and don't pile up at the end.
        # Use small batches (5) for better load balancing across workers.
        variants.sort(key=lambda v: max(len(a) for a in [v.ref] + v.alt), reverse=True)
        batch_size = 5
        args_in = []
        for i in range(0, len(variants), batch_size):
            args_in.append((ref, variants[i:i+batch_size], args.thresh))
        all_vars = []
        with Pool(n_threads) as p:
            for result in p.starmap(proc_batch, args_in):
                for vs, t in result:
                    all_vars += vs
    # Merge overlapping output calls (same array, different haplotypes) and write
    for c in contigs:
        out_c = sorted([v for v in all_vars if v.chrom == c], key=lambda v: v.pos)
        out_c = merge_output_vars(out_c)
        for v in sorted(out_c, key=lambda v: (v.pos, v.end)):
            fout.emit(v)
    t1 = time.time()
    mm, ut, st = 0, 0, 0
    mem_scale = 1 if platform.system() == 'Darwin' else 1024
    for who in (resource.RUSAGE_SELF, resource.RUSAGE_CHILDREN):
        ru = resource.getrusage(who)
        mm += ru.ru_maxrss * mem_scale
        ut += ru.ru_utime
        st += ru.ru_stime
    log.info('overall: %d mem %.3f user %.3f sys %.3f real', mm, ut, st, t1-t0)
    fout.close()
    input_vcf.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert long INDELs to CNV calls")
    parser.add_argument('ref', help='Reference FASTA')
    parser.add_argument('input_vcf', help='Input VCF with long INDELs')
    parser.add_argument('out_vcf', help='Output VCF file name')
    parser.add_argument('--contig', help="Contigs to process (comma-separated)")
    parser.add_argument('-t', '--threads', help="Concurrent processes (default: nproc)", default=None, type=int)
    parser.add_argument('--thresh', help="Minimum Levenshtein similarity ratio", default=0.9, type=float)
    parser.add_argument('-v', '--verbose', action='store_true', help="Verbose logging")
    args = parser.parse_args()
    main(args)
