import collections
import Levenshtein
import math
import io
import vcflib
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
import datetime

# Convert long INDELs (from assembly-based SV calls) into tandem CNV truth events.
# Pipeline: find_repeats (detect period) -> match_ref (locate repeat on reference)
#   -> extend array boundaries -> emit DEL/DUP with REFSTART/REFSTOP for cnv_eval.
# Handles segmental duplications: insertion may be distant from (POSSHIFT2 > 0) or
# a truncated copy of (SVLEN < PERIOD) an existing tandem array — both are valid CNVs
# because reads from the new copy map back to the reference array, increasing depth.

INDEL2CNV_VERSION = '2.1.0'  # combine_sv_cnv release this file belongs to

MIN_SEQ_LEN = 400
REF_SEARCH_RNG = 50000
REF_SEARCH_CNT = 6
MAX_EXTEND_ITERS = 200
MAX_EXTEND_SEQ = 20000
MAX_SEQ_LEN = 200000
MIN_PERIOD_RATIO = 0.1  # skip if matched portion < 10% of period (coincidental similarity)
MAX_N_FRAC = 0.5        # gate 16: a record at least half reference N is a gap artefact, not a CNV

def n_frac(seq):
    """Fraction of N in a sequence (0 for an empty one)."""
    return seq.upper().count('N') / len(seq) if seq else 0.0

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

# ---------------------------------------------------------------------------
# Repeat structure of an indel sequence
# ---------------------------------------------------------------------------
# A tandem expansion or contraction is k copies (k need not be an integer) of a
# unit of length P. Earlier versions scored fixed-offset chunks on a coarse grid
# and could neither land on an arbitrary P nor express a fractional k. Here the
# period candidates come from a k-mer gap histogram (single-base resolution) and
# are verified by aligning the sequence against itself shifted by the candidate:
# the alignment absorbs variable-length copies and partial copies.
KMER = 10               # k-mer size for the gap histogram
MAX_KMER_BINS = 24      # strongest gap bins (by support density) that get verified by alignment
MIN_OVERLAP = 50        # a partial tail counts as evidence from this length ...
MIN_OVERLAP_FRAC = 0.3  # ... and this fraction of the period (>= 1.3 copies). A weak fundamental cannot hurt the
                        # array walk (it tests the unit-step hypothesis too); it still becomes UNIT, the rung base
                        # when >= MIN_PERIOD, and the fit-ladder base, so it is a hypothesis, not a fact.
MIN_CORE_COVER = 0.5    # k-mer support must span at least half the sequence
SHIFT_WINDOW = 20000    # verify long sequences on a prefix of the core (>= 3 periods)
REFINE_WINDOW = 0.15    # search +-15% around a k-mer gap for the best period (a bin can sit on the shoulder of a broad peak)
LADDER_SLACK = 0.03     # rungs within this of the best self-identity are kept
ARRAY_PERIOD_MARGIN = 0.02  # an array-level period must beat a microsatellite fundamental by this (capped at 0.99)
MAX_RUNG = 16
MIN_PERIOD = MIN_SEQ_LEN  # reporting floor: a smaller unit is reported as one block

RepeatInfo = collections.namedtuple('RepeatInfo', 'period ratio core_start core_end candidates array_period')

_BASE_CODE = {ord('A'): 0, ord('C'): 1, ord('G'): 2, ord('T'): 3,
              ord('a'): 0, ord('c'): 1, ord('g'): 2, ord('t'): 3}


def _kmer_gap_bins(seq, k=KMER):
    """Gaps between consecutive occurrences of each k-mer, merged within 2%.

    Returns [(center, support, first, last)] sorted by support, where
    [first, last+k) is the span of positions that support the gap (the core).
    """
    n = len(seq)
    mask = (1 << (2 * k)) - 1
    last_pos = {}
    gaps = collections.Counter()
    first, last = {}, {}
    h = valid = 0
    for i, ch in enumerate(seq.encode()):
        c = _BASE_CODE.get(ch)
        if c is None:
            h = valid = 0
            continue
        h = ((h << 2) | c) & mask
        valid += 1
        if valid < k:
            continue
        start = i - k + 1
        j = last_pos.get(h)
        if j is not None:
            d = start - j
            gaps[d] += 1
            if d not in first:
                first[d] = j
            last[d] = start
        last_pos[h] = start
    min_support = max(3, n // 500)
    bins = []
    for d, c in sorted(gaps.items(), key=lambda t: -t[1]):
        if c < min_support:
            break
        tol = max(2, int(0.02 * d))
        for b in bins:
            if abs(b[0] - d) <= tol:
                b[1] += c
                b[2] = min(b[2], first[d])
                b[3] = max(b[3], last[d])
                break
        else:
            if len(bins) < 200:
                bins.append([d, c, first[d], last[d]])
    bins.sort(key=lambda b: -b[1] / max(1, n - b[0]))   # support per possible position
    return [tuple(b) for b in bins]


MAX_PAIRS = 40          # adjacent-copy pairs scored per period (evenly subsampled beyond this)


def _shift_ratio(core, d, cutoff=0.0):
    """Self-similarity of core at period d: length-weighted mean Levenshtein ratio of
    adjacent d-length chunks, plus the partial tail against the start of the previous
    chunk when the tail is long enough to mean anything.

    Chunks are distinct substrings, so the alignment cannot cheat by sliding one copy
    onto the next (aligning core[:L-d] with core[d:] scores 1-d/(L-d) for any small d).
    The core is the k-mer-supported span, so junk outside it does not enter the score.
    """
    L = len(core)
    if d < 1:
        return 0.0
    k, t = divmod(L, d)
    if k < 1:
        return 0.0
    pairs = [(i * d, (i + 1) * d, d) for i in range(k - 1)]
    if len(pairs) > MAX_PAIRS:
        step = len(pairs) / MAX_PAIRS
        pairs = [pairs[int(j * step)] for j in range(MAX_PAIRS)]
    if t >= max(MIN_OVERLAP, MIN_OVERLAP_FRAC * d):
        pairs.append(((k - 1) * d, k * d, t))
    if not pairs:
        return 0.0
    num = den = 0
    for a, b, w in pairs:
        num += w * Levenshtein.ratio(core[a:a + w], core[b:b + w], score_cutoff=cutoff)
        den += w
    return num / den


def _refine_period(core, d0, cutoff):
    """Local maximum of the shift ratio around d0, coarse-to-fine, re-centred if it
    lands on the window edge."""
    best_d, best_r = d0, _shift_ratio(core, d0, cutoff)
    w0 = max(2, int(REFINE_WINDOW * d0))
    res = max(1, d0 // 1000)      # final resolution: 1 bp up to 1 kb periods
    for _ in range(4):
        w, step = w0, max(res, w0 // 4)
        edge_lo, edge_hi = best_d - w, best_d + w
        while True:
            for d in range(best_d - w, best_d + w + 1, step):
                if d != best_d:
                    r = _shift_ratio(core, d, cutoff)
                    if r > best_r:
                        best_d, best_r = d, r
            if step <= res:
                break
            w, step = step, max(res, step // 4)
        if edge_lo < best_d < edge_hi:
            break
    return best_d, best_r


def _window(seq, lo, hi, d):
    """The core span, cut to a prefix of SHIFT_WINDOW (or 3 periods) for long sequences."""
    c = seq[lo:hi]
    lim = max(SHIFT_WINDOW, 3 * d)
    return c[:lim] if len(c) > lim else c


def analyze_repeats(seq, thresh):
    """Fundamental period of an indel sequence and the candidate units to try on
    the reference. The whole sequence (one block) is always a candidate."""
    n = len(seq)
    if n < MIN_SEQ_LEN:
        return RepeatInfo(None, 0.0, 0, n, [], None)
    fundamental, f_ratio, core = None, 0.0, (0, n)
    bins = sorted(_kmer_gap_bins(seq)[:MAX_KMER_BINS])
    for d0, support, lo, hi in bins:
        v = _verify_bin(seq, bins, d0, lo, hi + KMER, thresh)
        if v:
            fundamental, f_ratio, core = v[0], v[1], (v[2], v[3])
            break
    # A microsatellite fundamental hides the array-level structure (e.g. 837 bp units that are
    # mostly a 9-mer). Look for the verified period >= MIN_PERIOD whose self-identity is clearly
    # higher than the fundamental's; proc_variant tries it when the block does not fit the array.
    array_period = None
    if fundamental is not None and fundamental < MIN_PERIOD:
        # must beat the fundamental by the margin; above 0.97 the margin shrinks and at 0.99 a tie is enough,
        # since a near-perfect microsatellite leaves no room above it
        best = min(f_ratio + ARRAY_PERIOD_MARGIN, 0.99)
        for d0, support, lo, hi in bins:
            if d0 < MIN_PERIOD:
                continue
            v = _verify_bin(seq, bins, d0, lo, hi + KMER, thresh)
            if v and v[1] >= best:
                array_period, best = v[0], v[1]
    return RepeatInfo(fundamental, f_ratio, core[0], core[1],
                      _unit_candidates(seq, fundamental, f_ratio, core, thresh), array_period)


def _verify_bin(seq, bins, d0, lo, hi, thresh):
    """Refine a k-mer gap bin to the best nearby period and verify it on its core.
    Returns (period, ratio, core_lo, core_hi) or None."""
    n = len(seq)
    if hi - lo < MIN_CORE_COVER * n:
        return None
    c = _window(seq, lo, hi, d0)
    if _shift_ratio(c, d0, thresh - 0.15) < thresh - 0.1:
        return None
    d, r = _refine_period(c, d0, thresh - 0.15)
    # the refined period may belong to a neighbouring bin (e.g. a bin sitting on the
    # shoulder of the true peak): take the union of the supporting spans and rescore
    for e0, _, elo, ehi in bins:
        if e0 != d0 and abs(e0 - d) <= max(2, 0.02 * d):
            lo, hi = min(lo, elo), max(hi, ehi + KMER)
            r = _shift_ratio(_window(seq, lo, hi, d), d, thresh - 0.15)
            break
    return (d, r, lo, hi) if r >= thresh else None


def _unit_candidates(seq, period, ratio, core, thresh):
    """Units to match against the reference: the smallest rung of the period ladder
    within LADDER_SLACK of the best self-identity (a higher-order repeat wins when its
    copies are clearly more similar), the best rung, and the whole sequence as one
    block. A sub-MIN_PERIOD fundamental is reported as the block only."""
    n = len(seq)
    cands = []
    if period is not None and period >= MIN_PERIOD:
        c = seq[core[0]:core[1]]
        if len(c) > max(SHIFT_WINDOW, 3 * MAX_RUNG * period):
            c = c[:max(SHIFT_WINDOW, 3 * MAX_RUNG * period)]
        rungs = [(period, ratio)]
        for i in range(2, MAX_RUNG + 1):
            p = i * period
            if len(c) - p < p:
                break
            r = _shift_ratio(c, p, thresh - 0.05)
            if r >= thresh:
                rungs.append((p, r))
        best = max(r for p, r in rungs)
        keep = [(p, r) for p, r in rungs if r >= best - LADDER_SLACK]
        cands.append(keep[0][0])
        best_p = max(keep, key=lambda t: t[1])[0]
        if best_p != cands[0]:
            cands.append(best_p)
    if n not in cands:
        cands.append(n)
    return cands


def find_repeats(var_seq, thresh):
    """Candidate repeat units of an indel sequence (see analyze_repeats)."""
    return analyze_repeats(var_seq, thresh).candidates


# ---------------------------------------------------------------------------
# Tandem array extension on the reference
# ---------------------------------------------------------------------------
EXTEND_SLACK = 0.0      # array walk continues down to thresh - EXTEND_SLACK (set by --extend_slack)
MAX_EXTEND_STEP = 20000 # cap on the chunk compared per step when the fundamental is known (whole periods)
PARTIAL_WINDOW = 100    # local window that must still match at the end of a partial step
MIN_PARTIAL = 50


def _partial_step(ref_chunk, unit_ext, thresh):
    """Longest prefix of ref_chunk that still matches the tiled unit, trimmed so the
    last PARTIAL_WINDOW bases match on their own."""
    hi = min(len(ref_chunk), len(unit_ext))
    lo = 0
    while hi - lo > 8:
        mid = (lo + hi) // 2
        if Levenshtein.ratio(ref_chunk[:mid], unit_ext[:mid], score_cutoff=thresh) >= thresh:
            lo = mid
        else:
            hi = mid
    l, w = lo, PARTIAL_WINDOW
    while l >= w and Levenshtein.ratio(ref_chunk[l - w:l], unit_ext[l - w:l],
                                       score_cutoff=thresh - 0.05) < thresh - 0.05:
        l -= w // 2
    return l if l >= MIN_PARTIAL else 0


def _step_hypotheses(unit, unit_period):
    """Ways the reference array may continue past the matched unit: (chunk, E_right, E_left).
    The unit itself (step L) is always one hypothesis. When the indel has a fundamental shorter
    than the unit, its tiling in whole periods is another: it keeps the phase when the unit is a
    fractional number of copies, and it is cheaper for very long units. The walk picks per side
    by the first step, so a fundamental that does not describe the reference array is harmless."""
    L = len(unit)
    hyps = [(L, unit, unit)]
    if unit_period and 0 < unit_period < L:
        C = unit_period * max(1, min(L // unit_period, MAX_EXTEND_STEP // unit_period))
        if C != L:
            hyps.append((C, unit[L - C:], unit[:C]))
    return hyps


EXTEND_THRESH = None   # set by --match_thresh: the array walk keeps the main threshold

def _extend_array(ref, chrom, start, stop, unit_period, thresh):
    if EXTEND_THRESH is not None:
        thresh = EXTEND_THRESH
    """Walk outward from the matched reference unit [start, stop) while the reference keeps
    tiling it, then a partial step finds the ragged end of the array."""
    unit = ref.get(chrom, start, stop)
    L = len(unit)
    if L == 0:
        return start, stop
    thresh = thresh - EXTEND_SLACK
    hyps = _step_hypotheses(unit, unit_period)
    # right
    best = None
    for C, E, _ in hyps:
        r = Levenshtein.ratio(ref.get(chrom, stop, stop + C), E, score_cutoff=thresh - 0.1)
        if best is None or r > best[0]:
            best = (r, C, E)
    _, C, E = best
    pos = stop
    for _ in range(MAX_EXTEND_ITERS):
        chunk = ref.get(chrom, pos, pos + C)
        if len(chunk) < C or Levenshtein.ratio(chunk, E, score_cutoff=thresh) < thresh:
            break
        pos += C
    pos += _partial_step(ref.get(chrom, pos, pos + C), E, thresh)
    # left
    best = None
    for C, _, E in hyps:
        r = Levenshtein.ratio(ref.get(chrom, max(0, start - C), start), E, score_cutoff=thresh - 0.1)
        if best is None or r > best[0]:
            best = (r, C, E)
    _, C, E = best
    neg = start
    for _ in range(MAX_EXTEND_ITERS):
        nxt = neg - C
        if nxt < 0:
            break
        chunk = ref.get(chrom, nxt, neg)
        if len(chunk) < C or Levenshtein.ratio(chunk, E, score_cutoff=thresh) < thresh:
            break
        neg = nxt
    lo = max(0, neg - C)
    neg -= _partial_step(ref.get(chrom, lo, neg)[::-1], E[::-1], thresh)
    return neg, pos


# Align one repeat unit (seq) to a reference window around the INDEL position using mappy.
# If only part of seq aligned (partial match), extend the alignment boundary via binary
# search on Levenshtein ratio to handle split period boundaries.
# Then walk outward in steps of matched_ref_length to find the full tandem repeat array.
# Returns: (match_start, match_stop), (array_start, array_stop), (seq_start, seq_stop), ratio
_COMP = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')

def _revcomp(s):
    return s.translate(_COMP)[::-1]


def mm2_match(seq, ref, chrom, pos, half_win, thresh, unit_period=None):
    ref_search_start = max(0, pos - half_win)
    ref_search_stop = pos + half_win
    inverted = False
    # One identity scale for both matching paths. thresh is a Levenshtein ratio, 1 - edits/(2 * length),
    # so a ratio of thresh allows an edit rate of 2 * (1 - thresh); mappy's NM/blen is that edit rate.
    # (v1 used 1 - thresh here, which made this path twice as strict as the local one.)
    max_nm_ratio = 2 * (1 - thresh)
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
        if h.is_primary == False:
            continue
        nm = h.NM if hasattr(h, 'NM') else (h.blen - h.mlen)
        nm_rate = nm / h.blen if h.blen > 0 else 1.0
        if nm_rate >= max_nm_ratio or h.blen <= min_match_length:
            continue
        filtered.append((h, nm_rate))
    if not filtered:
        return None
    # an inverted copy (reverse strand) changes read depth like a direct one; use it only when
    # no forward-strand hit passes, so existing behaviour is unchanged elsewhere
    if any(h.strand == 1 for h, _ in filtered):
        filtered = [(h, nr) for h, nr in filtered if h.strand == 1]
    else:
        seq = _revcomp(seq)
        hits = list(aligner.map(seq))
        filtered = [(h, (h.NM if hasattr(h, 'NM') else h.blen - h.mlen) / h.blen) for h in hits
                    if h.is_primary and h.strand == 1 and h.blen > min_match_length]
        filtered = [(h, nr) for h, nr in filtered if nr < max_nm_ratio]
        if not filtered:
            return None
        inverted = True
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
    if not inverted and abs(len(seq) - matched_seq_length) > 50 and (len(seq)/matched_seq_length-1) < thresh and len(seq) <= MAX_EXTEND_SEQ:
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
                for l in range(max(start_i, seq_len_diff), end_i, step):
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
            if last_valid_ratio > 0 and best_ext > 0:
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
            if last_valid_ratio > 0 and best_ext > 0:
                best_ratio = last_valid_ratio
                ref_start -= best_ext
                best_matched_ref = ref_before[-best_ext:] + best_matched_ref
                matched_seq = seq[-best_ext+seq_start:] + seq[:seq_stop]
                matched_seq_length = len(matched_seq)
                seq_start = seq_stop - len(matched_seq)
                matched_ref_length += best_ext
    # Walk outward from the match to find the full tandem repeat array on reference.
    # REFSTART/REFSTOP define where read-depth change is expected (used by cnv_eval).
    ref_start_ext, ref_stop_ext = _extend_array(ref, chrom, ref_start, ref_stop, unit_period, thresh)
    if inverted:   # query coordinates back in the insertion's own orientation
        seq_start, seq_stop = len(seq) - seq_stop, len(seq) - seq_start
    mm2_match.last_inverted = inverted
    return (ref_start, ref_stop), (ref_start_ext, ref_stop_ext), (seq_start, seq_stop), best_ratio

# Try to match one repeat unit in the reference immediately adjacent to the INDEL.
# For insertions: align alt_seq[:period] against ref_after[:p] + ref_before[-(period-p):]
# at varying split positions p (binary search). Handles wrap-around period boundaries.
# For deletions: the deleted ref sequence itself is the repeat unit.
def local_match(v, ref, period, thresh, unit_period=None):
    if len(v.ref) < len(v.alt[0]):
        alt_seq = v.alt[0]
        ref_before = ref.get(v.chrom, v.pos-period, v.pos)
        ref_after = ref.get(v.chrom, v.pos, v.pos+period)
        step = max(period//40, 100)
        # v1 prepended the ALT's partial last copy to ref_after when the copy count was fractional.
        # That scored ALT bases against themselves (DIST 0.94 where the reference alone gives 0.78)
        # and reported a window shifted by the partial copy; only reference bases are compared now.
        alt_period_seq = alt_seq[:period]
        ratios = {p: Levenshtein.ratio(ref_after[:p] + ref_before[-(period-p):], alt_period_seq, score_cutoff=0.8) for p in range(0, period, step)}
        ratios[period] = Levenshtein.ratio(ref_after, alt_period_seq, score_cutoff=0.8)
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
                    ratios[p] = Levenshtein.ratio(ref_after[:p] + ref_before[-(period-p):], alt_period_seq, score_cutoff=0.8)
        if max_ratio < thresh:
            return None
        ref_start = v.pos - (period - best_pos)
        ref_end = v.pos + best_pos
    else:
        ref_start = v.pos + 1            # first deleted base (v.pos is the anchor)
        ref_end = v.pos + 1 + period
        max_ratio = 1.
    # the array walk compares the reference itself over [ref_start, ref_end) (v1 used a hybrid of
    # reference and ALT bases for insertions; they agree wherever the inserted copy matches)
    ref_start_ext, ref_end_ext = _extend_array(ref, v.chrom, ref_start, ref_end, unit_period, thresh)
    return (ref_start, ref_end), (ref_start_ext, ref_end_ext), (0, period), max_ratio

# Try local_match first (fast, exact adjacent match). Fall back to mm2_match
# (mappy alignment in a wider window) for segdups where the repeat may be distant.
def match_ref(v, ref, period, thresh, unit_period=None):
    mm2_match.last_inverted = False
    local_result = local_match(v, ref, period, thresh, unit_period)
    if local_result is not None:
        return local_result
    return mm2_match(v.alt[0][:period], ref, v.chrom, v.pos, min(max(period, len(v.alt[0]), len(v.ref))*REF_SEARCH_CNT, REF_SEARCH_RNG), thresh, unit_period)

RATIO_SLACK = 0.02   # match ratios within this of the best are treated as equal
NET_RECIPROCAL = 0.8 # a DUP and a DEL whose arrays overlap reciprocally by this much are one array and net
EXACT_FRAC = 0.05    # |round(n/p)*p - n| within this fraction of n counts as a whole copy count

def select_candidate(results, periods, indel_len):
    """Pick one matched unit. Among matches within RATIO_SLACK of the best ratio prefer the
    longest reference array, then a whole copy count (the block is exact by construction),
    then the smallest unit. results: (idx, match, array, seq_range, ratio)."""
    best_ratio = max(r[4] for r in results)
    pool = [r for r in results if best_ratio - r[4] <= RATIO_SLACK]
    max_ext = max(r[2][1] - r[2][0] for r in pool)
    ext_tol = 0.5 * min(periods[r[0]] for r in pool)
    def key(r):
        p = periods[r[0]]
        ext = r[2][1] - r[2][0]
        err = abs(int(indel_len / p + 0.5) * p - indel_len)
        return (max_ext - ext > ext_tol, err > EXACT_FRAC * indel_len, p, -r[4])
    return min(pool, key=key)

def _indel_len(v, is_insert):
    return len(v.alt[0]) if is_insert else len(v.ref) - len(v.alt[0])


def _fit_unit(sv, ref, unit, n, thresh, max_len=None):
    """An insertion longer than the reference array cannot match as one block. Find the
    largest multiple of the fundamental unit (below n, at least MIN_SEQ_LEN) that matches
    the reference next to the insertion. The match is roughly monotone in the unit length
    (it fails once the unit outgrows the array), but local_match's split search is a
    heuristic, so probe a coarse ladder first and binary search from the largest success.
    Returns (period, match) or None."""
    top = n - 1 if not max_len else min(n - 1, max_len)   # a unit cannot be longer than the array
    lo, hi = -(-MIN_SEQ_LEN // unit), top // unit
    if lo > hi:
        return None
    probes = sorted({lo + (hi - lo) * i // 6 for i in range(7)}, reverse=True)
    best = None
    fail = hi + 1
    for k in probes:
        r = local_match(sv, ref, k * unit, thresh, unit)
        if r is not None:
            best, lo = (k * unit, r), k
            break
        fail = k
    if best is None:
        return None
    while fail - lo > 1:
        mid = (lo + fail) // 2
        r = local_match(sv, ref, mid * unit, thresh, unit)
        if r is not None:
            lo, best = mid, (mid * unit, r)
        else:
            fail = mid
    return best


# Process one input INDEL variant into CNV call(s).
# Multi-allelic variants are split; nearby het INS with opposite phase are merged.
# For each allele: find_repeats -> match_ref -> pick best period -> emit DUP/DEL.
def proc_variant(v, ref, thresh, match_thresh=None):
    # thresh: period verification (and, unless --match_thresh is given, matching and the array walk);
    # match_thresh: the identity a local or mappy match must reach to be accepted
    mt = thresh if match_thresh is None else match_thresh
    stime = time.time()
    result_vs = []
    # Split multi-allelic into per-allele variants, trim shared prefix/suffix
    split_var = len(v.alt) > 1
    vs = []
    is_insert = []
    if split_var:
        # Only process alleles referenced by GT
        gt_indices = set(int(g) for g in re.split(r'[|/]', v.samples[0].get('GT', '0/0')) if g not in ('.', '0'))
        for i, a in enumerate(v.alt):
            if (i + 1) not in gt_indices:
                continue
            if a == '*':
                continue
            if len(a) < len(v.ref):
                if len(v.ref) > MIN_SEQ_LEN and len(v.ref) - len(a) > MIN_SEQ_LEN:
                    v1 = copy_variant(v, alt=[a])
                    if len(a) > 1:
                        if a == v1.ref[:len(a)] or len(a) < 10 or Levenshtein.ratio(a, v1.ref[:len(a)]) > thresh:
                            v1.pos += len(a) - 1
                            v1.ref = v.ref[len(a) - 1:]       # the new anchor is the allele's last base
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
    units = []
    aperiods = []
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
                               'REFSTART': vv.pos + 2, 'REFSTOP': end_pos, 'END': end_pos,
                               'INDEL': '%s:%d' % (vv.chrom, vv.pos + 1)}
                    vv.end = end_pos
                    vv.ref = ref.get(vv.chrom, vv.pos, vv.pos + 1)
                    vv.alt = ['<DEL>']
                    vv.line = None
                    result_vs.append(vv)
                    periods.append([])
                    units.append(None)
                    aperiods.append(None)
                else:
                    info = analyze_repeats(vv.ref[1:], thresh)
                    periods.append(info.candidates)
                    units.append(info.period)
                    aperiods.append(info.array_period)
            else:
                periods.append([])
                units.append(None)
                aperiods.append(None)
        else:
            if len(vv.alt[0]) - len(vv.ref) >= MIN_SEQ_LEN and len(vv.alt[0]) <= MAX_SEQ_LEN:
                info = analyze_repeats(vv.alt[0], thresh)
                periods.append(info.candidates)
                units.append(info.period)
                aperiods.append(info.array_period)
            else:
                periods.append([])
                units.append(None)
                aperiods.append(None)
    if len(sum(periods, [])) == 0:
        return result_vs, time.time() - stime
    # If exactly 2 alleles, same direction with matching periods, merge as homozygous
    if (all(is_insert) or not any(is_insert)) and len(periods) == 2 and min(len(p) for p in periods) > 0 \
            and abs(_indel_len(vs[0], is_insert[0]) / _indel_len(vs[1], is_insert[1]) - 1) <= 0.05:
        updated_periods = []
        for p0 in periods[0]:
            p1s = [p1 for p1 in periods[1] if abs(p0/p1-1) <= 0.05]
            if p1s and Levenshtein.ratio(vs[0].alt[0][:min(p0, p1s[0])], vs[1].alt[0][:min(p0, p1s[0])]) > 1 - (1-thresh)/2:
                updated_periods.append(min(p0, p1s[0]))
        if updated_periods:
            periods = [updated_periods]
            keep = 1 if len(vs[1].alt[0]) > len(vs[0].alt[0]) else 0
            vs = [vs[keep]]
            units = [units[keep] or units[1 - keep]]
            aperiods = [aperiods[keep] or aperiods[1 - keep]]
            vs[0].samples[0]['GT'] = '1/1'

    for si, sv in enumerate(vs):
        period = periods[si]
        if not period:
            continue
        results = []
        inverted = {}
        for j, p in enumerate(period):
            try:
                result = match_ref(sv, ref, p, mt, units[si])
            except Exception as e:
                log.warning("match_ref failed for %s:%d period=%d: %s", sv.chrom, sv.pos, p, e)
                result = None
            if result:
                results.append((j, *result))
                inverted[j] = mm2_match.last_inverted
        indel_len = len(sv.alt[0]) if is_insert[si] else len(sv.ref) - len(sv.alt[0])
        # The whole insertion did not match locally as one block (longer than the array,
        # or only a mappy partial match): try the largest unit that fits the array.
        if is_insert[si] and units[si] and not any(period[r[0]] == indel_len and r[3] == (0, indel_len) for r in results):
            fit = None
            if aperiods[si]:                    # the array-level period of a microsatellite-rich indel
                for p in (aperiods[si], 2 * aperiods[si]):
                    if p < indel_len:
                        m = local_match(sv, ref, p, mt, p)
                        if m:
                            fit = (p, m)
                            break
            if fit is None:
                # bound the ladder by twice the array already measured (a larger unit can
                # extend further than the smaller candidates did, but not without limit)
                max_ext = max((r[2][1] - r[2][0] for r in results), default=None)
                fit = _fit_unit(sv, ref, units[si], indel_len, mt, 2 * max_ext if max_ext else None)
            if fit:
                period.append(fit[0])
                results.append((len(period) - 1, *fit[1]))
        # Filter: matched portion must be large enough and a meaningful fraction of period
        results = [r for r in results if r[3][1] - r[3][0] >= MIN_SEQ_LEN
                   and (r[3][1] - r[3][0]) / period[r[0]] >= MIN_PERIOD_RATIO]
        if not results:
            continue
        result = select_candidate(results, period, indel_len)

        ri, best_match, ref_range, seq_range, best_ratio = result
        best_period = period[ri]
        if is_insert[si]:
            period_cnt = len(sv.alt[0]) / best_period
        else:
            overlap = max(0, min(sv.pos + len(sv.ref), ref_range[1]) - max(sv.pos + 1, ref_range[0]))
            period_cnt = overlap / best_period
        period_cnt = max(1, int(period_cnt + 0.5))
        seq_length = seq_range[1] - seq_range[0]
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
                'UNIT': units[si] or best_period,
                'REFSTART': ref_range[0] + 1,
                'REFSTOP': ref_range[1] + 1,
                'END': best_match[1],
                'INDEL': indel_pos}
        if inverted.get(ri):
            sv.info['INV'] = 1
        sv.pos = best_match[0] if is_dup else best_match[0] - 1   # a DEL's POS is the base before the deletion
        sv.end = best_match[1]
        sv.ref = ref.get(sv.chrom, sv.pos, sv.pos+1)
        sv.alt = ['<DUP>'] if is_dup else ['<DEL>']
        sv.line = None
        result_vs.append(sv)
    return result_vs, time.time() - stime

def proc_batch(ref, variants, thresh, extend_slack=0.0, match_thresh=None):
    global EXTEND_THRESH
    EXTEND_THRESH = thresh if match_thresh is not None else None
    global EXTEND_SLACK
    EXTEND_SLACK = extend_slack
    results = []
    for v in variants:
        try:
            results.append(proc_variant(v, ref, thresh, match_thresh))
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
    if len(v1.alt) != 1 or len(v2.alt) != 1:
        return None
    if len(v1.ref) > len(v1.alt[0]) or len(v2.ref) > len(v2.alt[0]):
        return None                      # both must be insertions
    ref_gap = ref.get(v1.chrom, v1.end, v2.pos+1)
    v1.ref += ref_gap
    v1.alt[0] += ref_gap
    v1.alt.append(v1.ref[:-1] + v2.alt[0])
    # allele 1 sits on the haplotype v1 was on, allele 2 on the other (v1 used to write 1|1,
    # which made proc_variant drop the second allele)
    v1.samples[0]['GT'] = '1|2' if v1.samples[0]['GT'] == '1|0' else '2|1'
    v1.end = v2.pos + 1
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
    num_periods = max(1, round(v.info['SVLEN'] / period))
    gt_sum = sum(gt_set(v.samples[0]['GT'])) or 1
    if cndiff % gt_sum:
        gt_sum = 1         # a merged count that is not per-haplotype: re-express it as a whole
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
    # keep the total change in base pairs: CNDIFF counts copies over all haplotypes
    v.info['CNDIFF'] = best_d * gt_sum if cndiff > 0 else -best_d * gt_sum
    v.info['SVLEN'] = svlen_periods * period
    v.pos = v.info['REFSTART'] - 1           # 0-based start of the array
    v.end = v.pos + v.info['SVLEN']
    v.info['END'] = v.end
    v.line = None
    return v

# Merge overlapping output CNV calls from the same tandem array.
# Handles: (1) opposite-phase same-type pairs -> hom, (2) opposite-type pairs -> net CNDIFF,
# (3) same-type DEL/DUP pairs with compatible periods -> combined CNDIFF.
def merge_output_vars(variants):
    """Records arrive sorted by POS. Each one is tried against every kept record whose array it
    overlaps, not only the previous one (a DUP formed by merging two DUPs must still net against a DEL
    that came before them), and a merged record is tried again until nothing more merges."""
    result = []
    for v in variants:
        cur = v
        while cur is not None:
            for i in range(len(result) - 1, -1, -1):
                merged = _try_merge_output(result[i], cur)
                if merged is None:
                    continue
                del result[i]
                cur = None if merged == 'drop' else merged
                break
            else:
                result.append(cur)
                cur = None
    # An inverted copy is used only to cancel the deletion side of an inversion (assemblies write an
    # inversion as a deletion plus an inverted insertion, sometimes in several pieces). On its own it is
    # not reported: read depth supports 19 of 45 such records on NA21102, i.e. none beyond chance.
    for v in result:
        v.info.pop('_BP', None)
    return [v for v in sorted(result, key=lambda v: v.pos) if not v.info.get('INV')]

def _event_bp(v):
    """Signed change in reference-matching base pairs over all haplotypes. CNDIFF counts copies
    of the matched unit: SVLEN, which equals PERIOD for a full local match, the matched part of
    the unit for a partial mappy match (the rest of the insertion is novel sequence and changes no
    depth), and the re-expressed unit after extend_v."""
    return v.info.get('_BP', v.info['CNDIFF'] * v.info['SVLEN'])


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
    # Opposite types (DUP+DEL or DEL+DUP) on one array net their change in base pairs. CNDIFF counts
    # copies of each record's own unit, so counts are only comparable through the unit. "One array"
    # is a reciprocal overlap of NET_RECIPROCAL: the two records' arrays were extended independently
    # and differ by more than 10% at ~30 pairs genome-wide that are plainly the same locus.
    if v1.alt[0] != v2.alt[0]:
        if min(r1, r2) < NET_RECIPROCAL:
            return None
        bp1, bp2 = _event_bp(v1), _event_bp(v2)
        net = bp1 + bp2
        if bool(v1.info.get('INV')) != bool(v2.info.get('INV')) and abs(net) > 0.1 * max(abs(bp1), abs(bp2)):
            return None       # an inverted copy cancels only a deletion of its own size
        # cancel when the residual is below the event floor or 5% of the larger event (an inversion
        # represented as a deletion plus an inverted insertion leaves a few unmatched bases per copy)
        if abs(net) < max(MIN_SEQ_LEN, EXACT_FRAC * max(abs(bp1), abs(bp2))):
            return 'drop'  # cancel out
        keep = v1 if abs(bp1) >= abs(bp2) else v2
        # count the residual in the finer unit of the two, bounded by the array: a one-copy deletion must
        # not vanish by rounding against a duplication re-expressed as one copy of seven units
        unit = min(v1.info['PERIOD'], v1.info['SVLEN'], v2.info['PERIOD'], v2.info['SVLEN'],
                   keep.info['REFSTOP'] - keep.info['REFSTART'])
        keep.info['PERIOD'] = keep.info['SVLEN'] = unit
        keep.info['CNDIFF'] = (1 if net > 0 else -1) * max(1, int(abs(net) / unit + 0.5))
        keep.info['_BP'] = net             # exact total, so a chain of merges rounds once
        keep = extend_v(keep)
        keep.line = None
        return keep
    # Same-type merges also need one array. A record nested inside a much larger array is a distinct event
    # (a 700 bp deletion of one sub-array inside a 720 kb satellite deletion); adding it would drive the
    # merged unit down to the smaller record's and count the large deletion in hundreds of tiny units.
    if min(r1, r2) < NET_RECIPROCAL:
        return None
    # Same type DEL: two deletions on one array add up in base pairs (the two alleles of a 1|2
    # site are two events, one per haplotype; the old rule kept one of them when the GTs matched)
    if v1.alt[0] == '<DEL>' and v2.alt[0] == '<DEL>':
        keep, other = (v1, v2) if len1 >= len2 else (v2, v1)
        # count in the finer unit so a nested unequal deletion is not rounded away
        unit = min(v1.info['PERIOD'], v1.info['SVLEN'], v2.info['PERIOD'], v2.info['SVLEN'],
                   keep.info['REFSTOP'] - keep.info['REFSTART'])
        total = _event_bp(keep) + _event_bp(other)
        gt1, gt2 = gt_set(v1.samples[0]['GT']), gt_set(v2.samples[0]['GT'])
        dot = '.' in v1.samples[0]['GT'] or '.' in v2.samples[0]['GT']
        if not dot and len(gt1) == len(gt2):
            keep.samples[0]['GT'] = '|'.join(str(g) for g in (gt1[i] | gt2[i] for i in range(len(gt1))))
        elif gt1 and gt2:
            keep.samples[0]['GT'] = str(gt1[0] | gt2[0]) + '/.'
        keep.info['PERIOD'] = keep.info['SVLEN'] = unit
        keep.info['CNDIFF'] = -max(1, int(abs(total) / unit + 0.5))
        keep.info['_BP'] = total
        keep.line = None
        return keep
    # Same type DUP with compatible periods (multiples of each other, or the same fundamental unit)
    if v1.alt[0] == '<DUP>' and v2.alt[0] == '<DUP>':
        if bool(v1.info.get('INV')) != bool(v2.info.get('INV')):
            return None       # an inverted copy is not reported; keep it out of reported records
        p1, p2 = v1.info['PERIOD'], v2.info['PERIOD']
        r = max(p1, p2) / min(p1, p2) if min(p1, p2) > 0 else 999
        u1, u2 = v1.info.get('UNIT', 0), v2.info.get('UNIT', 0)
        if abs(r - round(r)) > 0.1 and not (u1 and u2 and abs(u1 / u2 - 1) <= 0.05):
            log.debug("Incompatible periods %d vs %d at %s:%d", p1, p2, v1.chrom, v1.pos)
            return None
        keep, other = (v1, v2) if len1 >= len2 else (v2, v1)
        # the common unit: the finer of the two (SVLEN where a partial match carries PERIOD > SVLEN), or
        # the array when a partially matched block exceeds it (depth over an array smaller than the unit
        # changes by base pairs / array)
        period = min(p1, v1.info['SVLEN'], p2, v2.info['SVLEN'], keep.info['REFSTOP'] - keep.info['REFSTART'])
        cd1 = int(_event_bp(keep) / period + 0.5)
        cd2 = int(_event_bp(other) / period + 0.5)
        gt1, gt2 = gt_set(v1.samples[0]['GT']), gt_set(v2.samples[0]['GT'])
        if gt1 == gt2 or (cd1 == cd2 and [1, 1] not in (gt1, gt2)):
            dot = '.' in v1.samples[0]['GT'] or '.' in v2.samples[0]['GT']
            if not dot:
                gt = [gt1[i] | gt2[i] for i in range(min(len(gt1), len(gt2)))]
                keep.samples[0]['GT'] = '|'.join(str(g) for g in gt)
            else:
                keep.samples[0]['GT'] = str(gt1[0] | gt2[0]) + '/.'
            total = _event_bp(keep) + _event_bp(other)      # before the unit changes
            keep.info['PERIOD'] = keep.info['SVLEN'] = period
            keep.info['CNDIFF'] = max(1, int(total / period + 0.5))
            keep.info['_BP'] = total
            keep = extend_v(keep)
            keep.line = None
            return keep
    return None

def tool_version():
    """COMBINE_SV_CNV_VERSION when a driver set it, else combine_sv_cnv-v<INDEL2CNV_VERSION>, followed by
    ' (sentieon-cli-<version>)' when this file sits in the sentieon_cli package (sentieon_cli/scripts/)."""
    if os.environ.get('COMBINE_SV_CNV_VERSION'):
        return os.environ['COMBINE_SV_CNV_VERSION']
    ver = 'combine_sv_cnv-v' + INDEL2CNV_VERSION
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

def provenance_line():
    """The ##CommandLine header line, in the form of Sentieon's ##SentieonCommandLine."""
    date = datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
    cmd = ' '.join(sys.argv).replace('\\', '\\\\').replace('"', '\\"')
    return '##CommandLine.indel2cnv=<ID=indel2cnv,Version="%s",Date="%s",CommandLine="%s">' % (tool_version(), date, cmd)

def main(args):
    global EXTEND_SLACK, EXTEND_THRESH
    EXTEND_SLACK = args.extend_slack
    EXTEND_THRESH = args.thresh if args.match_thresh is not None else None
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=log_level, format='%(asctime)s %(levelname)s %(message)s',
                        datefmt='%H:%M:%S')
    ref = Reference(args.ref)
    input_vcf = vcflib.VCF(args.input_vcf, 'r')
    fout = vcflib.VCF(args.out_vcf, 'w')
    update = ('##INFO=<ID=POSSHIFT,Number=1,Type=Integer,Description="Best matching CNV sequence position shifted from INDEL position">',
              '##INFO=<ID=POSSHIFT2,Number=1,Type=Integer,Description="Distance from original INDEL position to CNV sequence reference boundary">',
              '##INFO=<ID=CNDIFF,Number=1,Type=Integer,Description="Copy number state difference">',
              '##INFO=<ID=DIST,Number=1,Type=Float,Description="Levenshtein similarity ratio of the matched unit to its reference copy, 1 - edits/(2 * length); identity is about 2 * DIST - 1">',
              '##INFO=<ID=PERIOD,Number=1,Type=Integer,Description="Length of the repeating unit">',
              '##INFO=<ID=UNIT,Number=1,Type=Integer,Description="Fundamental period of the INDEL sequence (PERIOD when none found)">',
              '##INFO=<ID=REFSTART,Number=1,Type=Integer,Description="Reference start locus of the tandem repeat array">',
              '##INFO=<ID=REFSTOP,Number=1,Type=Integer,Description="Reference stop locus of the tandem repeat array">',
              '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SVLEN">',
              '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SVTYPE">',
              '##ALT=<ID=DEL,Description="Deletion relative to the reference">',
              '##ALT=<ID=DUP,Description="Duplication relative to the reference">',
              '##INFO=<ID=END,Number=1,Type=Integer,Description="END location">',
              '##INFO=<ID=INDEL,Number=1,Type=String,Description="Source INDEL locus (chrom:pos)">',
              provenance_line(),
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
    n_skip_in = n_skip_out = 0
    skipped_spans = collections.defaultdict(list)   # gate 16: what the skipped records covered
    for c in contigs:
        vcf_c = input_vcf.range(c)
        for v in vcf_c:
            if (len(v.filter) == 0 or 'PASS' in v.filter) and max(len(a) for a in [v.ref] + v.alt) >= MIN_SEQ_LEN:
                if len(v.ref) >= max((len(a) for a in v.alt), default=0):
                    nf = n_frac(v.ref)
                    if nf >= MAX_N_FRAC:
                        log.info('skipped N-gap record %s:%d len=%d nfrac=%.2f', v.chrom, v.pos + 1, len(v.ref), nf)
                        n_skip_in += 1
                        # the deleted span: the REF allele without its anchor base, 0-based half-open
                        skipped_spans[v.chrom].append((v.pos + 1, v.pos + len(v.ref), 'in:%.2f' % nf))
                        continue
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
            try:
                vs, t = proc_variant(v, ref, args.thresh, args.match_thresh)
            except Exception as e:
                log.error("Failed processing %s:%d: %s", v.chrom, v.pos, e)
                continue
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
            args_in.append((ref, variants[i:i+batch_size], args.thresh, args.extend_slack, args.match_thresh))
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
            nf = n_frac(ref.get(v.chrom, v.info['REFSTART'] - 1, v.info['REFSTOP'] - 1))
            if nf >= MAX_N_FRAC:
                log.info('skipped N-gap record %s:%d len=%d nfrac=%.2f', v.chrom, v.pos + 1,
                         v.info['REFSTOP'] - v.info['REFSTART'], nf)
                n_skip_out += 1
                skipped_spans[v.chrom].append((v.info['REFSTART'] - 1, v.info['REFSTOP'] - 1, 'out:%.2f' % nf))
                continue
            fout.emit(v)
    t1 = time.time()
    mm, ut, st = 0, 0, 0
    mem_scale = 1 if platform.system() == 'Darwin' else 1024
    for who in (resource.RUSAGE_SELF, resource.RUSAGE_CHILDREN):
        ru = resource.getrusage(who)
        mm += ru.ru_maxrss * mem_scale
        ut += ru.ru_utime
        st += ru.ru_stime
    log.info('N-gap records skipped: %d on input, %d on output', n_skip_in, n_skip_out)
    if args.skipped_bed:
        with open(args.skipped_bed, 'w') as fb:
            for c in contigs:
                for bs, be, name in sorted(skipped_spans.get(c, ())):
                    fb.write('%s\t%d\t%d\t%s\n' % (c, bs, be, name))
        log.info('wrote the skipped spans to %s', args.skipped_bed)
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
    parser.add_argument('--thresh', help="Minimum Levenshtein similarity ratio of a match, 1 - edits/(2 * length), applied identically to the local and mappy paths (0.9 is about 80%% identity)", default=0.9, type=float)
    parser.add_argument('--extend_slack', help="Array extension continues down to thresh minus this", default=0.0, type=float)
    parser.add_argument('--match_thresh', help="Identity a match must reach (default: thresh); period verification and the array walk keep thresh", default=None, type=float)
    parser.add_argument('--skipped_bed', help="Write the spans of the records gate 16 skipped to this bed, "
                        "0-based half-open, name in|out:nfrac: the deleted span for an input skip and the "
                        "array for an output skip. Subtracting it from the high-confidence bed makes a "
                        "mostly-N record's real flank count as not assessed")
    parser.add_argument('-v', '--verbose', action='store_true', help="Verbose logging")
    args = parser.parse_args()
    main(args)
