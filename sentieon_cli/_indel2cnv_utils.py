"""Pure-Python replacements for ``Levenshtein.ratio`` and
``scipy.signal.find_peaks`` used by ``scripts/indel2cnv.py``.

Kept in the package (not ``scripts/``) so pytest can import them —
``scripts/`` is excluded from collection by ``conftest.py``.
"""

from __future__ import annotations

from collections import Counter
from typing import Sequence


def ratio(a: str, b: str, threshold: float = 0.0) -> float:
    """Return ``2 * LCS(a, b) / (len(a) + len(b))``.

    Equivalent to ``Levenshtein.ratio(a, b)`` (which uses insert/delete
    cost 1 and substitution cost 2, making the underlying distance
    equal to ``|a| + |b| - 2 * LCS(a, b)``).

    When ``threshold > 0`` and a cheap upper bound proves the ratio
    cannot reach ``threshold``, returns ``0.0`` as a sentinel without
    running the full LCS.  Otherwise returns the exact ratio.
    """
    la, lb = len(a), len(b)
    if la == 0 and lb == 0:
        return 1.0
    if la == 0 or lb == 0:
        return 0.0
    if a == b:
        return 1.0

    tot = la + lb

    # Length-based upper bound: 2*LCS <= 2*min(|a|,|b|).
    if threshold > 0.0 and 2 * min(la, lb) < threshold * tot:
        return 0.0

    # Character-multiset intersection is a true upper bound on 2*LCS.
    ub_matches = sum((Counter(a) & Counter(b)).values())
    if ub_matches == 0:
        return 0.0
    if threshold > 0.0 and (2.0 * ub_matches) < threshold * tot:
        return 0.0

    return 2.0 * _bit_parallel_lcs(a, b) / tot


def _bit_parallel_lcs(a: str, b: str) -> int:
    """LCS length via the Allison-Dix bit-parallel algorithm.

    Each iteration is O(len(a) / word_size) big-int operations; Python
    int ops exploit this implicitly.  Puts the shorter string in ``a``
    to minimize bitmap width.
    """
    if len(a) > len(b):
        a, b = b, a
    m = len(a)

    char_mask: dict[str, int] = {}
    for i, ch in enumerate(a):
        char_mask[ch] = char_mask.get(ch, 0) | (1 << i)

    full = (1 << m) - 1
    s = full
    for ch in b:
        t = char_mask.get(ch)
        if t is None:
            continue
        u = s & t
        if u == 0:
            continue
        s = ((s + u) | (s - u)) & full
    return m - bin(s).count("1")


def find_peaks(
    values: Sequence[float], *, distance: int, height: float
) -> list[int]:
    """Return indices of local maxima >= ``height`` with no two kept
    peaks closer than ``distance`` samples.

    Matches ``scipy.signal.find_peaks(x, distance=d, height=h)`` for the
    one call site that uses it.  Plateaus map to their midpoint index
    (``(start + end) // 2``), matching scipy's behavior.  Boundary
    values are never peaks.
    """
    n = len(values)
    candidates: list[int] = []
    i = 1
    while i < n - 1:
        if values[i] <= values[i - 1] or values[i] < height:
            i += 1
            continue
        j = i
        while j + 1 < n and values[j + 1] == values[i]:
            j += 1
        if j + 1 < n and values[j + 1] < values[i]:
            candidates.append((i + j) // 2)
        i = j + 1

    if distance <= 1 or len(candidates) <= 1:
        return candidates

    # Keep highest peaks first; drop any within `distance` of one
    # already kept.  Tie-break on index so earlier wins.
    ordered = sorted(candidates, key=lambda k: (-values[k], k))
    kept: list[int] = []
    for p in ordered:
        if all(abs(p - q) >= distance for q in kept):
            kept.append(p)
    kept.sort()
    return kept
