"""Session scoring for the BSRT normative tables.

This is a line-for-line port of `normativeSummary` in scoring.js. It exists so
the reference values are computed by exactly the rule the app applies to a
participant's own trials — anything else would put a systematic offset into
every z-score.

Raw-file conventions, verified against both workbooks:
  3000        no response before the next stimulus -> a MISS, not an RT
  RT < 100    false start; counted, then excluded from the RT statistics
  RT > 100    valid, whatever the latency (a late response is kept, as the app
              keeps it in normativeSummary)
"""

import math
import statistics as st

TIMEOUT = 3000
FALSE_START_MS = 100
LAPSE_MS = 500

VAR_KEYS = [
    'trials', 'hitRatio', 'hits', 'misses', 'lapses', 'falseStarts',
    'ep12', 'ep36', 'ep7',
    'rtMean', 'rtMedian', 'rtSd', 'rtFast10', 'rtSlow10', 'rtIpr',
    'rsMean', 'rsMedian', 'rsSd', 'rsFast10', 'rsSlow10', 'rsIpr',
]


def ceil_decile(sorted_asc, tail):
    """The workbook's decile size: ceil, not round (scoring.js ceilDecile)."""
    if not sorted_asc:
        return None
    k = max(1, math.ceil(len(sorted_asc) * 0.1))
    part = sorted_asc[:k] if tail == 'low' else sorted_asc[len(sorted_asc) - k:]
    return st.mean(part)


def error_profiles(outcomes):
    """Runs of consecutive misses, binned by run length (scoring.js)."""
    runs, cur = [], 0
    for o in outcomes:
        if o == 'miss':
            cur += 1
        elif cur:
            runs.append(cur)
            cur = 0
    if cur:
        runs.append(cur)
    ep = {'ep12': 0, 'ep36': 0, 'ep7': 0}
    for r in runs:
        if r <= 2:
            ep['ep12'] += 1
        elif r <= 6:
            ep['ep36'] += 1
        else:
            ep['ep7'] += 1
    return ep


def summarise(raw, minutes, trials_per_minute=20):
    """One session, over the first `minutes` minutes. None if it ran shorter.

    A session that stopped at sleep onset has no data past that point, so it
    contributes to the lengths it completed and to none after — returning None
    rather than a partial window keeps a shorter test from being averaged in as
    if it were a full one.
    """
    want = minutes * trials_per_minute
    if len(raw) < want:
        return None
    w = raw[:want]

    responded = [x for x in w if x != TIMEOUT]
    misses = len(w) - len(responded)
    false_starts = sum(1 for x in responded if x < FALSE_START_MS)
    valid = sorted(x for x in responded if x > FALSE_START_MS)

    ep = error_profiles(['miss' if x == TIMEOUT else 'hit' for x in w])

    out = {
        'trials': len(w),
        'hits': len(w) - misses,
        'misses': misses,
        'hitRatio': (len(w) - misses) / len(w) * 100.0,
        'falseStarts': false_starts,
        'lapses': sum(1 for v in valid if v > LAPSE_MS),
        'ep12': ep['ep12'], 'ep36': ep['ep36'], 'ep7': ep['ep7'],
    }
    if not valid:
        return out

    rs = sorted(1000.0 / v for v in valid)
    fast_rt = ceil_decile(valid, 'low')
    slow_rt = ceil_decile(valid, 'high')
    fast_rs = ceil_decile(rs, 'high')
    slow_rs = ceil_decile(rs, 'low')

    out.update({
        'rtMean': st.mean(valid),
        'rtMedian': st.median(valid),
        'rtSd': st.stdev(valid) if len(valid) > 1 else 0.0,
        'rtFast10': fast_rt,
        'rtSlow10': slow_rt,
        'rtIpr': abs(slow_rt - fast_rt),
        'rsMean': st.mean(rs),
        'rsMedian': st.median(rs),
        'rsSd': st.stdev(rs) if len(rs) > 1 else 0.0,
        'rsFast10': fast_rs,
        'rsSlow10': slow_rs,
        'rsIpr': abs(fast_rs - slow_rs),
    })
    return out


def read_sessions(path, sheet, first_col, last_col, hour_of):
    """Raw trials per session, trailing blanks trimmed.

    `hour_of(row)` returns the hour bin, so the two studies' different ways of
    recording when a session ran stay out of the scoring code.
    """
    import openpyxl
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb[sheet]
    out = []
    for r in ws.iter_rows(min_row=4, values_only=True):
        if r[0] is None:
            continue
        raw = list(r[first_col:last_col])
        k = len(raw)
        while k > 0 and raw[k - 1] is None:
            k -= 1
        if not k:
            continue
        assert all(x is not None for x in raw[:k]), 'gap inside a session'
        out.append({'id': r[0], 'hour': hour_of(r), 'raw': raw[:k]})
    return out
