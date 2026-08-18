"""BSRT scoring and scheduling, in pure Python.

PORT OF scoring.js. The browser, desktop and PsychoPy builds are meant to be
one instrument with three front ends, so this file is a deliberate line-for-line
port rather than an independent implementation — a second interpretation of the
spec would drift, and the drift would show up as a systematic difference between
data collected in different rooms.

`tests/test_equivalence.py` drives this module and scoring.js with identical
inputs and asserts the outputs match to 1e-9, including the seeded schedule
generator. If you change one, change the other and re-run it.

Deliberately stdlib-only: no numpy, no PsychoPy. Everything here can be tested
without a display, and the experiment script imports it rather than the reverse.

Trial classification, per the BSRT/OSLER specification:
    RT <= hitWindowMs (1000 ms)         -> HIT
    hitWindowMs < RT < isiMs (3000 ms)  -> MISS, but the RT is still recorded
    no response at all                  -> MISS, no RT

A LAPSE is a slow HIT: lapseMs < RT <= hitWindowMs. It is a property of a hit,
not a separate outcome, so lapses count inside the hit total and never
contribute to the sleep-onset criterion. Only misses do.
"""

import math

# ---------------- basic statistics ----------------
#
# All of these return None rather than raising on an empty input, matching the
# JS, because a minute with no hits in it is normal data and not an error.


def mean(a):
    return sum(a) / len(a) if a else None


def median(a):
    if not a:
        return None
    b = sorted(a)
    m = len(b) // 2
    return b[m] if len(b) % 2 else (b[m - 1] + b[m]) / 2.0


def sd(a):
    """Sample standard deviation (n-1). None below two observations."""
    if len(a) < 2:
        return None
    m = mean(a)
    return math.sqrt(sum((x - m) ** 2 for x in a) / (len(a) - 1))


def decile_mean(sorted_asc, tail):
    """Mean of the top or bottom decile of an ascending-sorted list.

    Named 'low'/'high' by POSITION, not 'fastest'/'slowest': for reaction time
    the fastest responses are the low tail, but for reaction speed they are the
    high tail, and naming by position removes the chance of silently inverting
    one of them. The decile is at least one observation.
    """
    if not sorted_asc:
        return None
    k = max(1, round_half_up(len(sorted_asc) * 0.1))
    part = sorted_asc[:k] if tail == 'low' else sorted_asc[len(sorted_asc) - k:]
    return mean(part)


def round_half_up(x):
    """JavaScript's Math.round: halves go up, not to even.

    Python's built-in round() is banker's rounding, so round(0.5) is 0 and
    round(2.5) is 2. A decile of exactly n*0.1 = 2.5 observations would then
    pick a different subset in the two implementations.
    """
    return int(math.floor(x + 0.5))


def ceil_decile(sorted_asc, tail):
    """The normative convention's decile size: ceil, not round."""
    if not sorted_asc:
        return None
    k = max(1, math.ceil(len(sorted_asc) * 0.1))
    part = sorted_asc[:k] if tail == 'low' else sorted_asc[len(sorted_asc) - k:]
    return mean(part)


def to_rs(rts):
    """Reaction speed, 1000/RT. Defined only for RT > 0, so anticipations at or
    below zero are excluded rather than producing an infinity."""
    return [1000.0 / r for r in rts if r > 0]


# ---------------- corrections ----------------


def apply_correction(rts, opts):
    """Exclusions, in this order: false starts, then outliers.

    Order matters. False starts are removed first so they cannot drag the mean
    and SD that define the outlier bounds.
    """
    kept = list(rts)
    n_false_starts = 0
    n_outliers = 0
    lo = hi = None

    if opts.get('removeFalseStarts'):
        before = len(kept)
        kept = [v for v in kept if v >= opts['falseStartMs']]
        n_false_starts = before - len(kept)

    if opts.get('removeOutliers') and len(kept) >= 3:
        m, s = mean(kept), sd(kept)
        if s is not None and s > 0:
            lo = m - opts['sdMultiplier'] * s
            hi = m + opts['sdMultiplier'] * s
            before2 = len(kept)
            kept = [v for v in kept if lo <= v <= hi]
            n_outliers = before2 - len(kept)

    return {'kept': kept, 'nFalseStarts': n_false_starts, 'nOutliers': n_outliers,
            'lo': lo, 'hi': hi}


# ---------------- one block of trials ----------------


def block_metrics(hit_rts, opts):
    """`hit_rts` are the raw reaction times of the HITS in this block.

    Misses contribute no RT, and late responses (RT > hit window) are excluded
    because they are misses. Miss and lapse counts are passed in separately.
    """
    srt = sorted(hit_rts)
    corr = apply_correction(hit_rts, opts)
    c_sorted = sorted(corr['kept'])

    rs = to_rs(hit_rts)
    rs_sorted = sorted(rs)
    c_rs = to_rs(corr['kept'])
    c_rs_sorted = sorted(c_rs)

    def ipr(a, b):
        """Interpercentile range: the distance between the two decile means.

        A within-test stability measure — a participant whose best and worst
        tenths sit far apart is performing erratically even when the average
        looks intact — and it is what the normative tables report.
        """
        return None if a is None or b is None else abs(a - b)

    fast_rt = decile_mean(srt, 'low')
    slow_rt = decile_mean(srt, 'high')
    c_fast_rt = decile_mean(c_sorted, 'low')
    c_slow_rt = decile_mean(c_sorted, 'high')
    fast_rs = decile_mean(rs_sorted, 'high')
    slow_rs = decile_mean(rs_sorted, 'low')
    c_fast_rs = decile_mean(c_rs_sorted, 'high')
    c_slow_rs = decile_mean(c_rs_sorted, 'low')

    return {
        'n': len(hit_rts),
        'nCorrected': len(corr['kept']),
        'nFalseStartsRemoved': corr['nFalseStarts'],
        'nOutliersRemoved': corr['nOutliers'],
        'outlierLoMs': corr['lo'],
        'outlierHiMs': corr['hi'],
        'nRsUndefined': len(hit_rts) - len(rs),

        'avgRt': mean(hit_rts),
        'medianRt': median(hit_rts),
        'sdRt': sd(hit_rts),
        'fastest10Rt': fast_rt,
        'slowest10Rt': slow_rt,
        'iprRt': ipr(slow_rt, fast_rt),

        'corrAvgRt': mean(corr['kept']),
        'corrMedianRt': median(corr['kept']),
        'corrSdRt': sd(corr['kept']),
        'corrFastest10Rt': c_fast_rt,
        'corrSlowest10Rt': c_slow_rt,
        'corrIprRt': ipr(c_slow_rt, c_fast_rt),

        'avgRs': mean(rs),
        'medianRs': median(rs),
        'sdRs': sd(rs),
        'fastest10Rs': fast_rs,
        'slowest10Rs': slow_rs,
        'iprRs': ipr(fast_rs, slow_rs),

        'corrAvgRs': mean(c_rs),
        'corrMedianRs': median(c_rs),
        'corrSdRs': sd(c_rs),
        'corrFastest10Rs': c_fast_rs,
        'corrSlowest10Rs': c_slow_rs,
        'corrIprRs': ipr(c_fast_rs, c_slow_rs),
    }


# ---------------- error profiles ----------------


def error_profiles(outcomes, miss_criterion):
    """Runs of consecutive misses, binned by run length.

    The frequency reported is the number of RUNS in each band, not the number
    of missed trials. A run of `miss_criterion` or more terminates the test, so
    ep7plus is normally 0 or 1.
    """
    runs, cur = [], 0
    for o in outcomes:
        if o == 'miss':
            cur += 1
        elif cur:
            runs.append(cur)
            cur = 0
    if cur:
        runs.append(cur)

    ep = {'ep1_2': 0, 'ep3_6': 0, 'ep7plus': 0, 'runs': runs, 'longestRun': 0}
    for r in runs:
        ep['longestRun'] = max(ep['longestRun'], r)
        if r <= 2:
            ep['ep1_2'] += 1
        elif r <= 6:
            ep['ep3_6'] += 1
        else:
            ep['ep7plus'] += 1
    # A criterion of 0 or less means 'no early termination' (PVT runs to time).
    ep['criterionReached'] = miss_criterion > 0 and ep['longestRun'] >= miss_criterion
    return ep


# ---------------- velocity and acceleration ----------------


def ols_slope(series):
    """Least-squares slope against index; None under two usable points."""
    xs, ys = [], []
    for i, v in enumerate(series):
        if v is not None:
            xs.append(i)
            ys.append(v)
    if len(xs) < 2:
        return None
    mx, my = mean(xs), mean(ys)
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den = sum((x - mx) ** 2 for x in xs)
    return None if den == 0 else num / den


def dynamics(per_minute, key):
    """Rate of change across minutes, and the rate of change of that.

    Both are undefined for the first minute (and acceleration for the second),
    reported as None rather than zero so the gap stays visible in the data.
    """
    series = [m[key] for m in per_minute]
    velocity = []
    for i in range(len(series)):
        velocity.append(None if i == 0 or series[i] is None or series[i - 1] is None
                        else series[i] - series[i - 1])
    acceleration = []
    for j in range(len(velocity)):
        acceleration.append(None if j < 2 or velocity[j] is None or velocity[j - 1] is None
                            else velocity[j] - velocity[j - 1])
    return {'velocity': velocity, 'acceleration': acceleration, 'slope': ols_slope(series)}


# ---------------- stimulus schedule ----------------


MASK32 = 0xFFFFFFFF


def make_rng(seed):
    """mulberry32, bit-identical to the JS.

    A seeded generator rather than random.random() so a schedule can be
    reproduced exactly from the seed recorded with the data, and so the same
    seed produces the same schedule in all three builds — which is what makes a
    participant's session reproducible from the CSV alone.

    JS semantics reproduced: `|0` and `>>>` are 32-bit, `Math.imul` is a 32-bit
    multiply, and the `+` before the final xor goes through ToInt32. Working in
    unsigned mod-2^32 throughout gives the same bit patterns.
    """
    a = seed & MASK32

    def rng():
        nonlocal a
        a = (a + 0x6D2B79F5) & MASK32
        t = ((a ^ (a >> 15)) * (1 | a)) & MASK32
        t = (t + (((t ^ (t >> 7)) * (61 | t)) & MASK32)) & MASK32 ^ t
        return ((t ^ (t >> 14)) & MASK32) / 4294967296.0

    return rng


def shuffle(arr, rng):
    """Fisher-Yates, unbiased, driven by the seeded generator."""
    a = list(arr)
    for i in range(len(a) - 1, 0, -1):
        j = int(math.floor(rng() * (i + 1)))
        a[i], a[j] = a[j], a[i]
    return a


def build_schedule(opts):
    """The inter-stimulus intervals for a whole trial.

    PVT mode: the specified set (2, 4, 6, 8, 10 s) sums to exactly 30 s, so a
    30-second block holds precisely one of each interval and the schedule is a
    fresh random PERMUTATION per block — random in order, balanced in
    composition, so no block can come out all-short or all-long. `method`
    records the rule used, because a custom set that does not sum to the block
    length cannot be balanced and falls back to sampling with replacement.
    """
    max_ms = opts['maxMs']

    if opts.get('mode') != 'pvt':
        n = max(1, math.ceil(max_ms / opts['isiMs']))
        return finish_schedule([opts['isiMs']] * n, max_ms, opts['blockMs'],
                               'fixed', opts['seed'], [opts['isiMs']])

    iset = list(opts['isiSetMs'])
    block_ms = opts['blockMs']
    balanced = sum(iset) == block_ms
    rng = make_rng(opts['seed'])

    isis = []
    total = 0
    guard = 0
    while total < max_ms and guard < 100000:
        guard += 1
        if balanced:
            for v in shuffle(iset, rng):
                if total >= max_ms:
                    break
                isis.append(v)
                total += v
        else:
            used = 0
            while used < block_ms and total < max_ms:
                v = iset[int(math.floor(rng() * len(iset)))]
                isis.append(v)
                used += v
                total += v

    return finish_schedule(isis, max_ms, block_ms,
                           'block_permutation' if balanced else 'sampled_with_replacement',
                           opts['seed'], iset)


def finish_schedule(isis, max_ms, block_ms, method, seed, interval_set):
    # Keep only whole intervals that fit inside the ceiling.
    kept = []
    total = 0
    for v in isis:
        if total + v > max_ms:
            break
        kept.append(v)
        total += v

    # Interval k spans [start_k, start_k + kept[k]) and its stimulus fires at
    # the END of it, so an M-interval schedule yields M stimuli, not M-1.
    #
    # epochIsi is this stimulus's own response window — normally the interval
    # leading into the NEXT stimulus. The last stimulus has no next interval to
    # borrow, so it gets the longest interval in the CONFIGURED set rather than
    # whatever was drawn last: fixed by the settings, not by the seed, so every
    # participant gets the same trial length whatever the shuffle produced.
    trailing_ms = max(interval_set)
    onsets, blocks, minutes, isi_before, epoch_isi = [], [], [], [], []
    acc = 0
    for k, v in enumerate(kept):
        blocks.append(int(acc // block_ms))
        # A minute is the same rule with a 60 s window. Both come from the
        # INTENDED schedule, never a measured onset: a millisecond of
        # presentation jitter would otherwise flip a boundary stimulus into the
        # neighbouring bucket and unbalance the counts.
        minutes.append(int(acc // 60000))
        acc += v
        onsets.append(acc)
        isi_before.append(v)
        epoch_isi.append(kept[k + 1] if k + 1 < len(kept) else trailing_ms)

    return {
        'isis': kept,
        'onsets': onsets,
        'blocks': blocks,
        'minutes': minutes,
        'isiBefore': isi_before,
        'epochIsi': epoch_isi,
        'method': method,
        'seed': seed,
        'blockMs': block_ms,
        'nStimuli': len(onsets),
        'leadInMs': kept[0] if kept else 0,
        'plannedDurationMs': onsets[-1] + epoch_isi[-1] if onsets else 0,
    }


# ---------------- response integrity ----------------


def detect_cheating(presses, n_trials, cfg):
    """Three independent signals that responding was not genuine.

    Continuous or repeated pressing is a way to fake alertness: hold the key
    down, or tap constantly, and every stimulus gets a fast "response". All
    three signals are reported so an experimenter can see WHY a trial was
    flagged. This never alters scoring and never excludes data by itself.
    """
    burst_window_ms = cfg.get('burstWindowMs', 1000)
    burst_limit = cfg.get('burstLimit', 5)
    rapid_gap_ms = cfg.get('rapidGapMs', 200)
    extra_rate_limit = cfg.get('extraRateLimit', 0.2)

    times = sorted(p['tMs'] for p in presses)

    burst_max = 0
    lo = 0
    for hi in range(len(times)):
        while times[hi] - times[lo] > burst_window_ms:
            lo += 1
        burst_max = max(burst_max, hi - lo + 1)

    rapid_pairs = sum(1 for i in range(1, len(times))
                      if times[i] - times[i - 1] < rapid_gap_ms)

    seen = set()
    extra_presses = 0
    for p in presses:
        k = p.get('epochIndex')
        if k is None:
            extra_presses += 1
        elif k in seen:
            extra_presses += 1
        else:
            seen.add(k)

    extra_rate = extra_presses / n_trials if n_trials else 0
    reasons = []
    if burst_max >= burst_limit:
        reasons.append('%d presses within %d ms' % (burst_max, burst_window_ms))
    if n_trials >= 5 and extra_rate > extra_rate_limit:
        reasons.append('%d%% more presses than stimuli' % round_half_up(extra_rate * 100))
    if len(times) > 10 and rapid_pairs / len(times) > 0.25:
        reasons.append('%d presses under %d ms apart' % (rapid_pairs, rapid_gap_ms))

    return {
        'totalPresses': len(times),
        'extraPresses': extra_presses,
        'extraRate': extra_rate,
        'burstMax': burst_max,
        'rapidPairs': rapid_pairs,
        'suspected': bool(reasons),
        'reasons': reasons,
        'thresholds': {'burstWindowMs': burst_window_ms, 'burstLimit': burst_limit,
                       'rapidGapMs': rapid_gap_ms, 'extraRateLimit': extra_rate_limit},
    }


# ---------------- main entry point ----------------


def score(epochs, cfg, presses=None):
    """Score a whole trial.

    epochs: [{'index', 'onsetMs', 'rtMs', ...}] where rtMs is the raw reaction
            time of whatever response occurred in that epoch (however late),
            or None.
    presses: every keypress as {'tMs', 'epochIndex'}, for the integrity check.
    """
    opts = {
        'removeFalseStarts': bool(cfg.get('removeFalseStarts')),
        'removeOutliers': bool(cfg.get('removeOutliers')),
        'falseStartMs': 100 if cfg.get('falseStartMs') is None else cfg['falseStartMs'],
        'sdMultiplier': 2 if cfg.get('sdMultiplier') is None else cfg['sdMultiplier'],
    }
    hit_window = cfg['hitWindowMs']

    def minute_of(e):
        """The minute comes from the schedule, so it is a planned bucket rather
        than something re-derived from a measured onset. The fallback is for
        epochs built without one: bucket by where the PRECEDING INTERVAL
        started, matching the schedule's rule. Bucketing on the onset alone
        pushes every boundary stimulus a minute late.
        """
        if e.get('minute') is not None:
            return e['minute']
        start = e['onsetMs'] if e.get('isiBeforeMs') is None else e['onsetMs'] - e['isiBeforeMs']
        return int(max(0, start) // 60000)

    trials = []
    for e in epochs:
        rt = e.get('rtMs')
        is_hit = rt is not None and rt <= hit_window
        trials.append({
            'index': e.get('index'),
            'onsetMs': e.get('onsetMs'),
            'minute': minute_of(e),
            'block': e.get('block'),
            'epochIsiMs': e.get('epochIsiMs'),
            'isiBeforeMs': e.get('isiBeforeMs'),
            'rtMs': rt,                                    # raw, always preserved
            'rsPerSec': 1000.0 / rt if rt is not None and rt > 0 else None,
            'outcome': 'hit' if is_hit else 'miss',
            'lateResponse': 1 if rt is not None and not is_hit else 0,
            'lapse': 1 if is_hit and rt > cfg['lapseMs'] else 0,
            'falseStart': 1 if rt is not None and rt < opts['falseStartMs'] else 0,
        })

    hits = [t for t in trials if t['outcome'] == 'hit']
    misses = [t for t in trials if t['outcome'] == 'miss']
    lapses = [t for t in trials if t['lapse']]
    late = [t for t in trials if t['lateResponse']]
    false_starts = [t for t in trials if t['falseStart']]

    last_minute = trials[-1]['minute'] if trials else -1
    per_minute = []
    for m in range(last_minute + 1):
        in_min = [t for t in trials if t['minute'] == m]
        min_hits = [t for t in in_min if t['outcome'] == 'hit']
        block = block_metrics([t['rtMs'] for t in min_hits], opts)
        block['minute'] = m + 1                            # 1-based for reporting
        block['trials'] = len(in_min)
        block['hits'] = len(min_hits)
        block['misses'] = len(in_min) - len(min_hits)
        block['lapses'] = sum(t['lapse'] for t in in_min)
        block['lateResponses'] = sum(t['lateResponse'] for t in in_min)
        block['hitRatio'] = len(min_hits) / len(in_min) if in_min else None
        per_minute.append(block)

    totals = block_metrics([t['rtMs'] for t in hits], opts)
    totals['trials'] = len(trials)
    totals['hits'] = len(hits)
    totals['misses'] = len(misses)
    totals['lapses'] = len(lapses)
    totals['lateResponses'] = len(late)
    totals['falseStarts'] = len(false_starts)
    totals['hitRatio'] = len(hits) / len(trials) if trials else None

    ep = error_profiles([t['outcome'] for t in trials], cfg['missCriterion'])

    return {
        'config': {
            'hitWindowMs': hit_window,
            'isiMs': cfg.get('isiMs'),
            'lapseMs': cfg['lapseMs'],
            'missCriterion': cfg['missCriterion'],
            'removeFalseStarts': opts['removeFalseStarts'],
            'removeOutliers': opts['removeOutliers'],
            'falseStartMs': opts['falseStartMs'],
            'sdMultiplier': opts['sdMultiplier'],
        },
        'trials': trials,
        'perMinute': per_minute,
        'totals': totals,
        'errorProfiles': ep,
        'integrity': detect_cheating(presses or [], len(trials), cfg),
        'dynamicsRs': dynamics(per_minute, 'avgRs'),
        'dynamicsRt': dynamics(per_minute, 'avgRt'),
    }


# ---------------- normative comparison ----------------


def normative_summary(trials, minutes, opts=None):
    """Re-score a trial the way the normative tables were computed.

    The app's own scoring is NOT reused here, because the two disagree in ways
    that would bias every z-score: the reference keeps every response however
    slow, where the app calls anything past the hit window a miss, and the
    reference sizes a decile with ceil rather than round.
    """
    opts = opts or {}
    in_window = [t for t in trials
                 if (t['minute'] < minutes if t.get('minute') is not None
                     else t['onsetMs'] <= minutes * 60000)]
    if not in_window:
        return None

    fs_ms = opts.get('falseStartMs', 100)
    lapse_ms = opts.get('lapseMs', 500)

    responded = [t for t in in_window if t['rtMs'] is not None]
    misses = len(in_window) - len(responded)
    false_starts = sum(1 for t in responded if t['rtMs'] < fs_ms)
    valid = sorted(t['rtMs'] for t in responded if t['rtMs'] > fs_ms)

    ep = error_profiles(['miss' if t['rtMs'] is None else 'hit' for t in in_window],
                        opts.get('missCriterion', 0))

    out = {
        'minutes': minutes,
        'trials': len(in_window),
        'hits': len(in_window) - misses,
        'misses': misses,
        'hitRatio': (len(in_window) - misses) / len(in_window) * 100.0,
        'falseStarts': false_starts,
        'lapses': sum(1 for v in valid if v > lapse_ms),
        'ep12': ep['ep1_2'], 'ep36': ep['ep3_6'], 'ep7': ep['ep7plus'],
        'nValid': len(valid),
    }
    if not valid:
        return out

    rs = sorted(1000.0 / v for v in valid)
    fast_rt = ceil_decile(valid, 'low')
    slow_rt = ceil_decile(valid, 'high')
    fast_rs = ceil_decile(rs, 'high')
    slow_rs = ceil_decile(rs, 'low')

    out.update({
        'rtMean': mean(valid), 'rtMedian': median(valid),
        'rtSd': sd(valid) if len(valid) > 1 else None,
        'rtFast10': fast_rt, 'rtSlow10': slow_rt, 'rtIpr': abs(slow_rt - fast_rt),
        'rsMean': mean(rs), 'rsMedian': median(rs),
        'rsSd': sd(rs) if len(rs) > 1 else None,
        'rsFast10': fast_rs, 'rsSlow10': slow_rs, 'rsIpr': abs(fast_rs - slow_rs),
    })
    return out


def compare_to_norm(value, norm):
    """Position of one value against its reference cell.

    `dir` is +1 when higher is better and -1 when higher is worse, so `zWorse`
    is always signed the same way: positive means worse than the reference,
    whatever the variable measures.
    """
    if value is None or not norm or norm.get('mean') is None:
        return None
    delta_worse = norm['mean'] - value if norm['dir'] >= 0 else value - norm['mean']
    out = {'value': value, 'mean': norm['mean'], 'sd': norm['sd'], 'n': norm.get('n'),
           'dir': norm['dir'], 'deltaWorse': delta_worse, 'z': None,
           'band': 'green', 'degenerate': False}
    if norm['dir'] == 0:
        out['band'] = 'neutral'
        return out
    if not norm['sd']:
        # No spread to measure against: worse at all means outside the whole
        # reference sample, which is worth flagging, but never as a z-score.
        out['degenerate'] = True
        out['band'] = 'red' if delta_worse > 0 else 'green'
        return out
    out['z'] = delta_worse / norm['sd']
    out['band'] = 'red' if out['z'] > 2 else ('orange' if out['z'] > 1 else 'green')
    return out
