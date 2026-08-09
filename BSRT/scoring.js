'use strict';

/* BSRT scoring.
 *
 * SHARED FILE — a byte-identical copy lives at BSRT-desktop/renderer/scoring.js.
 * Both apps must score identically, so change one and copy it to the other.
 * `npm run check:scoring` in BSRT-desktop asserts the two copies match.
 *
 * Trial classification, per the BSRT/OSLER specification:
 *
 *   RT <= hitWindowMs (1000 ms)          -> HIT
 *   hitWindowMs < RT < isiMs (3000 ms)   -> MISS, but the RT is still recorded
 *   no response at all                   -> MISS, no RT
 *
 * A LAPSE is a slow HIT: lapseMs < RT <= hitWindowMs. It is a property of a
 * hit, not a separate outcome, so lapses are counted inside the hit total and
 * never contribute to the sleep-onset criterion. Only misses do.
 *
 * Reaction speed is RS = 1000 / RT. It is defined only for RT > 0, so
 * anticipatory responses at or below zero are excluded from RS (and counted).
 *
 * Also builds the stimulus schedule, because the two modes differ only in how
 * the inter-stimulus intervals are generated:
 *   bsrt — a fixed interval, every epoch the same
 *   pvt  — intervals drawn from a set (2/4/6/8/10 s) varying within 30 s blocks
 */

(function (root, factory) {
  var api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.BSRTScoring = api;
})(typeof self !== 'undefined' ? self : this, function () {

  /* ---------------- basic statistics ---------------- */

  function mean(a) {
    if (!a.length) return null;
    var s = 0;
    for (var i = 0; i < a.length; i++) s += a[i];
    return s / a.length;
  }

  function median(a) {
    if (!a.length) return null;
    var b = a.slice().sort(function (x, y) { return x - y; });
    var m = b.length >> 1;
    return b.length % 2 ? b[m] : (b[m - 1] + b[m]) / 2;
  }

  /* Sample standard deviation (n-1). Null below two observations. */
  function sd(a) {
    if (a.length < 2) return null;
    var m = mean(a), s = 0;
    for (var i = 0; i < a.length; i++) s += (a[i] - m) * (a[i] - m);
    return Math.sqrt(s / (a.length - 1));
  }

  /*
   * Mean of the top or bottom decile of an ascending-sorted array.
   *
   * Deliberately named 'low'/'high' rather than 'fastest'/'slowest': for
   * reaction TIME the fastest responses are the low tail, but for reaction
   * SPEED they are the high tail. Naming the tail by position removes the
   * chance of silently inverting one of them.
   *
   * The decile size is at least one observation, so short blocks still report
   * a value rather than a null.
   */
  function decileMean(sortedAsc, tail) {
    if (!sortedAsc.length) return null;
    var k = Math.max(1, Math.round(sortedAsc.length * 0.1));
    var slice = tail === 'low' ? sortedAsc.slice(0, k) : sortedAsc.slice(sortedAsc.length - k);
    return mean(slice);
  }

  /* ---------------- corrections ---------------- */

  /*
   * Exclusions, applied in this order:
   *   1. false starts: RT < falseStartMs (default 100 ms)
   *   2. outliers: RT outside mean +/- sdMultiplier * SD
   *
   * Order matters. False starts are removed first so they cannot drag the
   * mean and SD that define the outlier bounds. The bounds are computed within
   * whatever set is being summarised, so a per-minute correction uses that
   * minute's own mean and SD, and a whole-test correction uses the whole test's.
   */
  function applyCorrection(rts, opts) {
    var kept = rts.slice();
    var nFalseStarts = 0, nOutliers = 0, lo = null, hi = null;

    if (opts.removeFalseStarts) {
      var before = kept.length;
      kept = kept.filter(function (v) { return v >= opts.falseStartMs; });
      nFalseStarts = before - kept.length;
    }

    if (opts.removeOutliers && kept.length >= 3) {
      var m = mean(kept), s = sd(kept);
      if (s !== null && s > 0) {
        lo = m - opts.sdMultiplier * s;
        hi = m + opts.sdMultiplier * s;
        var before2 = kept.length;
        kept = kept.filter(function (v) { return v >= lo && v <= hi; });
        nOutliers = before2 - kept.length;
      }
    }

    return { kept: kept, nFalseStarts: nFalseStarts, nOutliers: nOutliers, lo: lo, hi: hi };
  }

  function toRs(rts) {
    var out = [];
    for (var i = 0; i < rts.length; i++) if (rts[i] > 0) out.push(1000 / rts[i]);
    return out;
  }

  /* ---------------- one block of trials ---------------- */

  /*
   * `hitRts` are the raw reaction times of the HITS in this block — misses
   * contribute no RT, and late responses (RT > hit window) are excluded because
   * they are misses. Miss and lapse counts are passed in separately.
   */
  function blockMetrics(hitRts, opts) {
    var sorted = hitRts.slice().sort(function (a, b) { return a - b; });
    var corr = applyCorrection(hitRts, opts);
    var cSorted = corr.kept.slice().sort(function (a, b) { return a - b; });

    var rs = toRs(hitRts);
    var rsSorted = rs.slice().sort(function (a, b) { return a - b; });
    var cRs = toRs(corr.kept);
    var cRsSorted = cRs.slice().sort(function (a, b) { return a - b; });

    /*
     * Interpercentile range: the spread between the two decile means, as an
     * absolute distance. It is a within-test stability measure — a participant
     * whose best and worst tenths sit far apart is performing erratically even
     * when the average looks intact — and it is what the normative tables
     * report, so the app has to produce a comparable number.
     */
    function ipr(a, b) { return a === null || b === null ? null : Math.abs(a - b); }

    var fastRt = decileMean(sorted, 'low');
    var slowRt = decileMean(sorted, 'high');
    var cFastRt = decileMean(cSorted, 'low');
    var cSlowRt = decileMean(cSorted, 'high');
    var fastRs = decileMean(rsSorted, 'high');
    var slowRs = decileMean(rsSorted, 'low');
    var cFastRs = decileMean(cRsSorted, 'high');
    var cSlowRs = decileMean(cRsSorted, 'low');

    return {
      n: hitRts.length,
      nCorrected: corr.kept.length,
      nFalseStartsRemoved: corr.nFalseStarts,
      nOutliersRemoved: corr.nOutliers,
      outlierLoMs: corr.lo,
      outlierHiMs: corr.hi,
      nRsUndefined: hitRts.length - rs.length,

      /* reaction time, raw */
      avgRt: mean(hitRts),
      medianRt: median(hitRts),
      sdRt: sd(hitRts),
      fastest10Rt: fastRt,
      slowest10Rt: slowRt,
      iprRt: ipr(slowRt, fastRt),

      /* reaction time, corrected */
      corrAvgRt: mean(corr.kept),
      corrMedianRt: median(corr.kept),
      corrSdRt: sd(corr.kept),
      corrFastest10Rt: cFastRt,
      corrSlowest10Rt: cSlowRt,
      corrIprRt: ipr(cSlowRt, cFastRt),

      /* reaction speed, raw — fastest is the HIGH tail */
      avgRs: mean(rs),
      medianRs: median(rs),
      sdRs: sd(rs),
      fastest10Rs: fastRs,
      slowest10Rs: slowRs,
      iprRs: ipr(fastRs, slowRs),

      /* reaction speed, corrected */
      corrAvgRs: mean(cRs),
      corrMedianRs: median(cRs),
      corrSdRs: sd(cRs),
      corrFastest10Rs: cFastRs,
      corrSlowest10Rs: cSlowRs,
      corrIprRs: ipr(cFastRs, cSlowRs)
    };
  }

  /* ---------------- error profiles ---------------- */

  /*
   * Runs of consecutive misses, binned by run length. The frequency reported is
   * the number of RUNS in each band, not the number of missed trials.
   *
   * A run of `missCriterion` or more terminates the test, so ep7plus is
   * normally 0 or 1.
   */
  function errorProfiles(outcomes, missCriterion) {
    var runs = [], cur = 0;
    for (var i = 0; i < outcomes.length; i++) {
      if (outcomes[i] === 'miss') cur += 1;
      else if (cur) { runs.push(cur); cur = 0; }
    }
    if (cur) runs.push(cur);

    var ep = { ep1_2: 0, ep3_6: 0, ep7plus: 0, runs: runs, longestRun: 0 };
    for (var j = 0; j < runs.length; j++) {
      var r = runs[j];
      if (r > ep.longestRun) ep.longestRun = r;
      if (r <= 2) ep.ep1_2 += 1;
      else if (r <= 6) ep.ep3_6 += 1;
      else ep.ep7plus += 1;
    }
    // A criterion of 0 or less means 'no early termination' (PVT runs to time).
    ep.criterionReached = missCriterion > 0 && ep.longestRun >= missCriterion;
    return ep;
  }

  /* ---------------- velocity and acceleration ---------------- */

  /*
   * Rate of change of reaction speed across successive minutes, and the rate of
   * change of that. Positive velocity means the participant is speeding up.
   *
   * Both are undefined for the first minute (and acceleration for the second),
   * reported as null rather than zero so the gap is visible in the data.
   */
  function dynamics(perMinute, key) {
    var series = perMinute.map(function (m) { return m[key]; });
    var velocity = [], acceleration = [];

    for (var i = 0; i < series.length; i++) {
      velocity.push(i === 0 || series[i] === null || series[i - 1] === null
        ? null
        : series[i] - series[i - 1]);
    }
    for (var j = 0; j < velocity.length; j++) {
      acceleration.push(j < 2 || velocity[j] === null || velocity[j - 1] === null
        ? null
        : velocity[j] - velocity[j - 1]);
    }
    return { velocity: velocity, acceleration: acceleration, slope: olsSlope(series) };
  }

  /* Least-squares slope of a series against its index; null if under two points. */
  function olsSlope(series) {
    var xs = [], ys = [];
    for (var i = 0; i < series.length; i++) {
      if (series[i] !== null && series[i] !== undefined) { xs.push(i); ys.push(series[i]); }
    }
    if (xs.length < 2) return null;
    var mx = mean(xs), my = mean(ys), num = 0, den = 0;
    for (var j = 0; j < xs.length; j++) {
      num += (xs[j] - mx) * (ys[j] - my);
      den += (xs[j] - mx) * (xs[j] - mx);
    }
    return den === 0 ? null : num / den;
  }

  /* ---------------- stimulus schedule ---------------- */

  /*
   * mulberry32. A seeded generator rather than Math.random so a schedule can be
   * reproduced exactly from the seed recorded with the data — which matters when
   * a reviewer asks what the participant actually saw.
   */
  function makeRng(seed) {
    let a = seed >>> 0;
    return function () {
      a = (a + 0x6D2B79F5) | 0;
      let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  /* Fisher-Yates, unbiased, driven by the seeded generator. */
  function shuffle(arr, rng) {
    const a = arr.slice();
    for (let i = a.length - 1; i > 0; i--) {
      const j = Math.floor(rng() * (i + 1));
      const t = a[i]; a[i] = a[j]; a[j] = t;
    }
    return a;
  }

  /*
   * Builds the inter-stimulus intervals for a whole trial.
   *
   * PVT mode: the specified set (2, 4, 6, 8, 10 s) sums to exactly 30 s, so a
   * 30-second block holds precisely one of each interval. The schedule is
   * therefore a fresh random PERMUTATION per block: random in order, balanced in
   * composition, so every block delivers the same distribution of waits and no
   * block can happen to be all-short or all-long. `method` records which rule
   * was used, because if a custom set does not sum to the block length that
   * balance is impossible and the generator falls back to sampling with
   * replacement.
   *
   * Each interval PRECEDES its stimulus, so the participant always waits before
   * the first one rather than being surprised the instant the countdown ends.
   * A schedule of N intervals yields N stimuli — the last one gets a response
   * window as long as the longest interval in the configured set, so the
   * trial runs a little past the sum of the N intervals rather than cutting
   * the final stimulus off the instant it appears. That overshoot is a fixed
   * property of the settings, so it is identical for every participant
   * regardless of which interval the random schedule happened to draw last.
   */
  function buildSchedule(opts) {
    const maxMs = opts.maxMs;
    const isis = [];

    if (opts.mode !== 'pvt') {
      const n = Math.max(1, Math.ceil(maxMs / opts.isiMs));
      for (let i = 0; i < n; i++) isis.push(opts.isiMs);
      return finishSchedule(isis, maxMs, opts.blockMs, 'fixed', opts.seed, [opts.isiMs]);
    }

    const set = opts.isiSetMs.slice();
    const blockMs = opts.blockMs;
    const sum = set.reduce(function (a, b) { return a + b; }, 0);
    const balanced = sum === blockMs;
    const rng = makeRng(opts.seed);

    let total = 0;
    let guard = 0;
    while (total < maxMs && guard++ < 100000) {
      if (balanced) {
        const block = shuffle(set, rng);
        for (let i = 0; i < block.length && total < maxMs; i++) {
          isis.push(block[i]);
          total += block[i];
        }
      } else {
        let used = 0;
        while (used < blockMs && total < maxMs) {
          const v = set[Math.floor(rng() * set.length)];
          isis.push(v);
          used += v;
          total += v;
        }
      }
    }
    return finishSchedule(isis, maxMs, blockMs, balanced ? 'block_permutation' : 'sampled_with_replacement', opts.seed, set);
  }

  function finishSchedule(isis, maxMs, blockMs, method, seed, intervalSet) {
    // Keep only whole intervals that fit inside the ceiling.
    var kept = [];
    var total = 0;
    for (var i = 0; i < isis.length; i++) {
      if (total + isis[i] > maxMs) break;
      kept.push(isis[i]);
      total += isis[i];
    }

    /*
     * Interval k spans [start_k, start_k + kept[k]) and its stimulus fires at
     * the END of it — every kept interval gets a stimulus, including the last
     * one, so an M-interval schedule yields M stimuli, not M-1. (An earlier
     * version stopped one short, reserving the final interval as a trailing
     * response window with nothing to respond to: a 3-minute, 3000 ms BSRT
     * showed 59 flashes instead of 60.)
     *
     * A stimulus belongs to the block its PRECEDING interval started in, so
     * each block holds exactly one stimulus per interval in the set even though
     * the stimuli themselves sit at block-relative offsets.
     *
     * epochIsi is this stimulus's own response window — normally the interval
     * that leads into the NEXT stimulus, so a late response is still within
     * an active epoch until the next one begins. The last stimulus has no
     * next interval to borrow, so it gets the longest interval in the
     * CONFIGURED set instead of whatever happened to be drawn last. That
     * value is fixed by the settings, not by the random schedule, so trial
     * length stays identical for every participant regardless of seed — using
     * the actual last-drawn interval instead once made a 3-minute PVT run
     * anywhere from 602 to 610 seconds depending on which of 2/4/6/8/10 s
     * happened to be shuffled to the end.
     */
    var trailingMs = Math.max.apply(null, intervalSet);
    var onsets = [], blocks = [], minutes = [], isiBefore = [], epochIsi = [];
    var acc = 0;
    for (var k = 0; k < kept.length; k++) {
      blocks.push(Math.floor(acc / blockMs));
      // A minute is the same rule with a 60 s window. Both are taken from the
      // INTENDED schedule, never from a measured onset: presentation jitter of
      // a millisecond either way would otherwise flip a boundary stimulus into
      // the neighbouring bucket and unbalance the counts.
      minutes.push(Math.floor(acc / 60000));
      acc += kept[k];
      onsets.push(acc);
      isiBefore.push(kept[k]);
      epochIsi.push(k + 1 < kept.length ? kept[k + 1] : trailingMs);
    }

    return {
      isis: kept,
      onsets: onsets,
      blocks: blocks,
      minutes: minutes,
      isiBefore: isiBefore,
      epochIsi: epochIsi,
      method: method,
      seed: seed,
      blockMs: blockMs,
      nStimuli: onsets.length,
      leadInMs: kept.length ? kept[0] : 0,
      plannedDurationMs: onsets.length ? onsets[onsets.length - 1] + epochIsi[epochIsi.length - 1] : 0
    };
  }

  /* ---------------- response integrity ---------------- */

  /*
   * Continuous or repeated pressing is a way to fake alertness: hold the key
   * down, or tap constantly, and every stimulus gets a fast "response". It is
   * detectable because genuine responding produces roughly one press per
   * stimulus, spaced by the interval.
   *
   * Three independent signals, all reported so an experimenter can see WHY a
   * trial was flagged:
   *   - burstMax    the most presses inside any 1-second window
   *   - rapidPairs  consecutive presses closer together than rapidGapMs
   *   - extraPresses presses beyond the first in an epoch
   *
   * This flags a trial for human review. It never alters scoring, and it never
   * excludes data by itself — a startled double-tap is not cheating, and the
   * decision belongs to the experimenter.
   */
  function detectCheating(presses, nTrials, cfg) {
    var burstWindowMs = cfg.burstWindowMs == null ? 1000 : cfg.burstWindowMs;
    var burstLimit = cfg.burstLimit == null ? 5 : cfg.burstLimit;
    var rapidGapMs = cfg.rapidGapMs == null ? 200 : cfg.rapidGapMs;
    var extraRateLimit = cfg.extraRateLimit == null ? 0.2 : cfg.extraRateLimit;

    var times = presses.map(function (p) { return p.tMs; })
                       .sort(function (a, b) { return a - b; });

    var burstMax = 0;
    var lo = 0;
    for (var hi = 0; hi < times.length; hi++) {
      while (times[hi] - times[lo] > burstWindowMs) lo += 1;
      if (hi - lo + 1 > burstMax) burstMax = hi - lo + 1;
    }

    var rapidPairs = 0;
    for (var i = 1; i < times.length; i++) {
      if (times[i] - times[i - 1] < rapidGapMs) rapidPairs += 1;
    }

    var seen = {}, extraPresses = 0;
    for (var j = 0; j < presses.length; j++) {
      var k = presses[j].epochIndex;
      if (k === null || k === undefined) { extraPresses += 1; continue; }
      if (seen[k]) extraPresses += 1; else seen[k] = 1;
    }

    var extraRate = nTrials ? extraPresses / nTrials : 0;
    var reasons = [];
    if (burstMax >= burstLimit) {
      reasons.push(burstMax + ' presses within ' + burstWindowMs + ' ms');
    }
    if (nTrials >= 5 && extraRate > extraRateLimit) {
      reasons.push((extraRate * 100).toFixed(0) + '% more presses than stimuli');
    }
    if (times.length > 10 && rapidPairs / times.length > 0.25) {
      reasons.push(rapidPairs + ' presses under ' + rapidGapMs + ' ms apart');
    }

    return {
      totalPresses: times.length,
      extraPresses: extraPresses,
      extraRate: extraRate,
      burstMax: burstMax,
      rapidPairs: rapidPairs,
      suspected: reasons.length > 0,
      reasons: reasons,
      thresholds: {
        burstWindowMs: burstWindowMs, burstLimit: burstLimit,
        rapidGapMs: rapidGapMs, extraRateLimit: extraRateLimit
      }
    };
  }

  /* ---------------- main entry point ---------------- */

  /*
   * epochs: [{ index, onsetMs, rtMs }] where rtMs is the raw reaction time of
   *         whatever response occurred in that epoch (however late), or null.
   *
   * presses: every keypress in the trial as { tMs, epochIndex }, used only for
   *          the response-integrity check. Optional.
   *
   * cfg:    { hitWindowMs, isiMs, lapseMs, missCriterion,
   *           removeAnticipations, removeOutliers, anticipationMs, sdMultiplier }
   */
  function score(epochs, cfg, presses) {
    var opts = {
      removeFalseStarts: !!cfg.removeFalseStarts,
      removeOutliers: !!cfg.removeOutliers,
      falseStartMs: cfg.falseStartMs == null ? 100 : cfg.falseStartMs,
      sdMultiplier: cfg.sdMultiplier == null ? 2 : cfg.sdMultiplier
    };
    var hitWindow = cfg.hitWindowMs;

    /*
     * The minute comes from the schedule (see finishSchedule), so it is a
     * planned bucket rather than something re-derived from a measured onset.
     *
     * The fallback below is for epochs built without one. It bucket by where
     * the PRECEDING INTERVAL started, matching the schedule's rule: bucketing
     * on the onset alone pushes every boundary stimulus a minute late, because
     * onsets span (0, ceiling] rather than [0, ceiling). In a 2-minute,
     * 3000 ms BSRT that put the stimuli at 60 s and 120 s into minutes 2 and
     * 3, reading 19 / 20 / 1 across a two-minute test — a spurious third
     * minute holding one trial, whose mean, velocity and acceleration were all
     * meaningless. It is only a fallback because subtracting an intended
     * interval from a measured onset still lands within a millisecond of the
     * boundary, and jitter that small is enough to flip the bucket.
     */
    function minuteOf(e) {
      if (e.minute != null) return e.minute;
      var start = e.isiBeforeMs == null ? e.onsetMs : e.onsetMs - e.isiBeforeMs;
      if (start < 0) start = 0;
      return Math.floor(start / 60000);
    }

    /* --- classify every epoch --- */
    var trials = epochs.map(function (e) {
      var rt = (e.rtMs === undefined ? null : e.rtMs);
      var isHit = rt !== null && rt <= hitWindow;
      return {
        index: e.index,
        onsetMs: e.onsetMs,
        minute: minuteOf(e),
        block: e.block === undefined ? null : e.block,
        epochIsiMs: e.epochIsiMs === undefined ? null : e.epochIsiMs,
        isiBeforeMs: e.isiBeforeMs === undefined ? null : e.isiBeforeMs,
        rtMs: rt,                                   // raw, always preserved
        rsPerSec: rt !== null && rt > 0 ? 1000 / rt : null,
        outcome: isHit ? 'hit' : 'miss',
        lateResponse: rt !== null && !isHit ? 1 : 0, // responded, but too late
        lapse: isHit && rt > cfg.lapseMs ? 1 : 0,
        falseStart: rt !== null && rt < opts.falseStartMs ? 1 : 0
      };
    });

    var hits = trials.filter(function (t) { return t.outcome === 'hit'; });
    var misses = trials.filter(function (t) { return t.outcome === 'miss'; });
    var lapses = trials.filter(function (t) { return t.lapse; });
    var late = trials.filter(function (t) { return t.lateResponse; });
    var falseStarts = trials.filter(function (t) { return t.falseStart; });

    /* --- per minute --- */
    var lastMinute = trials.length ? trials[trials.length - 1].minute : -1;
    var perMinute = [];
    for (var m = 0; m <= lastMinute; m++) {
      var inMin = trials.filter(function (t) { return t.minute === m; });
      var minHits = inMin.filter(function (t) { return t.outcome === 'hit'; });
      var block = blockMetrics(minHits.map(function (t) { return t.rtMs; }), opts);
      block.minute = m + 1;                       // 1-based for reporting
      block.trials = inMin.length;
      block.hits = minHits.length;
      block.misses = inMin.length - minHits.length;
      block.lapses = inMin.filter(function (t) { return t.lapse; }).length;
      block.lateResponses = inMin.filter(function (t) { return t.lateResponse; }).length;
      block.hitRatio = inMin.length ? minHits.length / inMin.length : null;
      perMinute.push(block);
    }

    /* --- whole test --- */
    var totals = blockMetrics(hits.map(function (t) { return t.rtMs; }), opts);
    totals.trials = trials.length;
    totals.hits = hits.length;
    totals.misses = misses.length;
    totals.lapses = lapses.length;
    totals.lateResponses = late.length;
    totals.falseStarts = falseStarts.length;
    totals.hitRatio = trials.length ? hits.length / trials.length : null;

    var ep = errorProfiles(trials.map(function (t) { return t.outcome; }), cfg.missCriterion);

    return {
      config: {
        hitWindowMs: hitWindow,
        isiMs: cfg.isiMs,
        lapseMs: cfg.lapseMs,
        missCriterion: cfg.missCriterion,
        removeFalseStarts: opts.removeFalseStarts,
        removeOutliers: opts.removeOutliers,
        falseStartMs: opts.falseStartMs,
        sdMultiplier: opts.sdMultiplier
      },
      trials: trials,
      perMinute: perMinute,
      totals: totals,
      errorProfiles: ep,
      integrity: detectCheating(presses || [], trials.length, cfg),
      dynamicsRs: dynamics(perMinute, 'avgRs'),
      dynamicsRt: dynamics(perMinute, 'avgRt')
    };
  }

  /* ---------------- normative comparison ---------------- */

  /*
   * Recomputes a trial's summary the way the normative workbook computed its
   * reference values, over the first `minutes` minutes only.
   *
   * The app's own scoring is NOT reused here, because the two disagree in two
   * ways that would bias every z-score:
   *
   *   - The app treats any response slower than the hit window as a miss and
   *     drops it from the RT statistics. The norms keep every response that is
   *     not a timeout, however slow. Across the 291 reference sessions this
   *     alone moves the RT mean by a median of 4.8 ms.
   *   - The app sizes a decile with round(10% x n), the norms with ceil. Those
   *     pick different subset sizes in about one session in ten.
   *
   * So the participant is re-scored under the workbook's rules for the purpose
   * of comparison only. The app's own metrics, which follow the protocol the
   * user specified, are left untouched and remain what the rest of the results
   * screen and the raw export report.
   */
  function normativeSummary(trials, minutes, opts) {
    // The same minute buckets the per-minute table uses, so "the first three
    // minutes" means one thing across the whole app.
    var inWindow = trials.filter(function (t) {
      return t.minute != null ? t.minute < minutes : t.onsetMs <= minutes * 60000;
    });
    if (!inWindow.length) return null;

    // A response of any latency counts; only a timeout is a miss.
    var responded = inWindow.filter(function (t) { return t.rtMs !== null; });
    var misses = inWindow.length - responded.length;
    var fsMs = opts && opts.falseStartMs != null ? opts.falseStartMs : 100;
    var falseStarts = responded.filter(function (t) { return t.rtMs < fsMs; }).length;
    var lapseMs = opts && opts.lapseMs != null ? opts.lapseMs : 500;

    var valid = responded
      .filter(function (t) { return t.rtMs > fsMs; })
      .map(function (t) { return t.rtMs; })
      .sort(function (a, b) { return a - b; });

    var ep = errorProfiles(inWindow.map(function (t) {
      return t.rtMs === null ? 'miss' : 'hit';
    }), opts && opts.missCriterion != null ? opts.missCriterion : 0);

    var out = {
      minutes: minutes,
      trials: inWindow.length,
      hits: inWindow.length - misses,
      misses: misses,
      hitRatio: inWindow.length ? ((inWindow.length - misses) / inWindow.length) * 100 : null,
      falseStarts: falseStarts,
      lapses: valid.filter(function (v) { return v > lapseMs; }).length,
      ep12: ep.ep1_2,
      ep36: ep.ep3_6,
      ep7: ep.ep7plus,
      nValid: valid.length
    };

    if (!valid.length) return out;

    var rs = valid.map(function (v) { return 1000 / v; }).sort(function (a, b) { return a - b; });
    var fastRt = ceilDecile(valid, 'low');
    var slowRt = ceilDecile(valid, 'high');
    var fastRs = ceilDecile(rs, 'high');
    var slowRs = ceilDecile(rs, 'low');

    out.rtMean = mean(valid);
    out.rtMedian = median(valid);
    out.rtSd = sd(valid);
    out.rtFast10 = fastRt;
    out.rtSlow10 = slowRt;
    out.rtIpr = Math.abs(slowRt - fastRt);
    out.rsMean = mean(rs);
    out.rsMedian = median(rs);
    out.rsSd = sd(rs);
    out.rsFast10 = fastRs;
    out.rsSlow10 = slowRs;
    out.rsIpr = Math.abs(fastRs - slowRs);
    return out;
  }

  /* The workbook's decile size: ceil, not round. */
  function ceilDecile(sortedAsc, tail) {
    if (!sortedAsc.length) return null;
    var k = Math.max(1, Math.ceil(sortedAsc.length * 0.1));
    var slice = tail === 'low' ? sortedAsc.slice(0, k) : sortedAsc.slice(sortedAsc.length - k);
    return mean(slice);
  }

  /*
   * Position of one value against its reference cell.
   *
   * `dir` is +1 when a higher value is better and -1 when it is worse, so
   * `zWorse` is always signed the same way: positive means worse than the
   * reference mean, whatever the variable measures. The banding the user asked
   * for keys off that single number.
   *
   * An SD of exactly zero is common in the reference table — every control
   * session scored 20 hits out of 20 in the first minute, for instance. A
   * z-score is undefined there, so rather than dividing by zero and reporting
   * Infinity, the value is flagged as sitting outside a reference range that
   * had no spread at all, and z stays null.
   */
  function compareToNorm(value, norm) {
    if (value == null || !norm || norm.mean == null) return null;
    var deltaWorse = norm.dir >= 0 ? norm.mean - value : value - norm.mean;
    var out = {
      value: value,
      mean: norm.mean,
      sd: norm.sd,
      n: norm.n,
      dir: norm.dir,
      deltaWorse: deltaWorse,
      z: null,
      band: 'green',
      degenerate: false
    };
    if (norm.dir === 0) { out.band = 'neutral'; return out; }
    if (!norm.sd) {
      out.degenerate = true;
      // No spread to measure against: worse at all means outside the whole
      // reference sample, which is worth flagging, but never as a z-score.
      out.band = deltaWorse > 0 ? 'red' : 'green';
      return out;
    }
    out.z = deltaWorse / norm.sd;
    out.band = out.z > 2 ? 'red' : out.z > 1 ? 'orange' : 'green';
    return out;
  }

  /*
   * The whole normative comparison for one trial, or a refusal with a reason.
   *
   * Eligibility is deliberately strict. The reference sessions are BSRT runs
   * with a 3000 ms interval, which is 20 stimuli a minute; a PVT averages ten,
   * and a BSRT at some other interval delivers a different number again. Trial
   * counts, miss counts and lapse counts would then be measuring different
   * quantities, and reaction time itself moves with the interval. Rather than
   * emit a plausible-looking z-score against the wrong reference, the report
   * says which condition failed.
   */
  function normativeReport(o, norms) {
    var out = {
      available: false, reason: null, hour: o.hour, mode: o.mode, isiMs: o.isiMs,
      windowMinutes: null, testMinutes: null, truncated: false,
      belowProtocol: false, protocolMinutes: norms ? norms.protocolMinutes : null,
      hourMaxLength: null, provenance: null, protocols: [], mixedStudies: false,
      censored: 0,
      sessions: norms ? norms.sessions : null,
      participants: norms ? norms.data.participants : null,
      source: norms ? norms.data.source : null,
      n: null, summary: null, rows: []
    };
    if (!norms) { out.reason = 'no_norms'; return out; }
    if (o.mode !== 'bsrt') { out.reason = 'not_bsrt'; return out; }
    if (o.isiMs !== 3000) { out.reason = 'isi_mismatch'; return out; }
    if (!(o.hour >= 0 && o.hour < 24)) { out.reason = 'no_hour'; return out; }

    /*
     * How deep the reference goes depends on the hour. Only the four overnight
     * hours were ever tested past 8 minutes, so a 20-minute test at 14:00 is
     * compared over its first 8 minutes rather than against nothing at all;
     * windowMinutes says which window was actually used.
     */
    var hourMax = norms.maxLengthFor ? norms.maxLengthFor(o.hour) : norms.maxLength;
    if (!hourMax) { out.reason = 'no_hour'; return out; }
    out.hourMaxLength = hourMax;

    // Whole minutes actually completed, never more than the reference reaches.
    var completed = Math.floor((o.elapsedMs || 0) / 60000);
    var buckets = o.minuteBuckets == null ? completed : o.minuteBuckets;
    var window = Math.min(completed, buckets, hourMax);
    out.testMinutes = completed;
    if (window < 1) { out.reason = 'too_short'; return out; }
    out.windowMinutes = window;
    out.truncated = completed > window;

    /*
     * Which studies stand behind this particular cell, and how many sessions
     * had already stopped at sleep onset before reaching this length. Both
     * change what the comparison means, so both are carried out to the report
     * rather than being averaged away silently.
     */
    var prov = norms.provenance ? norms.provenance(o.hour, window) : null;
    if (prov) {
      out.provenance = prov;
      out.censored = prov[2] || 0;
      if (prov[0]) out.protocols.push(8);
      if (prov[1]) out.protocols.push(40);
      out.mixedStudies = prov[0] > 0 && prov[1] > 0;
    }
    // Shorter than the protocol the reference sessions were run under, so the
    // participants behind it were pacing themselves for a longer test.
    out.belowProtocol = out.protocols.some(function (p) { return window < p; });

    var summary = normativeSummary(o.trials, window, o);
    if (!summary) { out.reason = 'too_short'; return out; }

    out.summary = summary;
    out.n = norms.nFor(o.hour, window);
    out.available = true;

    for (var i = 0; i < norms.vars.length; i++) {
      var v = norms.vars[i];
      var norm = norms.lookup(o.hour, window, v.key);
      var cmp = compareToNorm(summary[v.key], norm);
      out.rows.push({
        key: v.key, label: v.label, section: v.section, unit: v.unit,
        dir: v.dir, comparison: cmp
      });
    }

    // A single headline count, so the panel can say how much stood out.
    out.nOrange = out.rows.filter(function (r) { return r.comparison && r.comparison.band === 'orange'; }).length;
    out.nRed = out.rows.filter(function (r) { return r.comparison && r.comparison.band === 'red'; }).length;
    return out;
  }

  return {
    score: score,
    normativeSummary: normativeSummary,
    normativeReport: normativeReport,
    compareToNorm: compareToNorm,
    detectCheating: detectCheating,
    buildSchedule: buildSchedule,
    makeRng: makeRng,
    shuffle: shuffle,
    blockMetrics: blockMetrics,
    errorProfiles: errorProfiles,
    applyCorrection: applyCorrection,
    dynamics: dynamics,
    olsSlope: olsSlope,
    stats: { mean: mean, median: median, sd: sd, decileMean: decileMean }
  };
});
