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
      fastest10Rt: decileMean(sorted, 'low'),
      slowest10Rt: decileMean(sorted, 'high'),

      /* reaction time, corrected */
      corrAvgRt: mean(corr.kept),
      corrMedianRt: median(corr.kept),
      corrSdRt: sd(corr.kept),
      corrFastest10Rt: decileMean(cSorted, 'low'),
      corrSlowest10Rt: decileMean(cSorted, 'high'),

      /* reaction speed, raw — fastest is the HIGH tail */
      avgRs: mean(rs),
      medianRs: median(rs),
      sdRs: sd(rs),
      fastest10Rs: decileMean(rsSorted, 'high'),
      slowest10Rs: decileMean(rsSorted, 'low'),

      /* reaction speed, corrected */
      corrAvgRs: mean(cRs),
      corrMedianRs: median(cRs),
      corrSdRs: sd(cRs),
      corrFastest10Rs: decileMean(cRsSorted, 'high'),
      corrSlowest10Rs: decileMean(cRsSorted, 'low')
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
   * The final interval carries no stimulus of its own — it is the response
   * window after the last one — so a schedule of N intervals yields N-1 stimuli
   * and lasts exactly the sum of its intervals. Total duration is therefore
   * unchanged, and identical for every participant, because each block sums to
   * the block length whatever order it was shuffled into.
   */
  function buildSchedule(opts) {
    const maxMs = opts.maxMs;
    const isis = [];

    if (opts.mode !== 'pvt') {
      const n = Math.max(1, Math.ceil(maxMs / opts.isiMs));
      for (let i = 0; i < n; i++) isis.push(opts.isiMs);
      return finishSchedule(isis, maxMs, opts.blockMs, 'fixed', opts.seed);
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
    return finishSchedule(isis, maxMs, blockMs, balanced ? 'block_permutation' : 'sampled_with_replacement', opts.seed);
  }

  function finishSchedule(isis, maxMs, blockMs, method, seed) {
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
     * the END of it. The last interval has no stimulus: it is the response
     * window for the one before, which is what keeps the trial from ending the
     * instant the final stimulus appears.
     *
     * A stimulus belongs to the block its PRECEDING interval started in, so
     * each block holds exactly one stimulus per interval in the set even though
     * the stimuli themselves sit at block-relative offsets.
     */
    var onsets = [], blocks = [], isiBefore = [], epochIsi = [];
    var acc = 0;
    for (var k = 0; k + 1 < kept.length; k++) {
      blocks.push(Math.floor(acc / blockMs));
      acc += kept[k];
      onsets.push(acc);
      isiBefore.push(kept[k]);
      epochIsi.push(kept[k + 1]);
    }

    return {
      isis: kept,
      onsets: onsets,
      blocks: blocks,
      isiBefore: isiBefore,
      epochIsi: epochIsi,
      method: method,
      seed: seed,
      blockMs: blockMs,
      nStimuli: onsets.length,
      leadInMs: kept.length ? kept[0] : 0,
      plannedDurationMs: total
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

    /* --- classify every epoch --- */
    var trials = epochs.map(function (e) {
      var rt = (e.rtMs === undefined ? null : e.rtMs);
      var isHit = rt !== null && rt <= hitWindow;
      return {
        index: e.index,
        onsetMs: e.onsetMs,
        minute: Math.floor(e.onsetMs / 60000),
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

  return {
    score: score,
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
