'use strict';

/* BSRT timing module.
 *
 * Three jobs:
 *   1. Measure the display's real frame interval (not the OS-declared one).
 *   2. Measure how long input events take to reach us after the OS stamped them.
 *   3. Give the task a frame-locked clock to schedule against.
 *
 * What this module can and cannot do is stated plainly in README.md. In short:
 * it removes quantisation and dispatch error, and it characterises what
 * remains. It cannot measure photon emission or physical key travel — no
 * software can. Those leave a roughly constant per-device offset that only a
 * photodiode can pin down.
 */

/* ---------------- robust statistics ---------------- */

function median(a) {
  if (!a.length) return null;
  const b = a.slice().sort((x, y) => x - y);
  const m = b.length >> 1;
  return b.length % 2 ? b[m] : (b[m - 1] + b[m]) / 2;
}

function quantile(a, q) {
  if (!a.length) return null;
  const b = a.slice().sort((x, y) => x - y);
  const pos = (b.length - 1) * q;
  const lo = Math.floor(pos);
  const hi = Math.ceil(pos);
  return lo === hi ? b[lo] : b[lo] + (b[hi] - b[lo]) * (pos - lo);
}

function mean(a) {
  if (!a.length) return null;
  return a.reduce((s, v) => s + v, 0) / a.length;
}

function sd(a) {
  if (a.length < 2) return null;
  const m = mean(a);
  return Math.sqrt(a.reduce((s, v) => s + (v - m) * (v - m), 0) / (a.length - 1));
}

/* Median absolute deviation, scaled to be comparable to an SD under normality. */
function mad(a) {
  if (a.length < 2) return null;
  const m = median(a);
  return 1.4826 * median(a.map((v) => Math.abs(v - m)));
}

/* ---------------- display calibration ---------------- */

/*
 * Runs a bare requestAnimationFrame loop and times every frame boundary.
 *
 * The rAF timestamp is the moment the browser began producing that frame, on
 * the same clock as performance.now(). Successive timestamps therefore give
 * the true compositor cadence, which is what the task schedules against.
 *
 * `warmupFrames` are discarded: the first frames after a loop starts are
 * routinely long while the compositor spins up, and including them would
 * inflate the jitter estimate and fake a dropped frame.
 */
function calibrateDisplay(nFrames, warmupFrames) {
  nFrames = nFrames || 240;
  warmupFrames = warmupFrames == null ? 10 : warmupFrames;

  return new Promise((resolve) => {
    const stamps = [];
    let seen = 0;

    function step(ts) {
      seen += 1;
      if (seen > warmupFrames) stamps.push(ts);
      if (seen < nFrames + warmupFrames) {
        requestAnimationFrame(step);
      } else {
        resolve(analyseFrames(stamps));
      }
    }
    requestAnimationFrame(step);
  });
}

function analyseFrames(stamps) {
  const intervals = [];
  for (let i = 1; i < stamps.length; i++) intervals.push(stamps[i] - stamps[i - 1]);

  const med = median(intervals);

  // A "dropped" frame shows up as an interval that is a whole multiple of the
  // nominal one. Count the missing frames, not merely the long intervals.
  let dropped = 0;
  const clean = [];
  for (const dt of intervals) {
    const mult = med ? Math.round(dt / med) : 1;
    if (mult >= 2) dropped += mult - 1;
    else clean.push(dt);
  }

  return {
    nFrames: stamps.length,
    nIntervals: intervals.length,
    frameIntervalMs: med,
    refreshHz: med ? 1000 / med : null,
    meanIntervalMs: mean(clean),
    sdIntervalMs: sd(clean),
    madIntervalMs: mad(clean),
    p05IntervalMs: quantile(clean, 0.05),
    p95IntervalMs: quantile(clean, 0.95),
    minIntervalMs: clean.length ? Math.min.apply(null, clean) : null,
    maxIntervalMs: clean.length ? Math.max.apply(null, clean) : null,
    droppedFrames: dropped,
    droppedRate: intervals.length ? dropped / intervals.length : null
  };
}

/*
 * Turns calibration output into a plain-language verdict, so an experimenter
 * who is not going to read an IQR still finds out their display is unstable.
 */
function gradeCalibration(cal, nominalHz) {
  const problems = [];

  if (!cal.frameIntervalMs) {
    return { grade: 'unusable', problems: ['No frames were timed.'] };
  }
  if (cal.droppedRate > 0.01) {
    problems.push(
      'Dropped ' + cal.droppedFrames + ' frames (' +
      (cal.droppedRate * 100).toFixed(1) + '%). Close other applications.'
    );
  }
  if (cal.madIntervalMs != null && cal.madIntervalMs > 1.0) {
    problems.push(
      'Frame timing is unstable (MAD ' + cal.madIntervalMs.toFixed(2) +
      ' ms). Reaction times will be noisier than usual.'
    );
  }
  if (nominalHz && cal.refreshHz) {
    const diff = Math.abs(cal.refreshHz - nominalHz);
    if (diff / nominalHz > 0.05) {
      problems.push(
        'Measured ' + cal.refreshHz.toFixed(1) + ' Hz but the system reports ' +
        nominalHz + ' Hz. The compositor may be limiting the frame rate.'
      );
    }
  }
  if (cal.refreshHz && cal.refreshHz < 50) {
    problems.push(
      'Refresh rate is only ' + cal.refreshHz.toFixed(1) +
      ' Hz, which coarsens stimulus onset resolution.'
    );
  }

  let grade = 'good';
  if (problems.length === 1) grade = 'fair';
  if (problems.length > 1) grade = 'poor';
  return { grade: grade, problems: problems };
}

/* ---------------- input dispatch calibration ---------------- */

/*
 * Chromium stamps input events with the time the *platform* delivered them,
 * on the performance.now() timebase. Reading event.timeStamp instead of
 * calling performance.now() inside the handler therefore skips whatever
 * queuing delay the event suffered inside the renderer.
 *
 * This probe measures that delay so we can (a) report it and (b) verify the
 * assumption. If a future Chromium changed the timebase, the deltas would go
 * wild and `usable` would come back false, at which point the task falls back
 * to handler-time and says so in the data rather than reporting nonsense.
 */
function makeInputProbe() {
  const deltas = [];
  return {
    record(evt) {
      const now = performance.now();
      if (typeof evt.timeStamp === 'number' && evt.timeStamp > 0) {
        deltas.push(now - evt.timeStamp);
      }
      return deltas.length;
    },
    result() {
      if (deltas.length < 3) return { usable: false, n: deltas.length, reason: 'too few samples' };
      const med = median(deltas);
      // Plausible dispatch delay: non-negative and under ~150 ms. Anything
      // outside that means event.timeStamp is not what we think it is.
      const usable = med >= -1 && med < 150;
      return {
        usable: usable,
        n: deltas.length,
        medianDelayMs: med,
        madDelayMs: mad(deltas),
        minDelayMs: Math.min.apply(null, deltas),
        maxDelayMs: Math.max.apply(null, deltas),
        reason: usable ? null : 'event.timeStamp is not on the performance clock'
      };
    }
  };
}

/*
 * The single place that decides when an input happened.
 * Falls back to handler time if the probe said the event clock is untrustworthy.
 */
function eventTime(evt, useEventStamp) {
  if (useEventStamp && typeof evt.timeStamp === 'number' && evt.timeStamp > 0) {
    return evt.timeStamp;
  }
  return performance.now();
}

window.BSRTTiming = {
  calibrateDisplay,
  analyseFrames,
  gradeCalibration,
  makeInputProbe,
  eventTime,
  stats: { median, mean, sd, mad, quantile }
};
