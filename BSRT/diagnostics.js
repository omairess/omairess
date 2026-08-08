'use strict';

/* Device and environment diagnostics.
 *
 * SHARED FILE — a byte-identical copy lives at BSRT-desktop/renderer/diagnostics.js.
 *
 * The desktop build knows what machine it is on because Electron tells it. A
 * browser has to find out by measuring, and some of what it finds is worse
 * news than the desktop ever has to report. That is the point: if a task is
 * going to be run on someone else's laptop for screening, the file should say
 * what that laptop was actually capable of.
 *
 *
 * ABOUT NETWORK LATENCY, WHICH IS THE THING PEOPLE EXPECT TO SEE HERE.
 *
 * It does not affect reaction times in this app, and reporting it as if it
 * did would be actively harmful — someone would exclude a participant for a
 * slow connection that had no bearing on their data.
 *
 * The whole task is client-side. Once the page has loaded, every stimulus is
 * drawn by the local compositor and every keypress is timed by the local
 * clock. Nothing crosses the network between the stimulus appearing and the
 * response being stamped. A participant on satellite internet and one on
 * fibre get identical timing, because the network is not in the loop.
 *
 * What the connection CAN do is fail to deliver the page, or deliver it
 * slowly, and a page still fetching assets can drop frames. That is why the
 * delivery figures below are recorded — as provenance for how the app reached
 * the participant, and to explain dropped frames in the first seconds — and
 * why they are labelled as provenance rather than as timing accuracy.
 *
 * The things that DO vary between remote devices, and are measured here:
 *   - the display's real refresh rate and stability (timing.js)
 *   - how coarse performance.now() is on that browser
 *   - how long input events take to reach the page — a Bluetooth keyboard can
 *     add tens of milliseconds that no software can recover
 *   - whether the window kept focus, and whether it ran full screen
 */

(function (root, factory) {
  var api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.BSRTDiagnostics = api;
})(typeof self !== 'undefined' ? self : this, function () {

  /*
   * The smallest difference performance.now() will actually report.
   *
   * Browsers deliberately coarsen this clock as a Spectre mitigation: 100 us
   * in Chromium, 1 ms in Firefox and Safari unless the page is
   * cross-origin-isolated. It is a hard floor on reaction-time resolution
   * that no amount of care in the task can get below, so it belongs in the
   * record next to the refresh rate.
   *
   * Measured rather than assumed, because the clamp depends on browser,
   * version and isolation state.
   */
  function timerResolutionMs(maxIterations) {
    var cap = maxIterations || 300000;
    var deltas = [];
    for (var i = 0; i < cap && deltas.length < 50; i++) {
      var a = performance.now();
      var b = performance.now();
      if (b > a) deltas.push(b - a);
    }
    if (!deltas.length) return null;
    deltas.sort(function (x, y) { return x - y; });
    return deltas[0];
  }

  /* How the page was delivered. Provenance, not timing — see the header. */
  function delivery() {
    var out = {
      protocol: null, host: null, local: null,
      loadMs: null, networkMs: null,
      connectionType: null, connectionRttMs: null, downlinkMbps: null, saveData: null
    };
    try {
      out.protocol = location.protocol.replace(':', '');
      out.host = location.hostname || null;
      // file:// and localhost never touched a network at all.
      out.local = out.protocol === 'file' ||
                  out.host === 'localhost' || out.host === '127.0.0.1';
    } catch (e) { /* not in a document */ }

    try {
      var nav = performance.getEntriesByType('navigation')[0];
      if (nav) {
        out.loadMs = nav.duration > 0 ? nav.duration : null;
        // Request to last byte: the part the connection is actually
        // responsible for.
        if (nav.responseEnd > 0 && nav.requestStart > 0) {
          out.networkMs = nav.responseEnd - nav.requestStart;
        }
      }
    } catch (e) { /* not supported */ }

    try {
      var c = navigator.connection;
      if (c) {
        out.connectionType = c.effectiveType || null;
        // Browser's own coarse estimate, rounded to 25 ms. Not a measurement
        // of anything this task depends on.
        out.connectionRttMs = typeof c.rtt === 'number' ? c.rtt : null;
        out.downlinkMbps = typeof c.downlink === 'number' ? c.downlink : null;
        out.saveData = typeof c.saveData === 'boolean' ? c.saveData : null;
      }
    } catch (e) { /* not supported */ }
    return out;
  }

  /* Everything obtainable without running a probe. */
  function snapshot() {
    var s = {
      userAgent: null, platform: null, language: null,
      screenW: null, screenH: null, windowW: null, windowH: null,
      devicePixelRatio: null, colorDepth: null,
      cores: null, memoryGb: null, maxTouchPoints: null,
      crossOriginIsolated: null, timerResolutionMs: null,
      delivery: delivery()
    };
    try {
      s.userAgent = navigator.userAgent || null;
      s.platform = navigator.platform || null;
      s.language = navigator.language || null;
      s.cores = navigator.hardwareConcurrency || null;
      s.memoryGb = typeof navigator.deviceMemory === 'number' ? navigator.deviceMemory : null;
      s.maxTouchPoints = typeof navigator.maxTouchPoints === 'number' ? navigator.maxTouchPoints : null;
    } catch (e) { /* ignore */ }
    try {
      s.screenW = screen.width; s.screenH = screen.height;
      s.colorDepth = screen.colorDepth;
      s.devicePixelRatio = window.devicePixelRatio;
      s.windowW = window.innerWidth; s.windowH = window.innerHeight;
    } catch (e) { /* ignore */ }
    try { s.crossOriginIsolated = !!window.crossOriginIsolated; } catch (e) { /* ignore */ }
    s.timerResolutionMs = timerResolutionMs();
    return s;
  }

  /* A short, readable identification of the browser and OS, for the setup line. */
  function describeBrowser(ua) {
    if (!ua) return 'unknown browser';
    var browser = 'unknown browser';
    var m;
    if ((m = ua.match(/Edg\/([\d.]+)/))) browser = 'Edge ' + m[1].split('.')[0];
    else if ((m = ua.match(/OPR\/([\d.]+)/))) browser = 'Opera ' + m[1].split('.')[0];
    else if ((m = ua.match(/Firefox\/([\d.]+)/))) browser = 'Firefox ' + m[1].split('.')[0];
    else if ((m = ua.match(/Chrome\/([\d.]+)/))) browser = 'Chrome ' + m[1].split('.')[0];
    else if ((m = ua.match(/Version\/([\d.]+).*Safari/))) browser = 'Safari ' + m[1].split('.')[0];

    var os = 'unknown OS';
    if (/Windows NT 10/.test(ua)) os = 'Windows 10/11';
    else if (/Windows/.test(ua)) os = 'Windows';
    else if (/Mac OS X/.test(ua)) os = 'macOS';
    else if (/Android/.test(ua)) os = 'Android';
    else if (/iPhone|iPad|iPod/.test(ua)) os = 'iOS';
    else if (/Linux/.test(ua)) os = 'Linux';
    return browser + ' on ' + os;
  }

  /*
   * A verdict aimed at someone deciding whether a remote device produced
   * usable data. Ordered worst-first, and every note says what the number
   * means for the reaction times rather than just quoting it.
   *
   * `cal` is the timing.js display calibration, `input` its input-probe
   * result; either may be null if that measurement was not made.
   */
  function screeningReport(env, cal, input, opts) {
    opts = opts || {};
    var notes = [];
    var blocking = 0;   // things that make the data hard to defend
    var caution = 0;    // things worth knowing but not disqualifying

    if (cal && cal.refreshHz) {
      if (cal.refreshHz < 50) {
        notes.push({ level: 'bad', key: 'diag.nLowHz',
                     vars: { hz: cal.refreshHz.toFixed(1), q: (cal.frameIntervalMs / 2).toFixed(1) } });
        blocking++;
      }
      if (cal.droppedRate != null && cal.droppedRate > 0.05) {
        notes.push({ level: 'bad', key: 'diag.nDropped',
                     vars: { pct: (cal.droppedRate * 100).toFixed(1) } });
        blocking++;
      } else if (cal.droppedRate != null && cal.droppedRate > 0.01) {
        notes.push({ level: 'warn', key: 'diag.nDroppedMild',
                     vars: { pct: (cal.droppedRate * 100).toFixed(1) } });
        caution++;
      }
      if (cal.madIntervalMs != null && cal.madIntervalMs > 2) {
        notes.push({ level: 'bad', key: 'diag.nJitter',
                     vars: { mad: cal.madIntervalMs.toFixed(2) } });
        blocking++;
      } else if (cal.madIntervalMs != null && cal.madIntervalMs > 1) {
        notes.push({ level: 'warn', key: 'diag.nJitterMild',
                     vars: { mad: cal.madIntervalMs.toFixed(2) } });
        caution++;
      }
    } else {
      notes.push({ level: 'warn', key: 'diag.nNoDisplay', vars: {} });
      caution++;
    }

    if (env && env.timerResolutionMs != null && env.timerResolutionMs >= 1) {
      notes.push({ level: 'warn', key: 'diag.nCoarseTimer',
                   vars: { res: env.timerResolutionMs.toFixed(2) } });
      caution++;
    }

    if (input && input.usable && input.medianDelayMs != null) {
      if (input.medianDelayMs > 50) {
        notes.push({ level: 'bad', key: 'diag.nInputSlow',
                     vars: { ms: input.medianDelayMs.toFixed(1) } });
        blocking++;
      } else if (input.medianDelayMs > 20) {
        notes.push({ level: 'warn', key: 'diag.nInputMild',
                     vars: { ms: input.medianDelayMs.toFixed(1) } });
        caution++;
      }
    } else if (input && !input.usable) {
      notes.push({ level: 'warn', key: 'diag.nInputUnknown', vars: {} });
      caution++;
    }

    if (env && env.maxTouchPoints > 0 && opts.touchUsed) {
      notes.push({ level: 'warn', key: 'diag.nTouch', vars: {} });
      caution++;
    }
    if (opts.blurCount > 0) {
      notes.push({ level: 'bad', key: 'diag.nBlur', vars: { n: opts.blurCount } });
      blocking++;
    }
    if (opts.fullscreen === false) {
      notes.push({ level: 'warn', key: 'diag.nWindowed', vars: {} });
      caution++;
    }

    var grade = blocking ? 'poor' : (caution ? 'fair' : 'good');
    return { grade: grade, blocking: blocking, caution: caution, notes: notes };
  }

  return {
    timerResolutionMs: timerResolutionMs,
    delivery: delivery,
    snapshot: snapshot,
    describeBrowser: describeBrowser,
    screeningReport: screeningReport
  };
});
