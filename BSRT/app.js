'use strict';

/* BSRT — Behavioral Sleep Resistance Task (browser build).
 *
 * Trial structure (see README):
 *   - a stimulus starts every isiMs, measured from the previous ONSET
 *   - it stays lit for up to stimMs, or until the participant responds
 *   - RT <= stimMs        -> HIT   (a slow hit, > lapseMs, is a LAPSE)
 *   - stimMs < RT < isiMs -> MISS, but the raw RT is still recorded
 *   - no response         -> MISS
 *   - the epoch clock never depends on behaviour
 *
 * All derived metrics come from scoring.js, shared byte-for-byte with the
 * desktop build so the two cannot drift apart.
 */

var S = window.BSRTScoring;

var LS_KEY = 'bsrt.sessions.v1';

/* ---------------- state ---------------- */

var cfg = null;
var meta = null;
var epochs = [];
var startTime = 0;
var startWall = null;
var currentEpoch = null;
var epochTimer = null;
var stimTimer = null;
var consecutiveMisses = 0;
var missRunStartIndex = -1;
var running = false;
var wakeLock = null;
var blurCount = 0;
var lastResult = null;

/* ---------------- helpers ---------------- */

function $(id) { return document.getElementById(id); }

function show(id) {
  var all = document.querySelectorAll('.screen');
  for (var i = 0; i < all.length; i++) all[i].classList.remove('active');
  $(id).classList.add('active');
}

function intVal(id, fb) { var n = parseInt($(id).value, 10); return isFinite(n) ? n : fb; }
function numVal(id, fb) { var n = parseFloat($(id).value); return isFinite(n) ? n : fb; }
function txtVal(id) { return $(id).value.trim(); }

function fmtClock(ms) {
  if (ms == null) return '—';
  var t = Math.round(ms / 1000), m = Math.floor(t / 60), s = t % 60;
  return m + ':' + (s < 10 ? '0' : '') + s;
}

function n1(v) { return v == null ? '—' : v.toFixed(1); }
function n2(v) { return v == null ? '—' : v.toFixed(2); }
function n3(v) { return v == null ? '—' : v.toFixed(3); }
function round(v, dp) { return v == null ? null : Number(v.toFixed(dp)); }

/* ---------------- CSV ---------------- */

function csvCell(v) {
  if (v === null || v === undefined) return '';
  var s = String(v);
  return /[",\n\r]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
}

function toCsv(header, rows) {
  var out = [header.map(csvCell).join(',')];
  for (var i = 0; i < rows.length; i++) out.push(rows[i].map(csvCell).join(','));
  return out.join('\r\n') + '\r\n';
}

function download(filename, text) {
  var blob = new Blob([text], { type: 'text/csv;charset=utf-8' });
  var url = URL.createObjectURL(blob);
  var a = document.createElement('a');
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  setTimeout(function () { URL.revokeObjectURL(url); }, 1000);
}

function safeName(s) { return String(s || 'NA').replace(/[^A-Za-z0-9_-]/g, '_'); }

function fileStem(r) {
  return 'bsrt_' + safeName(r.participant.participantId) + '_' +
         safeName(r.participant.sessionLabel) + '_t' + r.participant.trialNumber;
}

/* ---------------- storage ---------------- */

function loadSessions() {
  try {
    var a = JSON.parse(localStorage.getItem(LS_KEY) || '[]');
    return Array.isArray(a) ? a : [];
  } catch (e) { return []; }
}

function saveSession(r) {
  var all = loadSessions();
  all.push(r);
  try {
    localStorage.setItem(LS_KEY, JSON.stringify(all));
  } catch (e) {
    alert('Could not save this trial (storage may be full). Export it now.');
  }
  refreshCount();
}

function refreshCount() { $('sessionCount').textContent = loadSessions().length; }

/* ---------------- wake lock ---------------- */

function requestWakeLock() {
  if (!navigator.wakeLock || !navigator.wakeLock.request) return;
  navigator.wakeLock.request('screen').then(function (l) { wakeLock = l; }).catch(function () {});
}

function releaseWakeLock() {
  if (wakeLock) { try { wakeLock.release(); } catch (e) {} wakeLock = null; }
}

/* ---------------- trial ---------------- */

function validate() {
  var isi = intVal('isiMs', 3000), stim = intVal('stimMs', 1000);
  if (isi < 200) return 'Stimulus interval must be at least 200 ms.';
  if (stim >= isi) return 'Stimulus duration must be shorter than the stimulus interval.';
  if (intVal('missCriterion', 7) < 1) return 'The miss criterion must be at least 1.';
  if (intVal('maxMinutes', 40) < 1) return 'Maximum duration must be at least 1 minute.';
  if (intVal('lapseMs', 500) >= stim) return 'The lapse threshold must be below the hit window (stimulus duration).';
  return null;
}

function beginTrial() {
  var err = validate();
  if (err) { $('setupError').textContent = err; $('setupError').hidden = false; return; }
  $('setupError').hidden = true;

  var corr = $('correction').value; // none | anticipations | outliers | both
  var now = new Date();

  cfg = {
    isiMs: intVal('isiMs', 3000),
    stimMs: intVal('stimMs', 1000),
    hitWindowMs: intVal('stimMs', 1000),   // RT <= stimulus duration counts as a hit
    lapseMs: intVal('lapseMs', 500),
    missCriterion: intVal('missCriterion', 7),
    maxMinutes: intVal('maxMinutes', 40),
    maxMs: intVal('maxMinutes', 40) * 60000,
    correction: corr,
    removeAnticipations: corr === 'anticipations' || corr === 'both',
    removeOutliers: corr === 'outliers' || corr === 'both',
    anticipationMs: numVal('anticipationMs', 100),
    sdMultiplier: numVal('sdMultiplier', 2)
  };

  meta = {
    runId: 'bsrt-' + Date.now(),
    date: now.toISOString().slice(0, 10),
    time: now.toTimeString().slice(0, 8),
    participantId: txtVal('participantId') || 'NA',
    name: txtVal('pName'),
    address: txtVal('pAddress'),
    birthDate: txtVal('pBirth'),
    education: txtVal('pEducation'),
    sessionLabel: txtVal('sessionLabel') || 'NA',
    trialNumber: intVal('trialNumber', 1)
  };

  epochs = [];
  currentEpoch = null;
  consecutiveMisses = 0;
  missRunStartIndex = -1;
  blurCount = 0;
  running = false;

  requestWakeLock();
  if ($('useFullscreen').checked && document.documentElement.requestFullscreen) {
    document.documentElement.requestFullscreen().catch(function () {});
  }
  countdown();
}

function countdown() {
  show('screen-countdown');
  var n = 3;
  $('countdownNum').textContent = n;
  var iv = setInterval(function () {
    n -= 1;
    if (n <= 0) { clearInterval(iv); startTask(); }
    else $('countdownNum').textContent = n;
  }, 1000);
}

function startTask() {
  show('screen-task');
  ledOff();
  running = true;
  startWall = new Date().toISOString();
  startTime = performance.now();
  runEpoch(0);
}

/* Absolute targets, so presentation does not drift across a long trial. */
function scheduleEpoch(n) {
  var delay = (startTime + n * cfg.isiMs) - performance.now();
  epochTimer = setTimeout(function () { runEpoch(n); }, delay > 0 ? delay : 0);
}

function runEpoch(n) {
  finalizeEpoch();
  if (!running) return;
  if (n * cfg.isiMs >= cfg.maxMs) { endTask('max_duration'); return; }

  currentEpoch = {
    index: n,
    onsetMs: performance.now() - startTime,
    onsetPerf: performance.now(),
    responded: false,
    rtMs: null,
    extra: 0
  };
  epochs.push(currentEpoch);

  ledOn();
  stimTimer = setTimeout(ledOff, cfg.stimMs);
  scheduleEpoch(n + 1);
}

/*
 * An epoch counts toward the sleep-onset criterion unless it was a HIT.
 * A response slower than the hit window is a MISS, so it does NOT reset the
 * consecutive-miss run even though a key was pressed.
 */
function finalizeEpoch() {
  if (!currentEpoch) return;
  var ep = currentEpoch;
  currentEpoch = null;

  var isHit = ep.rtMs !== null && ep.rtMs <= cfg.hitWindowMs;
  if (isHit) { consecutiveMisses = 0; missRunStartIndex = -1; return; }

  if (consecutiveMisses === 0) missRunStartIndex = ep.index;
  consecutiveMisses += 1;
  if (consecutiveMisses >= cfg.missCriterion) endTask('sleep_onset');
}

/*
 * Every response is recorded with its raw RT, however late — the raw data must
 * survive for downstream analysis. Whether it counts as a hit is decided at
 * scoring time, not here.
 */
function handleResponse() {
  if (!running || !currentEpoch) return;
  if (currentEpoch.responded) { currentEpoch.extra += 1; return; }

  currentEpoch.responded = true;
  currentEpoch.rtMs = performance.now() - currentEpoch.onsetPerf;

  // Responding extinguishes the light immediately, but never shortens the
  // epoch: the next stimulus still starts a fixed interval after this onset.
  clearTimeout(stimTimer);
  ledOff();
}

function ledOn() { $('led').classList.add('on'); }
function ledOff() { $('led').classList.remove('on'); }

function endTask(reason) {
  running = false;
  clearTimeout(epochTimer);
  clearTimeout(stimTimer);
  ledOff();
  releaseWakeLock();
  if (document.fullscreenElement && document.exitFullscreen) document.exitFullscreen().catch(function () {});

  lastResult = buildResult(reason);
  saveSession(lastResult);
  renderResult(lastResult);
  show('screen-results');
}

/* ---------------- result assembly ---------------- */

function buildResult(reason) {
  var scored = S.score(epochs.map(function (e) {
    return { index: e.index, onsetMs: e.onsetMs, rtMs: e.rtMs };
  }), cfg);

  var sleepOnsetMs = null, sleepOnsetCriterionMs = null;
  if (reason === 'sleep_onset' && missRunStartIndex >= 0) {
    var first = epochs[missRunStartIndex];
    var last = epochs[missRunStartIndex + cfg.missCriterion - 1];
    sleepOnsetMs = first ? first.onsetMs : missRunStartIndex * cfg.isiMs;
    sleepOnsetCriterionMs = last ? last.onsetMs + cfg.isiMs
                                 : sleepOnsetMs + cfg.missCriterion * cfg.isiMs;
  }

  return {
    runId: meta.runId,
    startedAt: startWall,
    participant: meta,
    config: cfg,
    endReason: reason,
    sleptBeforeMax: reason === 'sleep_onset',
    sleepOnsetMs: sleepOnsetMs,
    sleepOnsetCriterionMs: sleepOnsetCriterionMs,
    elapsedMs: epochs.length ? epochs[epochs.length - 1].onsetMs + cfg.isiMs : 0,
    extraResponses: epochs.reduce(function (n, e) { return n + e.extra; }, 0),
    extras: epochs.map(function (e) { return e.extra; }),
    blurCount: blurCount,
    scored: scored
  };
}

/* ---------------- results rendering ---------------- */

function renderResult(r) {
  var p = r.participant, t = r.scored.totals, e = r.scored.errorProfiles;

  $('resultWho').textContent =
    (p.name || p.participantId) + ' · ' + p.sessionLabel + ' · trial ' + p.trialNumber +
    ' · ' + p.date + ' ' + p.time;

  if (r.sleptBeforeMax) {
    $('mLatency').textContent = fmtClock(r.sleepOnsetMs);
    $('mLatencyFoot').textContent =
      'to the first of ' + r.config.missCriterion + ' consecutive misses (criterion confirmed at ' +
      fmtClock(r.sleepOnsetCriterionMs) + ')';
    $('mOutcome').textContent = 'Sleep onset was scored before the ' + r.config.maxMinutes + '-minute ceiling.';
  } else {
    $('mLatency').textContent = '> ' + fmtClock(r.config.maxMs);
    $('mLatencyFoot').textContent = 'no sleep onset scored — censored at the maximum duration';
    $('mOutcome').textContent = r.endReason === 'aborted'
      ? 'The trial was ended early by the experimenter, so this is not a valid latency.'
      : 'The participant stayed awake for the full ' + r.config.maxMinutes + ' minutes.';
  }

  $('oTrials').textContent = t.trials;
  $('oHitRatio').textContent = t.hitRatio == null ? '—' : (t.hitRatio * 100).toFixed(1) + '%';
  $('oHits').textContent = t.hits;
  $('oMisses').textContent = t.misses + (t.lateResponses ? ' (' + t.lateResponses + ' late responses)' : '');
  $('oLapses').textContent = t.lapses + ' (> ' + r.config.lapseMs + ' ms)';
  $('oEp12').textContent = e.ep1_2;
  $('oEp36').textContent = e.ep3_6;
  $('oEp7').textContent = e.ep7plus;

  $('rtAvg').textContent = n1(t.avgRt);
  $('rtMedian').textContent = n1(t.medianRt);
  $('rtSd').textContent = n1(t.sdRt);
  $('rtFast').textContent = n1(t.fastest10Rt);
  $('rtSlow').textContent = n1(t.slowest10Rt);
  $('rtcAvg').textContent = n1(t.corrAvgRt);
  $('rtcMedian').textContent = n1(t.corrMedianRt);
  $('rtcFast').textContent = n1(t.corrFastest10Rt);
  $('rtcSlow').textContent = n1(t.corrSlowest10Rt);

  $('rsAvg').textContent = n3(t.avgRs);
  $('rsMedian').textContent = n3(t.medianRs);
  $('rsSd').textContent = n3(t.sdRs);
  $('rsFast').textContent = n3(t.fastest10Rs);
  $('rsSlow').textContent = n3(t.slowest10Rs);
  $('rscAvg').textContent = n3(t.corrAvgRs);
  $('rscMedian').textContent = n3(t.corrMedianRs);
  $('rscFast').textContent = n3(t.corrFastest10Rs);
  $('rscSlow').textContent = n3(t.corrSlowest10Rs);

  $('corrNote').textContent =
    r.config.correction === 'none'
      ? 'No correction applied — the corrected columns equal the raw ones.'
      : 'Correction: ' + describeCorrection(r.config) + '. Removed ' +
        t.nAnticipationsRemoved + ' anticipations and ' + t.nOutliersRemoved +
        ' outliers from ' + t.n + ' hits.';

  renderPerMinute(r);
}

function describeCorrection(c) {
  if (c.correction === 'anticipations') return 'RT < ' + c.anticipationMs + ' ms removed';
  if (c.correction === 'outliers') return 'RT beyond ' + c.sdMultiplier + ' SD removed';
  if (c.correction === 'both') return 'RT < ' + c.anticipationMs + ' ms removed, then RT beyond ' + c.sdMultiplier + ' SD';
  return 'none';
}

function renderPerMinute(r) {
  var pm = r.scored.perMinute;
  var vel = r.scored.dynamicsRs.velocity;
  var acc = r.scored.dynamicsRs.acceleration;
  var body = $('pmBody');
  body.innerHTML = '';

  for (var i = 0; i < pm.length; i++) {
    var m = pm[i];
    var tr = document.createElement('tr');
    var cells = [
      m.minute, m.trials, m.misses, m.lapses,
      n1(m.avgRt), n1(m.medianRt), n1(m.sdRt), n1(m.fastest10Rt), n1(m.slowest10Rt),
      n3(m.avgRs), n3(m.medianRs),
      vel[i] == null ? '—' : n2(vel[i]),
      acc[i] == null ? '—' : n2(acc[i])
    ];
    for (var c = 0; c < cells.length; c++) {
      var td = document.createElement('td');
      td.textContent = cells[c];
      tr.appendChild(td);
    }
    body.appendChild(tr);
  }

  $('pmSlope').textContent = r.scored.dynamicsRs.slope == null
    ? 'Not enough minutes to fit a trend.'
    : 'Reaction-speed trend across the test: ' + n3(r.scored.dynamicsRs.slope) +
      ' units per minute (negative = slowing down).';
}

/* ---------------- exports ---------------- */

var PARTICIPANT_COLS = ['run_id', 'date', 'time', 'participant_id', 'name', 'address',
                        'birth_date', 'educational_level', 'session_label', 'trial_number'];

function participantVals(r) {
  var p = r.participant;
  return [r.runId, p.date, p.time, p.participantId, p.name, p.address,
          p.birthDate, p.education, p.sessionLabel, p.trialNumber];
}

var RAW_HEADER = PARTICIPANT_COLS.concat([
  'epoch_index', 'minute', 'onset_ms', 'responded', 'rt_ms', 'rs_per_sec',
  'outcome', 'lapse', 'late_response', 'anticipation', 'extra_responses'
]);

function rawRows(r) {
  var pv = participantVals(r);
  var extras = r.extras || [];
  return r.scored.trials.map(function (t, i) {
    return pv.concat([
      t.index, t.minute + 1, round(t.onsetMs, 2),
      t.rtMs === null ? 0 : 1,
      round(t.rtMs, 3), round(t.rsPerSec, 5),
      t.outcome, t.lapse, t.lateResponse, t.anticipation,
      extras[i] == null ? '' : extras[i]
    ]);
  });
}

var PM_HEADER = PARTICIPANT_COLS.concat([
  'minute', 'trials', 'hits', 'misses', 'lapses', 'late_responses', 'hit_ratio',
  'avg_rt', 'median_rt', 'stdev_rt', 'fastest10_rt', 'slowest10_rt',
  'corr_avg_rt', 'corr_median_rt', 'corr_stdev_rt', 'corr_fastest10_rt', 'corr_slowest10_rt',
  'avg_rs', 'median_rs', 'stdev_rs', 'fastest10_rs', 'slowest10_rs',
  'corr_avg_rs', 'corr_median_rs', 'corr_stdev_rs', 'corr_fastest10_rs', 'corr_slowest10_rs',
  'n_rt', 'n_rt_corrected', 'anticipations_removed', 'outliers_removed',
  'velocity_rs', 'acceleration_rs', 'velocity_rt', 'acceleration_rt'
]);

function pmRows(r) {
  var pv = participantVals(r);
  var vRs = r.scored.dynamicsRs, vRt = r.scored.dynamicsRt;
  return r.scored.perMinute.map(function (m, i) {
    return pv.concat([
      m.minute, m.trials, m.hits, m.misses, m.lapses, m.lateResponses, round(m.hitRatio, 4),
      round(m.avgRt, 2), round(m.medianRt, 2), round(m.sdRt, 2), round(m.fastest10Rt, 2), round(m.slowest10Rt, 2),
      round(m.corrAvgRt, 2), round(m.corrMedianRt, 2), round(m.corrSdRt, 2), round(m.corrFastest10Rt, 2), round(m.corrSlowest10Rt, 2),
      round(m.avgRs, 5), round(m.medianRs, 5), round(m.sdRs, 5), round(m.fastest10Rs, 5), round(m.slowest10Rs, 5),
      round(m.corrAvgRs, 5), round(m.corrMedianRs, 5), round(m.corrSdRs, 5), round(m.corrFastest10Rs, 5), round(m.corrSlowest10Rs, 5),
      m.n, m.nCorrected, m.nAnticipationsRemoved, m.nOutliersRemoved,
      round(vRs.velocity[i], 5), round(vRs.acceleration[i], 5),
      round(vRt.velocity[i], 3), round(vRt.acceleration[i], 3)
    ]);
  });
}

var SUMMARY_HEADER = PARTICIPANT_COLS.concat([
  'sleep_onset_ms', 'sleep_onset_criterion_ms', 'slept_before_max', 'end_reason', 'elapsed_ms',
  'total_trialrun', 'hit_ratio', 'total_hits', 'total_miss', 'total_lapse', 'late_responses',
  'ep_1_2', 'ep_3_6', 'ep_7plus', 'longest_miss_run',
  'total_avg_rt', 'total_median_rt', 'total_stdev_rt', 'fastest10_rt', 'slowest10_rt',
  'corr_total_avg_rt', 'corr_total_median_rt', 'corr_total_stdev_rt', 'corr_fastest10_rt', 'corr_slowest10_rt',
  'total_avg_rs', 'total_median_rs', 'total_stdev_rs', 'fastest10_rs', 'slowest10_rs',
  'corr_total_avg_rs', 'corr_total_median_rs', 'corr_total_stdev_rs', 'corr_fastest10_rs', 'corr_slowest10_rs',
  'n_rt', 'n_rt_corrected', 'anticipations_removed', 'outliers_removed', 'rs_undefined',
  'rs_slope_per_min', 'rt_slope_per_min',
  'isi_ms', 'stim_ms', 'hit_window_ms', 'lapse_threshold_ms', 'miss_criterion', 'max_minutes',
  'correction', 'anticipation_threshold_ms', 'sd_multiplier',
  'extra_responses', 'page_blur_count'
]);

function summaryRow(r) {
  var t = r.scored.totals, e = r.scored.errorProfiles, c = r.config;
  return participantVals(r).concat([
    round(r.sleepOnsetMs, 1), round(r.sleepOnsetCriterionMs, 1), r.sleptBeforeMax ? 1 : 0,
    r.endReason, round(r.elapsedMs, 1),
    t.trials, round(t.hitRatio, 4), t.hits, t.misses, t.lapses, t.lateResponses,
    e.ep1_2, e.ep3_6, e.ep7plus, e.longestRun,
    round(t.avgRt, 2), round(t.medianRt, 2), round(t.sdRt, 2), round(t.fastest10Rt, 2), round(t.slowest10Rt, 2),
    round(t.corrAvgRt, 2), round(t.corrMedianRt, 2), round(t.corrSdRt, 2), round(t.corrFastest10Rt, 2), round(t.corrSlowest10Rt, 2),
    round(t.avgRs, 5), round(t.medianRs, 5), round(t.sdRs, 5), round(t.fastest10Rs, 5), round(t.slowest10Rs, 5),
    round(t.corrAvgRs, 5), round(t.corrMedianRs, 5), round(t.corrSdRs, 5), round(t.corrFastest10Rs, 5), round(t.corrSlowest10Rs, 5),
    t.n, t.nCorrected, t.nAnticipationsRemoved, t.nOutliersRemoved, t.nRsUndefined,
    round(r.scored.dynamicsRs.slope, 5), round(r.scored.dynamicsRt.slope, 3),
    c.isiMs, c.stimMs, c.hitWindowMs, c.lapseMs, c.missCriterion, c.maxMinutes,
    c.correction, c.anticipationMs, c.sdMultiplier,
    r.extraResponses, r.blurCount
  ]);
}

function allRows(builder, header) {
  var all = loadSessions();
  if (!all.length) { alert('No saved trials to export.'); return null; }
  var rows = [];
  for (var i = 0; i < all.length; i++) rows = rows.concat(builder(all[i]));
  return toCsv(header, rows);
}

/* ---------------- wiring ---------------- */

$('btnStart').addEventListener('click', beginTrial);

$('btnAgain').addEventListener('click', function () {
  $('trialNumber').value = intVal('trialNumber', 1) + 1;
  show('screen-setup');
});

$('btnExportRaw').addEventListener('click', function () {
  if (lastResult) download(fileStem(lastResult) + '_raw.csv', toCsv(RAW_HEADER, rawRows(lastResult)));
});
$('btnExportPm').addEventListener('click', function () {
  if (lastResult) download(fileStem(lastResult) + '_perminute.csv', toCsv(PM_HEADER, pmRows(lastResult)));
});
$('btnExportSummary').addEventListener('click', function () {
  if (lastResult) download(fileStem(lastResult) + '_summary.csv', toCsv(SUMMARY_HEADER, [summaryRow(lastResult)]));
});

$('btnExportAllRaw').addEventListener('click', function () {
  var csv = allRows(rawRows, RAW_HEADER);
  if (csv) download('bsrt_all_raw.csv', csv);
});
$('btnExportAllPm').addEventListener('click', function () {
  var csv = allRows(pmRows, PM_HEADER);
  if (csv) download('bsrt_all_perminute.csv', csv);
});
$('btnExportAllSummary').addEventListener('click', function () {
  var all = loadSessions();
  if (!all.length) { alert('No saved trials to export.'); return; }
  download('bsrt_all_summaries.csv', toCsv(SUMMARY_HEADER, all.map(summaryRow)));
});

$('btnAbort').addEventListener('click', function (e) {
  e.stopPropagation();
  if (running) endTask('aborted');
});

$('screen-task').addEventListener('pointerdown', function (e) {
  if (e.target && e.target.id === 'btnAbort') return;
  handleResponse();
});

document.addEventListener('keydown', function (e) {
  if (!running) return;
  if (e.code === 'Space' || e.key === ' ' || e.code === 'Enter' || e.key === 'Enter') {
    e.preventDefault();
    handleResponse();
  }
});

document.addEventListener('visibilitychange', function () {
  if (running && document.hidden) blurCount += 1;
  if (!document.hidden && running) requestWakeLock();
});

/* Suggest the short-protocol lapse threshold when the ceiling drops below 5 min. */
$('maxMinutes').addEventListener('change', function () {
  var short = intVal('maxMinutes', 40) < 5;
  $('lapseHint').textContent = short
    ? 'Short protocol (< 5 min): 355 ms is the conventional lapse threshold.'
    : 'Regular protocol: 500 ms is the conventional lapse threshold.';
});

refreshCount();
