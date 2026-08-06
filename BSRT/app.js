'use strict';

/* BSRT — Behavioral Sleep Resistance Task
 * A software implementation of the OSLER (Oxford Sleep Resistance Test) paradigm.
 * Vanilla JS, no build step, no network. Data stays in this browser.
 */

var LS_KEY = 'bsrt.sessions.v1';
/* A "lapse" is a response that arrived, but too slowly to count as alert
 * responding. The conventional PVT threshold is 500 ms and OSLER work often
 * uses 1000 ms, but short protocols need a tighter cut (355 ms has been used
 * for 3-minute tests), so this is a per-trial setting rather than a constant. */
var DEFAULT_LAPSE_MS = 1000;

/* ---------------- state ---------------- */

var cfg = null;          // frozen config for the running trial
var meta = null;         // participant / session identifiers
var epochs = [];         // one record per stimulus presentation
var startTime = 0;       // performance.now() at trial start
var startWall = null;    // wall-clock ISO string at trial start
var currentEpoch = null;
var epochTimer = null;
var stimTimer = null;
var consecutiveMisses = 0;
var longestMissRun = 0;
var missRunStartIndex = -1;
var running = false;
var blurCount = 0;
var wakeLock = null;
var lastResult = null;

/* ---------------- small helpers ---------------- */

function $(id) { return document.getElementById(id); }

function show(screenId) {
  var screens = document.querySelectorAll('.screen');
  for (var i = 0; i < screens.length; i++) screens[i].classList.remove('active');
  $(screenId).classList.add('active');
}

function intVal(id, fallback) {
  var n = parseInt($(id).value, 10);
  return isFinite(n) ? n : fallback;
}

function mean(a) {
  if (!a.length) return null;
  var s = 0;
  for (var i = 0; i < a.length; i++) s += a[i];
  return s / a.length;
}

function sd(a) {
  if (a.length < 2) return null;
  var m = mean(a), s = 0;
  for (var i = 0; i < a.length; i++) s += (a[i] - m) * (a[i] - m);
  return Math.sqrt(s / (a.length - 1));
}

function median(a) {
  if (!a.length) return null;
  var b = a.slice().sort(function (x, y) { return x - y; });
  var mid = Math.floor(b.length / 2);
  return b.length % 2 ? b[mid] : (b[mid - 1] + b[mid]) / 2;
}

function fmtClock(ms) {
  if (ms == null) return '—';
  var total = Math.round(ms / 1000);
  var m = Math.floor(total / 60), s = total % 60;
  return m + ':' + (s < 10 ? '0' : '') + s;
}

function fmtMs(v) { return v == null ? '—' : Math.round(v) + ' ms'; }

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

/* ---------------- storage ---------------- */

function loadSessions() {
  try {
    var raw = localStorage.getItem(LS_KEY);
    var arr = raw ? JSON.parse(raw) : [];
    return Array.isArray(arr) ? arr : [];
  } catch (e) {
    return [];
  }
}

function saveSession(result) {
  var all = loadSessions();
  all.push(result);
  try {
    localStorage.setItem(LS_KEY, JSON.stringify(all));
  } catch (e) {
    alert('Could not save this trial to browser storage (it may be full). Export the CSV now so the data is not lost.');
  }
  refreshSessionCount();
}

function refreshSessionCount() {
  $('sessionCount').textContent = loadSessions().length;
}

/* ---------------- screen wake / fullscreen ---------------- */

function requestWakeLock() {
  if (!navigator.wakeLock || !navigator.wakeLock.request) return;
  navigator.wakeLock.request('screen').then(function (lock) {
    wakeLock = lock;
  }).catch(function () { /* not critical */ });
}

function releaseWakeLock() {
  if (wakeLock) {
    try { wakeLock.release(); } catch (e) { /* ignore */ }
    wakeLock = null;
  }
}

function enterFullscreen() {
  var el = document.documentElement;
  if (el.requestFullscreen) {
    el.requestFullscreen().catch(function () { /* user may decline */ });
  }
}

function exitFullscreen() {
  if (document.fullscreenElement && document.exitFullscreen) {
    document.exitFullscreen().catch(function () { /* ignore */ });
  }
}

/* ---------------- the task ---------------- */

function beginTrial() {
  var maxMinutes = intVal('maxMinutes', 40);
  var isiMs = intVal('isiMs', 3000);
  var stimMs = intVal('stimMs', 1000);
  var missCriterion = intVal('missCriterion', 7);

  var err = null;
  if (isiMs < 200) err = 'Stimulus interval must be at least 200 ms.';
  else if (stimMs >= isiMs) err = 'Stimulus duration must be shorter than the stimulus interval.';
  else if (missCriterion < 1) err = 'The miss criterion must be at least 1.';
  else if (maxMinutes < 1) err = 'Maximum duration must be at least 1 minute.';

  if (err) {
    $('setupError').textContent = err;
    $('setupError').hidden = false;
    return;
  }
  $('setupError').hidden = true;

  cfg = {
    isiMs: isiMs,
    stimMs: stimMs,
    missCriterion: missCriterion,
    maxMinutes: maxMinutes,
    maxMs: maxMinutes * 60000,
    lapseMs: intVal('lapseMs', DEFAULT_LAPSE_MS),
    responseWindow: $('responseWindow').value
  };
  meta = {
    runId: 'bsrt-' + Date.now(),
    participantId: $('participantId').value.trim() || 'NA',
    sessionLabel: $('sessionLabel').value.trim() || 'NA',
    trialNumber: intVal('trialNumber', 1),
    userAgent: navigator.userAgent
  };

  epochs = [];
  currentEpoch = null;
  consecutiveMisses = 0;
  longestMissRun = 0;
  missRunStartIndex = -1;
  blurCount = 0;
  running = false;

  if ($('useFullscreen').checked) enterFullscreen();
  requestWakeLock();

  runCountdown();
}

function runCountdown() {
  show('screen-countdown');
  var n = 3;
  $('countdownNum').textContent = n;
  var iv = setInterval(function () {
    n -= 1;
    if (n <= 0) {
      clearInterval(iv);
      startTask();
    } else {
      $('countdownNum').textContent = n;
    }
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

/* Schedule against absolute targets so timing does not drift over 40 minutes. */
function scheduleEpoch(n) {
  var target = startTime + n * cfg.isiMs;
  var delay = target - performance.now();
  epochTimer = setTimeout(function () { runEpoch(n); }, delay > 0 ? delay : 0);
}

function runEpoch(n) {
  finalizeEpoch();
  if (!running) return;

  if (n * cfg.isiMs >= cfg.maxMs) {
    endTask('max_duration');
    return;
  }

  currentEpoch = {
    index: n,
    scheduledMs: n * cfg.isiMs,
    onsetMs: performance.now() - startTime,
    onsetPerf: performance.now(),
    responded: false,
    rtMs: null,
    extra: 0,
    late: 0
  };
  epochs.push(currentEpoch);

  ledOn();
  stimTimer = setTimeout(ledOff, cfg.stimMs);

  scheduleEpoch(n + 1);
}

function finalizeEpoch() {
  if (!currentEpoch) return;
  var ep = currentEpoch;
  currentEpoch = null;

  if (ep.responded) {
    consecutiveMisses = 0;
    missRunStartIndex = -1;
    return;
  }

  if (consecutiveMisses === 0) missRunStartIndex = ep.index;
  consecutiveMisses += 1;
  if (consecutiveMisses > longestMissRun) longestMissRun = consecutiveMisses;

  if (consecutiveMisses >= cfg.missCriterion) endTask('sleep_onset');
}

function handleResponse() {
  if (!running) return;
  if (!currentEpoch) return; // response outside any epoch (should not occur)

  var elapsed = performance.now() - currentEpoch.onsetPerf;

  // Under the 'stimulus' response window, a press after the light has gone out
  // is not a response to that stimulus, so the epoch stays a miss.
  if (cfg.responseWindow === 'stimulus' && elapsed > cfg.stimMs) {
    currentEpoch.late += 1;
    return;
  }

  if (currentEpoch.responded) {
    currentEpoch.extra += 1;
    return;
  }
  currentEpoch.responded = true;
  currentEpoch.rtMs = elapsed;

  // Responding extinguishes the stimulus immediately. The epoch clock is
  // untouched: the next stimulus still starts at a fixed offset from this
  // one's onset, so total task duration does not depend on responding.
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
  exitFullscreen();

  lastResult = buildResult(reason);
  saveSession(lastResult);
  renderResult(lastResult);
  show('screen-results');
}

/* ---------------- scoring ---------------- */

function buildResult(reason) {
  var rts = [];
  var hits = 0, misses = 0, extra = 0;
  for (var i = 0; i < epochs.length; i++) {
    var ep = epochs[i];
    extra += ep.extra;
    if (ep.responded) { hits++; rts.push(ep.rtMs); } else { misses++; }
  }

  var lapses = 0;
  for (var j = 0; j < rts.length; j++) if (rts[j] > cfg.lapseMs) lapses++;

  var late = 0;
  for (var k = 0; k < epochs.length; k++) late += epochs[k].late;

  // Sleep latency. Primary definition: time from trial start to the onset of the
  // FIRST stimulus in the terminating run of consecutive misses. Also reported:
  // the time at which the criterion was actually confirmed (end of the last
  // missed epoch), which is `missCriterion * isi` later.
  var latencyFirstMissMs = null;
  var latencyCriterionMs = null;
  if (reason === 'sleep_onset' && missRunStartIndex >= 0) {
    latencyFirstMissMs = missRunStartIndex * cfg.isiMs;
    latencyCriterionMs = latencyFirstMissMs + cfg.missCriterion * cfg.isiMs;
  }

  var elapsedMs = epochs.length ? (epochs[epochs.length - 1].onsetMs + cfg.isiMs) : 0;

  return {
    runId: meta.runId,
    startedAt: startWall,
    participantId: meta.participantId,
    sessionLabel: meta.sessionLabel,
    trialNumber: meta.trialNumber,
    userAgent: meta.userAgent,
    config: cfg,
    endReason: reason,
    sleptBeforeMax: reason === 'sleep_onset',
    latencyFirstMissMs: latencyFirstMissMs,
    latencyCriterionMs: latencyCriterionMs,
    elapsedMs: elapsedMs,
    nEpochs: epochs.length,
    hits: hits,
    misses: misses,
    hitRate: epochs.length ? hits / epochs.length : null,
    longestMissRun: longestMissRun,
    meanRtMs: mean(rts),
    medianRtMs: median(rts),
    sdRtMs: sd(rts),
    lapses: lapses,
    lapseThresholdMs: cfg.lapseMs,
    lateResponses: late,
    extraResponses: extra,
    responseWindowLabel: cfg.responseWindow,
    blurCount: blurCount,
    epochs: epochs.map(function (e) {
      return {
        index: e.index,
        scheduledMs: Math.round(e.scheduledMs),
        onsetMs: Math.round(e.onsetMs),
        responded: e.responded ? 1 : 0,
        rtMs: e.rtMs == null ? null : Math.round(e.rtMs),
        lapse: e.rtMs == null ? null : (e.rtMs > cfg.lapseMs ? 1 : 0),
        extra: e.extra,
        late: e.late
      };
    })
  };
}

/* ---------------- results rendering ---------------- */

function renderResult(r) {
  $('resultWho').textContent =
    r.participantId + ' · ' + r.sessionLabel + ' · trial ' + r.trialNumber;

  if (r.sleptBeforeMax) {
    $('mLatency').textContent = fmtClock(r.latencyFirstMissMs);
    $('mLatencyFoot').textContent =
      'to the first of ' + r.config.missCriterion +
      ' consecutive missed stimuli (criterion confirmed at ' + fmtClock(r.latencyCriterionMs) + ')';
    $('mOutcome').textContent =
      'Sleep onset was scored before the ' + r.config.maxMinutes + '-minute ceiling.';
  } else {
    $('mLatency').textContent = '> ' + fmtClock(r.config.maxMs);
    $('mLatencyFoot').textContent = 'no sleep onset scored — censored at the maximum duration';
    $('mOutcome').textContent = r.endReason === 'aborted'
      ? 'The trial was ended early by the experimenter, so this is not a valid latency.'
      : 'The participant stayed awake for the full ' + r.config.maxMinutes + ' minutes.';
  }

  $('mEpochs').textContent = r.nEpochs;
  $('mHits').textContent = r.hits;
  $('mMisses').textContent = r.misses;
  $('mHitRate').textContent = r.hitRate == null ? '—' : (r.hitRate * 100).toFixed(1) + '%';
  $('mLongestRun').textContent = r.longestMissRun;
  $('mMeanRt').textContent = fmtMs(r.meanRtMs);
  $('mMedianRt').textContent = fmtMs(r.medianRtMs);
  $('mSdRt').textContent = fmtMs(r.sdRtMs);
  $('lapseLabel').textContent = 'Lapses (>' + r.lapseThresholdMs + ' ms)';
  $('mLapses').textContent = r.lapses;
  $('mLate').textContent = r.responseWindowLabel === 'stimulus'
    ? r.lateResponses + ' (after the light went out — not counted)'
    : 'n/a (full-epoch window)';
  $('mExtra').textContent = r.extraResponses;
  $('mBlur').textContent = r.blurCount ? r.blurCount + '×' : 'no';
}

/* ---------------- exports ---------------- */

var EPOCH_HEADER = [
  'run_id', 'participant_id', 'session_label', 'trial_number', 'started_at',
  'epoch_index', 'scheduled_onset_ms', 'actual_onset_ms', 'responded', 'rt_ms', 'lapse',
  'extra_responses', 'late_responses'
];

function epochRows(r) {
  return r.epochs.map(function (e) {
    return [
      r.runId, r.participantId, r.sessionLabel, r.trialNumber, r.startedAt,
      e.index, e.scheduledMs, e.onsetMs, e.responded, e.rtMs, e.lapse, e.extra, e.late
    ];
  });
}

var SUMMARY_HEADER = [
  'run_id', 'participant_id', 'session_label', 'trial_number', 'started_at',
  'isi_ms', 'stim_ms', 'miss_criterion', 'max_minutes',
  'lapse_threshold_ms', 'response_window',
  'end_reason', 'slept_before_max',
  'latency_first_miss_ms', 'latency_criterion_ms', 'elapsed_ms',
  'n_epochs', 'hits', 'misses', 'hit_rate', 'longest_miss_run',
  'mean_rt_ms', 'median_rt_ms', 'sd_rt_ms', 'lapses', 'extra_responses', 'late_responses',
  'page_blur_count'
];

function round(v, dp) { return v == null ? null : Number(v.toFixed(dp)); }

function summaryRow(r) {
  return [
    r.runId, r.participantId, r.sessionLabel, r.trialNumber, r.startedAt,
    r.config.isiMs, r.config.stimMs, r.config.missCriterion, r.config.maxMinutes,
    r.config.lapseMs, r.config.responseWindow,
    r.endReason, r.sleptBeforeMax ? 1 : 0,
    r.latencyFirstMissMs, r.latencyCriterionMs, Math.round(r.elapsedMs),
    r.nEpochs, r.hits, r.misses, round(r.hitRate, 4), r.longestMissRun,
    round(r.meanRtMs, 1), round(r.medianRtMs, 1), round(r.sdRtMs, 1),
    r.lapses, r.extraResponses, r.lateResponses, r.blurCount
  ];
}

function exportEpochs(r) {
  download(
    'bsrt_' + safeName(r.participantId) + '_' + safeName(r.sessionLabel) + '_t' + r.trialNumber + '_epochs.csv',
    toCsv(EPOCH_HEADER, epochRows(r))
  );
}

function exportSummary(r) {
  download(
    'bsrt_' + safeName(r.participantId) + '_' + safeName(r.sessionLabel) + '_t' + r.trialNumber + '_summary.csv',
    toCsv(SUMMARY_HEADER, [summaryRow(r)])
  );
}

function exportAllEpochs() {
  var all = loadSessions();
  if (!all.length) { alert('No saved trials to export.'); return; }
  var rows = [];
  for (var i = 0; i < all.length; i++) rows = rows.concat(epochRows(all[i]));
  download('bsrt_all_epochs.csv', toCsv(EPOCH_HEADER, rows));
}

function exportAllSummaries() {
  var all = loadSessions();
  if (!all.length) { alert('No saved trials to export.'); return; }
  download('bsrt_all_summaries.csv', toCsv(SUMMARY_HEADER, all.map(summaryRow)));
}

/* ---------------- sessions screen ---------------- */

function renderSessions() {
  var all = loadSessions();
  var box = $('sessionList');
  box.innerHTML = '';

  if (!all.length) {
    var p = document.createElement('p');
    p.className = 'empty';
    p.textContent = 'No trials saved yet.';
    box.appendChild(p);
    return;
  }

  for (var i = all.length - 1; i >= 0; i--) {
    var r = all[i];
    var row = document.createElement('div');
    row.className = 'session-row';

    var who = document.createElement('span');
    who.className = 'who';
    who.textContent = r.participantId;

    var m = document.createElement('span');
    m.className = 'meta';
    m.textContent = r.sessionLabel + ' · trial ' + r.trialNumber + ' · ' +
      new Date(r.startedAt).toLocaleString();

    var lat = document.createElement('span');
    lat.className = 'lat';
    lat.textContent = r.sleptBeforeMax
      ? fmtClock(r.latencyFirstMissMs)
      : '> ' + fmtClock(r.config.maxMs);

    row.appendChild(who);
    row.appendChild(m);
    row.appendChild(lat);
    box.appendChild(row);
  }
}

/* ---------------- wiring ---------------- */

$('btnStart').addEventListener('click', beginTrial);

$('btnSessions').addEventListener('click', function () {
  renderSessions();
  show('screen-sessions');
});

$('btnBack').addEventListener('click', function () { show('screen-setup'); });

$('btnAgain').addEventListener('click', function () {
  $('trialNumber').value = intVal('trialNumber', 1) + 1;
  show('screen-setup');
});

$('btnExportEpochs').addEventListener('click', function () {
  if (lastResult) exportEpochs(lastResult);
});
$('btnExportSummary').addEventListener('click', function () {
  if (lastResult) exportSummary(lastResult);
});
$('btnExportAll').addEventListener('click', exportAllEpochs);
$('btnExportAllSummary').addEventListener('click', exportAllSummaries);

$('btnClearAll').addEventListener('click', function () {
  if (!loadSessions().length) return;
  if (confirm('Delete all saved trials from this browser? Export first — this cannot be undone.')) {
    localStorage.removeItem(LS_KEY);
    renderSessions();
    refreshSessionCount();
  }
});

$('btnAbort').addEventListener('click', function (e) {
  e.stopPropagation();
  if (running) endTask('aborted');
});

/* Responses: pointer anywhere on the task screen, or the spacebar. */
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

/* Flag if the participant navigates away mid-trial — it invalidates timing. */
document.addEventListener('visibilitychange', function () {
  if (running && document.hidden) blurCount += 1;
  if (!document.hidden && running) requestWakeLock();
});

refreshSessionCount();
