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
var L = window.BSRTi18n;

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
var pendingResult = null;     // held while the post-task KSS is answered
var presses = [];             // every keypress, for the integrity check
var kssBefore = null;
var kssAfter = null;
var kssStage = null;          // 'before' | 'after'
var clockRaf = null;
var audioCtx = null;

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
  function applyLanguage() {
  L.setLanguage($('language').value);
  L.applyTranslations(document);
  applyMode();
  if (lastResult) renderResult(lastResult);
}

function applyMode() {
  var mode = $('mode').value;
  document.body.setAttribute('data-mode', mode);
  // A PVT runs for a fixed duration; the sleep-onset criterion is an OSLER idea.
  $('criterionOn').checked = mode !== 'pvt';
  $('modeNote').textContent = mode === 'pvt'
    ? 'PVT: intervals vary within each block, drawn so that every block contains one of each. The stimulus is a millisecond counter; the sleep-onset criterion is off by default.'
    : 'BSRT / OSLER: a fixed interval between stimuli, with sleep onset scored from consecutive misses.';
  $('instructionsText').textContent = L.t(mode === 'pvt' ? 'instructions.pvt' : 'instructions.bsrt');
  $('taskHint').textContent = L.t(mode === 'pvt' ? 'task.hintPvt' : 'task.hintBsrt');
}

$('mode').addEventListener('change', applyMode);
$('language').addEventListener('change', applyLanguage);
applyLanguage();
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

/* ---------------- alarm ---------------- */

/*
 * Synthesised rather than a sound file, so the app stays dependency-free and
 * works offline. The context is created on the Start click, which is the user
 * gesture browsers require before audio may play.
 */
function initAudio() {
  if (audioCtx) return;
  try {
    var Ctx = window.AudioContext || window.webkitAudioContext;
    if (Ctx) audioCtx = new Ctx();
  } catch (e) { audioCtx = null; }
}

function playAlarm() {
  if (!cfg || !cfg.alarm || !audioCtx) return;
  try {
    if (audioCtx.state === 'suspended') audioCtx.resume();
    var t0 = audioCtx.currentTime;
    for (var i = 0; i < 3; i++) {
      var start = t0 + i * 0.35;
      var osc = audioCtx.createOscillator();
      var gain = audioCtx.createGain();
      osc.type = 'square';
      osc.frequency.setValueAtTime(880, start);
      gain.gain.setValueAtTime(0.0001, start);
      gain.gain.exponentialRampToValueAtTime(0.25, start + 0.02);
      gain.gain.exponentialRampToValueAtTime(0.0001, start + 0.30);
      osc.connect(gain);
      gain.connect(audioCtx.destination);
      osc.start(start);
      osc.stop(start + 0.32);
    }
  } catch (e) { /* audio is a convenience, never a reason to lose a trial */ }
}

/* ---------------- trial ---------------- */

function parseIsiSet() {
  return $('isiSetMs').value.split(/[,\s]+/)
    .map(function (v) { return Math.round(parseFloat(v) * 1000); })
    .filter(function (v) { return isFinite(v) && v > 0; });
}

function validate() {
  var mode = $('mode').value;
  var stim = intVal('stimMs', 1000);
  var shortestIsi;

  if (mode === 'pvt') {
    var set = parseIsiSet();
    if (set.length < 2) return 'Give at least two inter-stimulus intervals, in seconds (e.g. 2, 4, 6, 8, 10).';
    if (intVal('blockMs', 30) < 1) return 'Block length must be at least 1 second.';
    shortestIsi = Math.min.apply(null, set);
  } else {
    shortestIsi = intVal('isiMs', 3000);
    if (shortestIsi < 200) return 'Stimulus interval must be at least 200 ms.';
  }

  if (stim >= shortestIsi) return 'Stimulus duration must be shorter than the shortest stimulus interval (' + shortestIsi + ' ms).';
  if (intVal('maxMinutes', 40) < 1) return 'Maximum duration must be at least 1 minute.';
  if (intVal('lapseMs', 500) >= stim) return 'The lapse threshold must be below the hit window (stimulus duration).';
  return null;
}

function beginTrial() {
  var err = validate();
  if (err) { $('setupError').textContent = err; $('setupError').hidden = false; return; }
  $('setupError').hidden = true;

  var corr = $('correction').value; // none | falseStarts | outliers | both
  var mode = $('mode').value;
  var now = new Date();
  var seedTyped = parseInt($('seed').value, 10);
  var seed = isFinite(seedTyped) ? seedTyped : Math.floor(Math.random() * 2147483647);

  cfg = {
    mode: mode,
    isiSetMs: mode === 'pvt' ? parseIsiSet() : null,
    blockMs: intVal('blockMs', 30) * 1000,
    seed: seed,
    isiMs: intVal('isiMs', 3000),
    stimMs: intVal('stimMs', 1000),
    hitWindowMs: intVal('stimMs', 1000),   // RT <= stimulus duration counts as a hit
    lapseMs: intVal('lapseMs', 500),
    /* 0 disables early termination; a PVT runs for its full duration. */
    missCriterion: $('criterionOn').checked ? intVal('missCriterion', 7) : 0,
    maxMinutes: intVal('maxMinutes', 40),
    maxMs: intVal('maxMinutes', 40) * 60000,
    correction: corr,
    removeFalseStarts: corr === 'falseStarts' || corr === 'both',
    removeOutliers: corr === 'outliers' || corr === 'both',
    falseStartMs: numVal('falseStartMs', 100),
    sdMultiplier: numVal('sdMultiplier', 2),
    alarm: $('alarmOn').checked,
    kssWhen: $('kssWhen').value,
    feedbackMs: 1000,
    language: L.getLanguage()
  };

  cfg.schedule = S.buildSchedule({
    mode: cfg.mode,
    isiMs: cfg.isiMs,
    isiSetMs: cfg.isiSetMs,
    blockMs: cfg.blockMs,
    maxMs: cfg.maxMs,
    seed: cfg.seed
  });

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
  presses = [];
  currentEpoch = null;
  consecutiveMisses = 0;
  missRunStartIndex = -1;
  blurCount = 0;
  running = false;
  kssBefore = null;
  kssAfter = null;

  initAudio();
  requestWakeLock();
  if ($('useFullscreen').checked && document.documentElement.requestFullscreen) {
    document.documentElement.requestFullscreen().catch(function () {});
  }
  // The KSS, when requested, is answered before anything else happens.
  if (cfg.kssWhen === 'before' || cfg.kssWhen === 'both') showKss('before');
  else countdown();
}

/* ---------------- Karolinska Sleepiness Scale ---------------- */

function showKss(stage) {
  kssStage = stage;
  $('kssTitle').textContent = L.t(stage === 'before' ? 'kss.beforeTitle' : 'kss.afterTitle');
  $('kssQuestion').textContent = L.t('kss.question');
  $('kssInstruction').textContent = L.t('kss.instruction');

  var anchors = L.kssAnchors();
  var box = $('kssOptions');
  box.innerHTML = '';
  anchors.forEach(function (label, i) {
    var b = document.createElement('button');
    b.className = 'kss-option';
    b.type = 'button';
    var n = document.createElement('span');
    n.className = 'kss-num';
    n.textContent = i + 1;
    var tx = document.createElement('span');
    tx.textContent = label;
    b.appendChild(n);
    b.appendChild(tx);
    b.addEventListener('click', function () { answerKss(i + 1); });
    box.appendChild(b);
  });
  show('screen-kss');
}

function answerKss(value) {
  if (kssStage === 'before') {
    kssBefore = value;
    countdown();
  } else {
    kssAfter = value;
    finishResult();
  }
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
  startClock();
  runEpoch(0);
}

/*
 * Absolute targets taken from the pre-built schedule, so presentation does not
 * drift across a long trial and a variable-interval (PVT) schedule is delivered
 * exactly as generated.
 */
function scheduleEpoch(n) {
  var last = cfg.schedule.nStimuli - 1;
  // Past the final stimulus there is no onset to aim at, but the trial still
  // has to be ended — so schedule one last call at the end of the last epoch.
  // Returning here instead would leave the task running with nothing to fire.
  var target = n < cfg.schedule.nStimuli
    ? cfg.schedule.onsets[n]
    : cfg.schedule.onsets[last] + cfg.schedule.isis[last];
  var delay = (startTime + target) - performance.now();
  epochTimer = setTimeout(function () { runEpoch(n); }, delay > 0 ? delay : 0);
}

function runEpoch(n) {
  finalizeEpoch();
  if (!running) return;
  if (n >= cfg.schedule.nStimuli) { endTask('max_duration'); return; }

  currentEpoch = {
    index: n,
    onsetMs: performance.now() - startTime,
    onsetPerf: performance.now(),
    block: cfg.schedule.blocks[n],
    epochIsiMs: cfg.schedule.isis[n],
    isiBeforeMs: n === 0 ? null : cfg.schedule.isis[n - 1],
    responded: false,
    rtMs: null,
    extra: 0,
    frozen: false
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
  if (cfg.missCriterion > 0 && consecutiveMisses >= cfg.missCriterion) endTask('sleep_onset');
}

/*
 * Every response is recorded with its raw RT, however late — the raw data must
 * survive for downstream analysis. Whether it counts as a hit is decided at
 * scoring time, not here.
 */
function handleResponse() {
  if (!running) return;
  var now = performance.now();

  // Every press is logged, including presses outside any epoch, because the
  // integrity check needs the full pattern — not just the scored responses.
  presses.push({ tMs: now - startTime, epochIndex: currentEpoch ? currentEpoch.index : null });

  if (!currentEpoch) return;
  if (currentEpoch.responded) { currentEpoch.extra += 1; return; }

  currentEpoch.responded = true;
  currentEpoch.rtMs = now - currentEpoch.onsetPerf;

  clearTimeout(stimTimer);
  if (cfg.mode === 'pvt') {
    // PVT convention: the counter stops on the response and the achieved
    // reaction time is shown back as feedback, then clears.
    currentEpoch.frozen = true;
    $('clock').textContent = Math.max(0, Math.round(currentEpoch.rtMs));
    $('clock').classList.add('frozen');
    stimTimer = setTimeout(ledOff, Math.min(cfg.feedbackMs, cfg.hitWindowMs));
  } else {
    // BSRT: responding extinguishes the light immediately. Either way the
    // epoch clock is untouched.
    ledOff();
  }
}

function ledOn() {
  if (cfg.mode === 'pvt') {
    var el = $('clock');
    el.classList.remove('frozen');
    el.textContent = '0';
    el.hidden = false;
  } else {
    $('led').classList.add('on');
  }
}

function ledOff() {
  $('led').classList.remove('on');
  var el = $('clock');
  el.hidden = true;
  el.classList.remove('frozen');
}

/*
 * The PVT counter is display only — it never feeds the reaction time, which is
 * still measured from stimulus onset to the keypress. It runs on its own
 * animation frame so redrawing digits cannot disturb epoch scheduling.
 */
function startClock() {
  if (cfg.mode !== 'pvt') return;
  function tick() {
    if (!running) return;
    clockRaf = requestAnimationFrame(tick);
    if (!currentEpoch || currentEpoch.frozen) return;
    var el = $('clock');
    if (el.hidden) return;
    var elapsed = performance.now() - currentEpoch.onsetPerf;
    if (elapsed < 0) elapsed = 0;
    if (elapsed > cfg.hitWindowMs) elapsed = cfg.hitWindowMs;
    el.textContent = Math.round(elapsed);
  }
  clockRaf = requestAnimationFrame(tick);
}

function endTask(reason) {
  running = false;
  clearTimeout(epochTimer);
  clearTimeout(stimTimer);
  if (clockRaf) cancelAnimationFrame(clockRaf);
  clockRaf = null;
  ledOff();
  releaseWakeLock();
  if (document.fullscreenElement && document.exitFullscreen) document.exitFullscreen().catch(function () {});

  if (reason === 'sleep_onset') playAlarm();

  pendingResult = buildResult(reason);
  if (cfg.kssWhen === 'after' || cfg.kssWhen === 'both') showKss('after');
  else finishResult();
}

/* Assembled only once the post-task KSS (if any) has been answered. */
function finishResult() {
  pendingResult.kssBefore = kssBefore;
  pendingResult.kssAfter = kssAfter;
  lastResult = pendingResult;
  pendingResult = null;
  saveSession(lastResult);
  renderResult(lastResult);
  show('screen-results');
}

/* ---------------- result assembly ---------------- */

function buildResult(reason) {
  var scored = S.score(epochs.map(function (e) {
    return {
      index: e.index, onsetMs: e.onsetMs, rtMs: e.rtMs,
      block: e.block, epochIsiMs: e.epochIsiMs, isiBeforeMs: e.isiBeforeMs
    };
  }), cfg, presses);

  var sleepOnsetMs = null, sleepOnsetCriterionMs = null;
  if (reason === 'sleep_onset' && missRunStartIndex >= 0) {
    var first = epochs[missRunStartIndex];
    var last = epochs[missRunStartIndex + cfg.missCriterion - 1];
    sleepOnsetMs = first ? first.onsetMs : cfg.schedule.onsets[missRunStartIndex];
    // With a variable schedule the criterion is confirmed one INTERVAL after the
    // last missed stimulus, and that interval differs per epoch.
    sleepOnsetCriterionMs = last ? last.onsetMs + last.epochIsiMs
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
    elapsedMs: epochs.length ? epochs[epochs.length - 1].onsetMs + epochs[epochs.length - 1].epochIsiMs : 0,
    extraResponses: epochs.reduce(function (n, e) { return n + e.extra; }, 0),
    kssBefore: null,
    kssAfter: null,
    language: cfg.language,
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
  $('oFalseStarts').textContent = t.falseStarts + ' (< ' + r.config.falseStartMs + ' ms)';
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
        t.nFalseStartsRemoved + ' false starts and ' + t.nOutliersRemoved +
        ' outliers from ' + t.n + ' hits.';

  renderIntegrity(r);
  renderKss(r);
  renderPerMinute(r);
}

function renderIntegrity(r) {
  var g = r.scored.integrity;
  $('inPresses').textContent = g.totalPresses;
  $('inExtra').textContent = g.extraPresses + (g.extraRate ? ' (' + (g.extraRate * 100).toFixed(0) + '% of trials)' : '');
  $('inBurst').textContent = g.burstMax + ' in ' + g.thresholds.burstWindowMs + ' ms';
  $('inRapid').textContent = g.rapidPairs + ' under ' + g.thresholds.rapidGapMs + ' ms apart';

  var box = $('integrityVerdict');
  if (g.suspected) {
    box.className = 'note warnbox';
    box.textContent = L.t('results.integrityFlag') + ' ' + g.reasons.join('; ') + '.';
  } else {
    box.className = 'note tight';
    box.textContent = L.t('results.integrityOk');
  }
}

function renderKss(r) {
  var card = $('kssResult');
  if (r.kssBefore == null && r.kssAfter == null) { card.hidden = true; return; }
  card.hidden = false;
  var anchors = L.kssAnchors(r.language);
  $('kssBeforeVal').textContent = r.kssBefore == null
    ? '—' : r.kssBefore + ' — ' + anchors[r.kssBefore - 1];
  $('kssAfterVal').textContent = r.kssAfter == null
    ? '—' : r.kssAfter + ' — ' + anchors[r.kssAfter - 1];
}

function describeCorrection(c) {
  if (c.correction === 'falseStarts') return 'RT < ' + c.falseStartMs + ' ms removed';
  if (c.correction === 'outliers') return 'RT beyond ' + c.sdMultiplier + ' SD removed';
  if (c.correction === 'both') return 'RT < ' + c.falseStartMs + ' ms removed, then RT beyond ' + c.sdMultiplier + ' SD';
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

var PARTICIPANT_COLS = ['run_id', 'date', 'time', 'language', 'participant_id', 'name', 'address',
                        'birth_date', 'educational_level', 'session_label', 'trial_number'];

function participantVals(r) {
  var p = r.participant;
  return [r.runId, p.date, p.time, r.language || 'en', p.participantId, p.name, p.address,
          p.birthDate, p.education, p.sessionLabel, p.trialNumber];
}

var RAW_HEADER = PARTICIPANT_COLS.concat([
  'epoch_index', 'block', 'minute', 'onset_ms', 'epoch_isi_ms', 'isi_before_ms',
  'responded', 'rt_ms', 'rs_per_sec',
  'outcome', 'lapse', 'late_response', 'false_start', 'extra_responses'
]);

function rawRows(r) {
  var pv = participantVals(r);
  var extras = r.extras || [];
  return r.scored.trials.map(function (t, i) {
    return pv.concat([
      t.index, t.block === null ? '' : t.block + 1, t.minute + 1, round(t.onsetMs, 2),
      t.epochIsiMs, t.isiBeforeMs,
      t.rtMs === null ? 0 : 1,
      round(t.rtMs, 3), round(t.rsPerSec, 5),
      t.outcome, t.lapse, t.lateResponse, t.falseStart,
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
  'n_rt', 'n_rt_corrected', 'false_starts_removed', 'outliers_removed',
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
      m.n, m.nCorrected, m.nFalseStartsRemoved, m.nOutliersRemoved,
      round(vRs.velocity[i], 5), round(vRs.acceleration[i], 5),
      round(vRt.velocity[i], 3), round(vRt.acceleration[i], 3)
    ]);
  });
}

var SUMMARY_HEADER = PARTICIPANT_COLS.concat([
  'sleep_onset_ms', 'sleep_onset_criterion_ms', 'slept_before_max', 'end_reason', 'elapsed_ms',
  'total_trialrun', 'hit_ratio', 'total_hits', 'total_miss', 'total_lapse', 'late_responses',
  'total_false_starts',
  'ep_1_2', 'ep_3_6', 'ep_7plus', 'longest_miss_run',
  'total_avg_rt', 'total_median_rt', 'total_stdev_rt', 'fastest10_rt', 'slowest10_rt',
  'corr_total_avg_rt', 'corr_total_median_rt', 'corr_total_stdev_rt', 'corr_fastest10_rt', 'corr_slowest10_rt',
  'total_avg_rs', 'total_median_rs', 'total_stdev_rs', 'fastest10_rs', 'slowest10_rs',
  'corr_total_avg_rs', 'corr_total_median_rs', 'corr_total_stdev_rs', 'corr_fastest10_rs', 'corr_slowest10_rs',
  'n_rt', 'n_rt_corrected', 'false_starts_removed', 'outliers_removed', 'rs_undefined',
  'rs_slope_per_min', 'rt_slope_per_min',
  'mode', 'isi_ms', 'isi_set_s', 'block_s', 'schedule_method', 'schedule_seed', 'scheduled_stimuli',
  'stim_ms', 'hit_window_ms', 'lapse_threshold_ms', 'miss_criterion', 'max_minutes',
  'correction', 'false_start_threshold_ms', 'sd_multiplier',
  'kss_when', 'kss_before', 'kss_after',
  'total_presses', 'extra_presses', 'burst_max', 'rapid_pairs', 'cheating_suspected', 'cheating_reasons',
  'alarm_enabled', 'extra_responses', 'page_blur_count'
]);

function summaryRow(r) {
  var t = r.scored.totals, e = r.scored.errorProfiles, c = r.config;
  var g = r.scored.integrity;
  return participantVals(r).concat([
    round(r.sleepOnsetMs, 1), round(r.sleepOnsetCriterionMs, 1), r.sleptBeforeMax ? 1 : 0,
    r.endReason, round(r.elapsedMs, 1),
    t.trials, round(t.hitRatio, 4), t.hits, t.misses, t.lapses, t.lateResponses,
    t.falseStarts,
    e.ep1_2, e.ep3_6, e.ep7plus, e.longestRun,
    round(t.avgRt, 2), round(t.medianRt, 2), round(t.sdRt, 2), round(t.fastest10Rt, 2), round(t.slowest10Rt, 2),
    round(t.corrAvgRt, 2), round(t.corrMedianRt, 2), round(t.corrSdRt, 2), round(t.corrFastest10Rt, 2), round(t.corrSlowest10Rt, 2),
    round(t.avgRs, 5), round(t.medianRs, 5), round(t.sdRs, 5), round(t.fastest10Rs, 5), round(t.slowest10Rs, 5),
    round(t.corrAvgRs, 5), round(t.corrMedianRs, 5), round(t.corrSdRs, 5), round(t.corrFastest10Rs, 5), round(t.corrSlowest10Rs, 5),
    t.n, t.nCorrected, t.nFalseStartsRemoved, t.nOutliersRemoved, t.nRsUndefined,
    round(r.scored.dynamicsRs.slope, 5), round(r.scored.dynamicsRt.slope, 3),
    c.mode, c.mode === 'pvt' ? '' : c.isiMs,
    c.isiSetMs ? c.isiSetMs.map(function (v) { return v / 1000; }).join(' ') : '',
    c.blockMs / 1000, c.schedule.method, c.schedule.seed, c.schedule.nStimuli,
    c.stimMs, c.hitWindowMs, c.lapseMs, c.missCriterion, c.maxMinutes,
    c.correction, c.falseStartMs, c.sdMultiplier,
    c.kssWhen, r.kssBefore, r.kssAfter,
    g.totalPresses, g.extraPresses, g.burstMax, g.rapidPairs,
    g.suspected ? 1 : 0, g.reasons.join('; '),
    c.alarm ? 1 : 0, r.extraResponses, r.blurCount
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

function applyLanguage() {
  L.setLanguage($('language').value);
  L.applyTranslations(document);
  applyMode();
  if (lastResult) renderResult(lastResult);
}

function applyMode() {
  var mode = $('mode').value;
  document.body.setAttribute('data-mode', mode);
  // A PVT runs for a fixed duration; the sleep-onset criterion is an OSLER idea.
  $('criterionOn').checked = mode !== 'pvt';
  $('modeNote').textContent = mode === 'pvt'
    ? 'PVT: intervals vary within each block, drawn so that every block contains one of each. The stimulus is a millisecond counter; the sleep-onset criterion is off by default.'
    : 'BSRT / OSLER: a fixed interval between stimuli, with sleep onset scored from consecutive misses.';
  $('instructionsText').textContent = L.t(mode === 'pvt' ? 'instructions.pvt' : 'instructions.bsrt');
  $('taskHint').textContent = L.t(mode === 'pvt' ? 'task.hintPvt' : 'task.hintBsrt');
}

$('mode').addEventListener('change', applyMode);
$('language').addEventListener('change', applyLanguage);
applyLanguage();
refreshCount();
