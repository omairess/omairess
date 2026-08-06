'use strict';

/* BSRT desktop — task logic.
 *
 * Trial structure and scoring are identical to the browser build (scoring.js is
 * shared byte-for-byte). What differs is timing: stimuli are scheduled by
 * counting real vsync intervals, and onsets come from the frame clock.
 * See timing.js and README.md.
 */

const T = window.BSRTTiming;
const S = window.BSRTScoring;
const LS_KEY = 'bsrt.desktop.sessions.v1';
const PROBE_TAPS = 5;

/* ---------------- state ---------------- */

let displayInfo = null;
let calibration = null;
let calGrade = null;
let inputProbeResult = null;
let useEventStamp = true;

let cfg = null;
let meta = null;
let epochs = [];
let currentEpoch = null;

let running = false;
let rafId = null;
let frameIndex = -1;
let lastFrameTs = null;
let trialStartTs = null;
let startWall = null;
let droppedFrames = 0;

let framesPerEpoch = 0;
let framesStimOn = 0;
let achievedIsiMs = 0;
let maxEpochs = 0;

let consecutiveMisses = 0;
let missRunStartIndex = -1;
let blurCount = 0;
let lastResult = null;

/* ---------------- helpers ---------------- */

const $ = (id) => document.getElementById(id);

function show(id) {
  document.querySelectorAll('.screen').forEach((s) => s.classList.remove('active'));
  $(id).classList.add('active');
}

const intVal = (id, fb) => { const n = parseInt($(id).value, 10); return isFinite(n) ? n : fb; };
const numVal = (id, fb) => { const n = parseFloat($(id).value); return isFinite(n) ? n : fb; };
const txtVal = (id) => $(id).value.trim();

function fmtClock(ms) {
  if (ms == null) return '—';
  const t = Math.round(ms / 1000), m = Math.floor(t / 60), s = t % 60;
  return m + ':' + (s < 10 ? '0' : '') + s;
}

const n1 = (v) => (v == null ? '—' : v.toFixed(1));
const n2 = (v) => (v == null ? '—' : v.toFixed(2));
const n3 = (v) => (v == null ? '—' : v.toFixed(3));
const round = (v, dp) => (v == null ? null : Number(v.toFixed(dp)));

/* ---------------- CSV ---------------- */

function csvCell(v) {
  if (v === null || v === undefined) return '';
  const s = String(v);
  return /[",\n\r]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
}

const toCsv = (header, rows) =>
  [header.map(csvCell).join(',')].concat(rows.map((r) => r.map(csvCell).join(','))).join('\r\n') + '\r\n';

async function saveCsv(name, text) {
  const res = await window.bsrt.saveCsv(name, text);
  if (!res.ok && res.reason !== 'canceled') alert('Could not save: ' + res.reason);
}

const safeName = (s) => String(s || 'NA').replace(/[^A-Za-z0-9_-]/g, '_');
const fileStem = (r) => 'bsrt_' + safeName(r.participant.participantId) + '_' +
                        safeName(r.participant.sessionLabel) + '_t' + r.participant.trialNumber;

/* ---------------- storage ---------------- */

function loadSessions() {
  try {
    const a = JSON.parse(localStorage.getItem(LS_KEY) || '[]');
    return Array.isArray(a) ? a : [];
  } catch (e) { return []; }
}

function saveSession(r) {
  const all = loadSessions();
  all.push(r);
  try { localStorage.setItem(LS_KEY, JSON.stringify(all)); }
  catch (e) { alert('Could not save to local storage. Export this trial now.'); }
  $('sessionCount').textContent = loadSessions().length;
}

/* ---------------- startup ---------------- */

async function init() {
  $('sessionCount').textContent = loadSessions().length;
  try { displayInfo = await window.bsrt.getDisplayInfo(); }
  catch (e) { displayInfo = null; }
  renderEnvironment();
}

function renderEnvironment() {
  if (!displayInfo) { $('envInfo').textContent = 'Display information unavailable.'; return; }
  const d = displayInfo.current;
  $('envInfo').textContent = [
    d.displayFrequency ? d.displayFrequency + ' Hz (reported)' : 'refresh rate unknown',
    d.bounds.width + '×' + d.bounds.height,
    'scale ' + d.scaleFactor + '×',
    displayInfo.platform === 'darwin' ? 'macOS' : displayInfo.platform === 'win32' ? 'Windows' : displayInfo.platform,
    'Electron ' + displayInfo.electron
  ].join(' · ');
  if (displayInfo.displayCount > 1) {
    $('envWarn').hidden = false;
    $('envWarn').textContent = displayInfo.displayCount +
      ' displays detected. Run the task on one screen only — mirroring or extending can force the compositor to the slowest display.';
  }
}

/* ---------------- calibration ---------------- */

function validateSetup() {
  const isi = intVal('isiMs', 3000), stim = intVal('stimMs', 1000);
  if (isi < 200) return 'Stimulus interval must be at least 200 ms.';
  if (stim >= isi) return 'Stimulus duration must be shorter than the stimulus interval.';
  if (intVal('missCriterion', 7) < 1) return 'The miss criterion must be at least 1.';
  if (intVal('maxMinutes', 40) < 1) return 'Maximum duration must be at least 1 minute.';
  if (intVal('lapseMs', 500) >= stim) return 'The lapse threshold must be below the hit window (stimulus duration).';
  return null;
}

async function startCalibration() {
  const err = validateSetup();
  if (err) { $('setupError').textContent = err; $('setupError').hidden = false; return; }
  $('setupError').hidden = true;

  show('screen-calibrate');
  $('calStatus').textContent = 'Measuring frame timing…';
  $('calReport').hidden = true;
  $('calProbe').hidden = true;
  $('btnBegin').hidden = true;

  await new Promise((r) => setTimeout(r, 250));

  calibration = await T.calibrateDisplay(intVal('calFrames', 240), 10);
  calGrade = T.gradeCalibration(calibration, displayInfo ? displayInfo.current.displayFrequency : null);
  renderCalReport();
  runInputProbe();
}

function renderCalReport() {
  const c = calibration;
  $('calStatus').textContent = '';
  $('calReport').hidden = false;
  $('calHz').textContent = c.refreshHz ? c.refreshHz.toFixed(2) + ' Hz' : '—';
  $('calInterval').textContent = n2(c.frameIntervalMs) + ' ms';
  $('calJitter').textContent = n2(c.madIntervalMs) + ' ms';
  $('calRange').textContent = n2(c.p05IntervalMs) + ' – ' + n2(c.p95IntervalMs) + ' ms';
  $('calDropped').textContent = c.droppedFrames +
    (c.droppedRate != null ? ' (' + (c.droppedRate * 100).toFixed(2) + '%)' : '');

  const badge = $('calGrade');
  badge.textContent = calGrade.grade;
  badge.className = 'badge badge-' + calGrade.grade;

  const list = $('calProblems');
  list.innerHTML = '';
  if (!calGrade.problems.length) {
    const li = document.createElement('li');
    li.textContent = 'Frame timing is stable. Onset times are quantised to ' +
      n2(c.frameIntervalMs) + ' ms, the physical limit of this display.';
    list.appendChild(li);
  } else {
    calGrade.problems.forEach((p) => {
      const li = document.createElement('li');
      li.className = 'warn';
      li.textContent = p;
      list.appendChild(li);
    });
  }
}

let probe = null, probeCount = 0;

function runInputProbe() {
  probe = T.makeInputProbe();
  probeCount = 0;
  inputProbeResult = null;
  $('calProbe').hidden = false;
  $('probeStatus').textContent = 'Press the space bar ' + PROBE_TAPS + ' more times.';
  document.addEventListener('keydown', probeKey);
}

function probeKey(e) {
  if (e.code !== 'Space') return;
  e.preventDefault();
  probe.record(e);
  probeCount += 1;
  const left = PROBE_TAPS - probeCount;
  if (left > 0) {
    $('probeStatus').textContent = 'Press the space bar ' + left + ' more time' + (left === 1 ? '' : 's') + '.';
    return;
  }
  document.removeEventListener('keydown', probeKey);
  inputProbeResult = probe.result();
  useEventStamp = inputProbeResult.usable;

  if (inputProbeResult.usable) {
    $('probeStatus').innerHTML =
      'Input dispatch delay: <strong>' + n2(inputProbeResult.medianDelayMs) + ' ms</strong> (median). ' +
      'Measured and removed — responses are timed from the OS event stamp, not from when this app handled them.';
    $('probeStatus').className = 'ok';
  } else {
    $('probeStatus').textContent =
      'Could not verify the input clock (' + inputProbeResult.reason + '). Falling back to handler time; ' +
      'reaction times carry a few ms of extra dispatch noise. This is recorded in the data.';
    $('probeStatus').className = 'warn';
  }
  $('btnBegin').hidden = false;
}

/* ---------------- trial ---------------- */

function beginTrial() {
  const isiMs = intVal('isiMs', 3000);
  const stimMs = intVal('stimMs', 1000);
  const corr = $('correction').value;
  const now = new Date();

  framesPerEpoch = Math.max(1, Math.round(isiMs / calibration.frameIntervalMs));
  framesStimOn = Math.max(1, Math.round(stimMs / calibration.frameIntervalMs));
  if (framesStimOn >= framesPerEpoch) framesStimOn = framesPerEpoch - 1;
  achievedIsiMs = framesPerEpoch * calibration.frameIntervalMs;

  cfg = {
    isiMs,
    stimMs,
    /* The hit window is the NOMINAL stimulus duration, not the frame-quantised
     * one, so classification is identical across machines with different
     * refresh rates. Only the light's actual duration is quantised. */
    hitWindowMs: stimMs,
    lapseMs: intVal('lapseMs', 500),
    missCriterion: intVal('missCriterion', 7),
    maxMinutes: intVal('maxMinutes', 40),
    maxMs: intVal('maxMinutes', 40) * 60000,
    correction: corr,
    removeAnticipations: corr === 'anticipations' || corr === 'both',
    removeOutliers: corr === 'outliers' || corr === 'both',
    anticipationMs: numVal('anticipationMs', 100),
    sdMultiplier: numVal('sdMultiplier', 2),
    presentationOffsetFrames: numVal('presOffsetFrames', 1),
    hardwareOffsetMs: numVal('hwOffsetMs', 0),
    photodiode: $('photodiode').checked,
    framesPerEpoch,
    framesStimOn,
    achievedIsiMs,
    achievedStimMs: framesStimOn * calibration.frameIntervalMs
  };
  maxEpochs = Math.floor(cfg.maxMs / achievedIsiMs);

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
  droppedFrames = 0;
  frameIndex = -1;
  lastFrameTs = null;

  $('patch').hidden = !cfg.photodiode;
  window.bsrt.preventDisplaySleep(true);
  if ($('useFullscreen').checked) window.bsrt.setFullscreen(true);
  countdown();
}

function countdown() {
  show('screen-countdown');
  let n = 3;
  $('countdownNum').textContent = n;
  const iv = setInterval(() => {
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
  trialStartTs = null;
  rafId = requestAnimationFrame(frameLoop);
}

function frameLoop(ts) {
  if (!running) return;
  rafId = requestAnimationFrame(frameLoop);

  if (trialStartTs === null) trialStartTs = ts;

  if (lastFrameTs !== null) {
    const mult = Math.round((ts - lastFrameTs) / calibration.frameIntervalMs);
    if (mult >= 2) droppedFrames += mult - 1;
  }
  lastFrameTs = ts;

  frameIndex += 1;
  const epochIdx = Math.floor(frameIndex / framesPerEpoch);
  const phase = frameIndex - epochIdx * framesPerEpoch;

  if (phase === 0) {
    finalizeEpoch();
    if (!running) return;
    if (epochIdx >= maxEpochs) { endTask('max_duration'); return; }
    startEpoch(epochIdx, ts);
  } else if (phase === framesStimOn) {
    ledOff();
  } else if (currentEpoch && currentEpoch.extinguished) {
    // Extinguish here rather than in the event handler so the offset lands on a
    // frame boundary like every other transition — keeps a photodiode trace
    // interpretable.
    ledOff();
    currentEpoch.extinguished = false;
  }
}

function startEpoch(idx, frameTs) {
  ledOn();
  currentEpoch = {
    index: idx,
    frameTs,
    presentationMs: frameTs + cfg.presentationOffsetFrames * calibration.frameIntervalMs,
    onsetMs: frameTs - trialStartTs,
    responded: false,
    rtRawMs: null,
    extra: 0,
    extinguished: false
  };
  epochs.push(currentEpoch);
}

/*
 * An epoch counts toward the sleep-onset criterion unless it was a HIT.
 * A response slower than the hit window is a MISS, so it does NOT reset the
 * consecutive-miss run even though a key was pressed.
 */
function finalizeEpoch() {
  if (!currentEpoch) return;
  const ep = currentEpoch;
  currentEpoch = null;

  const rt = ep.rtRawMs === null ? null : ep.rtRawMs - cfg.hardwareOffsetMs;
  const isHit = rt !== null && rt <= cfg.hitWindowMs;
  if (isHit) { consecutiveMisses = 0; missRunStartIndex = -1; return; }

  if (consecutiveMisses === 0) missRunStartIndex = ep.index;
  consecutiveMisses += 1;
  if (consecutiveMisses >= cfg.missCriterion) endTask('sleep_onset');
}

/*
 * Every response is recorded with its raw RT, however late. Hit/miss
 * classification happens at scoring time, not here.
 */
function handleResponse(evt) {
  if (!running || !currentEpoch) return;
  if (currentEpoch.responded) { currentEpoch.extra += 1; return; }

  const t = T.eventTime(evt, useEventStamp);
  currentEpoch.responded = true;
  currentEpoch.rtRawMs = t - currentEpoch.presentationMs;

  // Extinguishes the light, but never shortens the epoch.
  currentEpoch.extinguished = true;
}

function ledOn() {
  $('led').classList.add('on');
  if (cfg.photodiode) $('patch').classList.add('on');
}

function ledOff() {
  $('led').classList.remove('on');
  $('patch').classList.remove('on');
}

function endTask(reason) {
  running = false;
  if (rafId) cancelAnimationFrame(rafId);
  rafId = null;
  ledOff();
  $('patch').hidden = true;
  window.bsrt.preventDisplaySleep(false);
  window.bsrt.setFullscreen(false);

  lastResult = buildResult(reason);
  saveSession(lastResult);
  renderResult(lastResult);
  show('screen-results');
}

/* ---------------- result assembly ---------------- */

function buildResult(reason) {
  /*
   * The measured hardware offset is applied BEFORE scoring, because it is a
   * measurement correction rather than a statistical exclusion: with it, rtMs
   * is the participant's actual reaction time, so hit/miss classification uses
   * the right number. Both raw and offset-applied values are exported.
   */
  const scored = S.score(epochs.map((e) => ({
    index: e.index,
    onsetMs: e.onsetMs,
    rtMs: e.rtRawMs === null ? null : e.rtRawMs - cfg.hardwareOffsetMs
  })), cfg);

  let sleepOnsetMs = null, sleepOnsetCriterionMs = null;
  if (reason === 'sleep_onset' && missRunStartIndex >= 0) {
    const first = epochs[missRunStartIndex];
    const last = epochs[missRunStartIndex + cfg.missCriterion - 1];
    sleepOnsetMs = first ? first.onsetMs : missRunStartIndex * achievedIsiMs;
    sleepOnsetCriterionMs = last ? last.onsetMs + achievedIsiMs
                                 : sleepOnsetMs + cfg.missCriterion * achievedIsiMs;
  }

  return {
    runId: meta.runId,
    startedAt: startWall,
    participant: meta,
    config: cfg,
    display: displayInfo,
    calibration,
    calibrationGrade: calGrade ? calGrade.grade : null,
    inputProbe: inputProbeResult,
    inputTimeSource: useEventStamp ? 'os_event_stamp' : 'handler_time',
    endReason: reason,
    sleptBeforeMax: reason === 'sleep_onset',
    sleepOnsetMs,
    sleepOnsetCriterionMs,
    elapsedMs: epochs.length ? epochs[epochs.length - 1].onsetMs + achievedIsiMs : 0,
    extraResponses: epochs.reduce((n, e) => n + e.extra, 0),
    extras: epochs.map((e) => e.extra),
    rawRts: epochs.map((e) => e.rtRawMs),
    droppedFramesDuringTrial: droppedFrames,
    blurCount,
    scored
  };
}

/* ---------------- rendering ---------------- */

function renderResult(r) {
  const p = r.participant, t = r.scored.totals, e = r.scored.errorProfiles;

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

  $('tRefresh').textContent = r.calibration.refreshHz.toFixed(2) + ' Hz';
  $('tInterval').textContent = n2(r.calibration.frameIntervalMs) + ' ms';
  $('tAchievedIsi').textContent = n2(r.config.achievedIsiMs) + ' ms (requested ' + r.config.isiMs + ')';
  $('tQuant').textContent = '± ' + n2(r.calibration.frameIntervalMs / 2) + ' ms';
  $('tDropped').textContent = r.droppedFramesDuringTrial;
  $('tInputSrc').textContent = r.inputTimeSource === 'os_event_stamp' ? 'OS event stamp' : 'handler time (fallback)';
  $('tInputDelay').textContent = r.inputProbe && r.inputProbe.usable
    ? n2(r.inputProbe.medianDelayMs) + ' ms (removed)' : '—';
  $('tHwOffset').textContent = r.config.hardwareOffsetMs
    ? r.config.hardwareOffsetMs + ' ms (applied before scoring)' : 'not measured';
}

function describeCorrection(c) {
  if (c.correction === 'anticipations') return 'RT < ' + c.anticipationMs + ' ms removed';
  if (c.correction === 'outliers') return 'RT beyond ' + c.sdMultiplier + ' SD removed';
  if (c.correction === 'both') return 'RT < ' + c.anticipationMs + ' ms removed, then RT beyond ' + c.sdMultiplier + ' SD';
  return 'none';
}

function renderPerMinute(r) {
  const pm = r.scored.perMinute;
  const vel = r.scored.dynamicsRs.velocity;
  const acc = r.scored.dynamicsRs.acceleration;
  const body = $('pmBody');
  body.innerHTML = '';

  pm.forEach((m, i) => {
    const tr = document.createElement('tr');
    [m.minute, m.trials, m.misses, m.lapses,
     n1(m.avgRt), n1(m.medianRt), n1(m.sdRt), n1(m.fastest10Rt), n1(m.slowest10Rt),
     n3(m.avgRs), n3(m.medianRs),
     vel[i] == null ? '—' : n2(vel[i]),
     acc[i] == null ? '—' : n2(acc[i])
    ].forEach((v) => {
      const td = document.createElement('td');
      td.textContent = v;
      tr.appendChild(td);
    });
    body.appendChild(tr);
  });

  $('pmSlope').textContent = r.scored.dynamicsRs.slope == null
    ? 'Not enough minutes to fit a trend.'
    : 'Reaction-speed trend across the test: ' + n3(r.scored.dynamicsRs.slope) +
      ' units per minute (negative = slowing down).';
}

/* ---------------- exports ---------------- */

const PARTICIPANT_COLS = ['run_id', 'date', 'time', 'participant_id', 'name', 'address',
                          'birth_date', 'educational_level', 'session_label', 'trial_number'];

function participantVals(r) {
  const p = r.participant;
  return [r.runId, p.date, p.time, p.participantId, p.name, p.address,
          p.birthDate, p.education, p.sessionLabel, p.trialNumber];
}

const RAW_HEADER = PARTICIPANT_COLS.concat([
  'epoch_index', 'minute', 'onset_ms', 'responded',
  'rt_raw_ms', 'rt_ms', 'rs_per_sec',
  'outcome', 'lapse', 'late_response', 'anticipation', 'extra_responses'
]);

function rawRows(r) {
  const pv = participantVals(r);
  const extras = r.extras || [];
  const raws = r.rawRts || [];
  return r.scored.trials.map((t, i) => pv.concat([
    t.index, t.minute + 1, round(t.onsetMs, 2),
    t.rtMs === null ? 0 : 1,
    round(raws[i] == null ? null : raws[i], 3),   // before the hardware offset
    round(t.rtMs, 3),                              // after it — what was scored
    round(t.rsPerSec, 5),
    t.outcome, t.lapse, t.lateResponse, t.anticipation,
    extras[i] == null ? '' : extras[i]
  ]));
}

const PM_HEADER = PARTICIPANT_COLS.concat([
  'minute', 'trials', 'hits', 'misses', 'lapses', 'late_responses', 'hit_ratio',
  'avg_rt', 'median_rt', 'stdev_rt', 'fastest10_rt', 'slowest10_rt',
  'corr_avg_rt', 'corr_median_rt', 'corr_stdev_rt', 'corr_fastest10_rt', 'corr_slowest10_rt',
  'avg_rs', 'median_rs', 'stdev_rs', 'fastest10_rs', 'slowest10_rs',
  'corr_avg_rs', 'corr_median_rs', 'corr_stdev_rs', 'corr_fastest10_rs', 'corr_slowest10_rs',
  'n_rt', 'n_rt_corrected', 'anticipations_removed', 'outliers_removed',
  'velocity_rs', 'acceleration_rs', 'velocity_rt', 'acceleration_rt'
]);

function pmRows(r) {
  const pv = participantVals(r);
  const vRs = r.scored.dynamicsRs, vRt = r.scored.dynamicsRt;
  return r.scored.perMinute.map((m, i) => pv.concat([
    m.minute, m.trials, m.hits, m.misses, m.lapses, m.lateResponses, round(m.hitRatio, 4),
    round(m.avgRt, 2), round(m.medianRt, 2), round(m.sdRt, 2), round(m.fastest10Rt, 2), round(m.slowest10Rt, 2),
    round(m.corrAvgRt, 2), round(m.corrMedianRt, 2), round(m.corrSdRt, 2), round(m.corrFastest10Rt, 2), round(m.corrSlowest10Rt, 2),
    round(m.avgRs, 5), round(m.medianRs, 5), round(m.sdRs, 5), round(m.fastest10Rs, 5), round(m.slowest10Rs, 5),
    round(m.corrAvgRs, 5), round(m.corrMedianRs, 5), round(m.corrSdRs, 5), round(m.corrFastest10Rs, 5), round(m.corrSlowest10Rs, 5),
    m.n, m.nCorrected, m.nAnticipationsRemoved, m.nOutliersRemoved,
    round(vRs.velocity[i], 5), round(vRs.acceleration[i], 5),
    round(vRt.velocity[i], 3), round(vRt.acceleration[i], 3)
  ]));
}

const SUMMARY_HEADER = PARTICIPANT_COLS.concat([
  'sleep_onset_ms', 'sleep_onset_criterion_ms', 'slept_before_max', 'end_reason', 'elapsed_ms',
  'total_trialrun', 'hit_ratio', 'total_hits', 'total_miss', 'total_lapse', 'late_responses',
  'ep_1_2', 'ep_3_6', 'ep_7plus', 'longest_miss_run',
  'total_avg_rt', 'total_median_rt', 'total_stdev_rt', 'fastest10_rt', 'slowest10_rt',
  'corr_total_avg_rt', 'corr_total_median_rt', 'corr_total_stdev_rt', 'corr_fastest10_rt', 'corr_slowest10_rt',
  'total_avg_rs', 'total_median_rs', 'total_stdev_rs', 'fastest10_rs', 'slowest10_rs',
  'corr_total_avg_rs', 'corr_total_median_rs', 'corr_total_stdev_rs', 'corr_fastest10_rs', 'corr_slowest10_rs',
  'n_rt', 'n_rt_corrected', 'anticipations_removed', 'outliers_removed', 'rs_undefined',
  'rs_slope_per_min', 'rt_slope_per_min',
  'isi_requested_ms', 'isi_achieved_ms', 'stim_requested_ms', 'stim_achieved_ms',
  'hit_window_ms', 'lapse_threshold_ms', 'miss_criterion', 'max_minutes',
  'correction', 'anticipation_threshold_ms', 'sd_multiplier',
  'refresh_hz_measured', 'refresh_hz_reported', 'frame_interval_ms', 'frame_jitter_mad_ms',
  'onset_quantisation_ms', 'dropped_frames_calibration', 'dropped_frames_trial', 'calibration_grade',
  'input_time_source', 'input_dispatch_median_ms',
  'presentation_offset_frames', 'hardware_offset_ms',
  'platform', 'os_release', 'electron_version', 'chrome_version', 'display_scale_factor',
  'extra_responses', 'page_blur_count'
]);

function summaryRow(r) {
  const t = r.scored.totals, e = r.scored.errorProfiles, c = r.config;
  const cal = r.calibration, d = r.display ? r.display.current : {};
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
    c.isiMs, round(c.achievedIsiMs, 3), c.stimMs, round(c.achievedStimMs, 3),
    c.hitWindowMs, c.lapseMs, c.missCriterion, c.maxMinutes,
    c.correction, c.anticipationMs, c.sdMultiplier,
    round(cal.refreshHz, 3), d.displayFrequency == null ? null : d.displayFrequency,
    round(cal.frameIntervalMs, 4), round(cal.madIntervalMs, 4),
    round(cal.frameIntervalMs / 2, 3), cal.droppedFrames, r.droppedFramesDuringTrial, r.calibrationGrade,
    r.inputTimeSource, r.inputProbe && r.inputProbe.usable ? round(r.inputProbe.medianDelayMs, 3) : null,
    c.presentationOffsetFrames, c.hardwareOffsetMs,
    r.display ? r.display.platform : null, r.display ? r.display.osRelease : null,
    r.display ? r.display.electron : null, r.display ? r.display.chrome : null,
    d.scaleFactor == null ? null : d.scaleFactor,
    r.extraResponses, r.blurCount
  ]);
}

function bulk(builder, header, name) {
  const all = loadSessions();
  if (!all.length) { alert('No saved trials to export.'); return; }
  let rows = [];
  all.forEach((r) => { rows = rows.concat(builder(r)); });
  saveCsv(name, toCsv(header, rows));
}

/* ---------------- wiring ---------------- */

$('btnStart').addEventListener('click', startCalibration);
$('btnBegin').addEventListener('click', beginTrial);
$('btnCalBack').addEventListener('click', () => show('screen-setup'));

$('btnAgain').addEventListener('click', () => {
  $('trialNumber').value = intVal('trialNumber', 1) + 1;
  show('screen-setup');
});

$('btnExportRaw').addEventListener('click', () => {
  if (lastResult) saveCsv(fileStem(lastResult) + '_raw.csv', toCsv(RAW_HEADER, rawRows(lastResult)));
});
$('btnExportPm').addEventListener('click', () => {
  if (lastResult) saveCsv(fileStem(lastResult) + '_perminute.csv', toCsv(PM_HEADER, pmRows(lastResult)));
});
$('btnExportSummary').addEventListener('click', () => {
  if (lastResult) saveCsv(fileStem(lastResult) + '_summary.csv', toCsv(SUMMARY_HEADER, [summaryRow(lastResult)]));
});

$('btnExportAllRaw').addEventListener('click', () => bulk(rawRows, RAW_HEADER, 'bsrt_all_raw.csv'));
$('btnExportAllPm').addEventListener('click', () => bulk(pmRows, PM_HEADER, 'bsrt_all_perminute.csv'));
$('btnExportAllSummary').addEventListener('click', () => {
  const all = loadSessions();
  if (!all.length) { alert('No saved trials to export.'); return; }
  saveCsv('bsrt_all_summaries.csv', toCsv(SUMMARY_HEADER, all.map(summaryRow)));
});

$('btnAbort').addEventListener('click', (e) => {
  e.stopPropagation();
  if (running) endTask('aborted');
});

$('screen-task').addEventListener('pointerdown', (e) => {
  if (e.target && e.target.id === 'btnAbort') return;
  handleResponse(e);
});

document.addEventListener('keydown', (e) => {
  if (!running) return;
  if (e.code === 'Space' || e.code === 'Enter') { e.preventDefault(); handleResponse(e); }
});

document.addEventListener('visibilitychange', () => {
  if (running && document.hidden) blurCount += 1;
});

$('maxMinutes').addEventListener('change', () => {
  const short = intVal('maxMinutes', 40) < 5;
  $('lapseHint').textContent = short
    ? 'Short protocol (< 5 min): 355 ms is the conventional lapse threshold.'
    : 'Regular protocol: 500 ms is the conventional lapse threshold.';
});

init();
