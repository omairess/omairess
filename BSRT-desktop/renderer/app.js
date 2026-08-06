'use strict';

/* BSRT desktop — task logic.
 *
 * Differs from the browser version in one fundamental way: stimuli are
 * scheduled by counting real vsync intervals rather than by setTimeout, and
 * onset times come from the frame clock rather than from "whenever the class
 * got added". See timing.js and README.md.
 */

const T = window.BSRTTiming;
const LS_KEY = 'bsrt.desktop.sessions.v1';
/* Conventional PVT lapse thresholds are 500 ms; OSLER work often uses 1000 ms;
 * short protocols need tighter cuts (355 ms has been used for 3-minute tests).
 * Per-trial setting rather than a constant. */
const DEFAULT_LAPSE_MS = 1000;
const PROBE_TAPS = 5;

/*
 * Responses faster than this are anticipations, not reactions — the participant
 * pressed before they could have processed the stimulus. Standard PVT practice
 * is to count them separately and keep them out of the reaction-time
 * distribution. Negative values land here too: they mean the response beat our
 * estimate of when the frame reached the panel, which is a sign the
 * presentation-pipeline assumption is off for this machine.
 */
const ANTICIPATION_MS = 100;

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
let longestMissRun = 0;
let missRunStartIndex = -1;
let blurCount = 0;
let lastResult = null;

/* ---------------- helpers ---------------- */

const $ = (id) => document.getElementById(id);

function show(id) {
  document.querySelectorAll('.screen').forEach((s) => s.classList.remove('active'));
  $(id).classList.add('active');
}

function intVal(id, fb) {
  const n = parseInt($(id).value, 10);
  return isFinite(n) ? n : fb;
}

function numVal(id, fb) {
  const n = parseFloat($(id).value);
  return isFinite(n) ? n : fb;
}

function fmtClock(ms) {
  if (ms == null) return '—';
  const t = Math.round(ms / 1000);
  const m = Math.floor(t / 60);
  const s = t % 60;
  return m + ':' + (s < 10 ? '0' : '') + s;
}

const fmtMs = (v) => (v == null ? '—' : Math.round(v) + ' ms');
const fmt2 = (v) => (v == null ? '—' : v.toFixed(2));

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

/* ---------------- storage ---------------- */

function loadSessions() {
  try {
    const a = JSON.parse(localStorage.getItem(LS_KEY) || '[]');
    return Array.isArray(a) ? a : [];
  } catch (e) {
    return [];
  }
}

function saveSession(r) {
  const all = loadSessions();
  all.push(r);
  try {
    localStorage.setItem(LS_KEY, JSON.stringify(all));
  } catch (e) {
    alert('Could not save to local storage. Export this trial now.');
  }
  $('sessionCount').textContent = loadSessions().length;
}

/* ---------------- startup ---------------- */

async function init() {
  $('sessionCount').textContent = loadSessions().length;
  try {
    displayInfo = await window.bsrt.getDisplayInfo();
  } catch (e) {
    displayInfo = null;
  }
  renderEnvironment();
}

function renderEnvironment() {
  if (!displayInfo) {
    $('envInfo').textContent = 'Display information unavailable.';
    return;
  }
  const d = displayInfo.current;
  const bits = [
    d.displayFrequency ? d.displayFrequency + ' Hz (reported)' : 'refresh rate unknown',
    d.bounds.width + '×' + d.bounds.height,
    'scale ' + d.scaleFactor + '×',
    displayInfo.platform === 'darwin' ? 'macOS' : displayInfo.platform === 'win32' ? 'Windows' : displayInfo.platform,
    'Electron ' + displayInfo.electron
  ];
  $('envInfo').textContent = bits.join(' · ');
  if (displayInfo.displayCount > 1) {
    $('envWarn').hidden = false;
    $('envWarn').textContent =
      displayInfo.displayCount + ' displays detected. Run the task on one screen only — ' +
      'mirroring or extending can force the compositor to the slowest display.';
  }
}

/* ---------------- calibration flow ---------------- */

async function startCalibration() {
  const err = validateSetup();
  if (err) {
    $('setupError').textContent = err;
    $('setupError').hidden = false;
    return;
  }
  $('setupError').hidden = true;

  show('screen-calibrate');
  $('calStatus').textContent = 'Measuring frame timing…';
  $('calReport').hidden = true;
  $('calProbe').hidden = true;
  $('btnBegin').hidden = true;

  // Let the screen settle before timing it.
  await new Promise((r) => setTimeout(r, 250));

  const frames = intVal('calFrames', 240);
  calibration = await T.calibrateDisplay(frames, 10);
  calGrade = T.gradeCalibration(calibration, displayInfo ? displayInfo.current.displayFrequency : null);

  renderCalReport();
  runInputProbe();
}

function renderCalReport() {
  const c = calibration;
  $('calStatus').textContent = '';
  $('calReport').hidden = false;

  $('calHz').textContent = c.refreshHz ? c.refreshHz.toFixed(2) + ' Hz' : '—';
  $('calInterval').textContent = fmt2(c.frameIntervalMs) + ' ms';
  $('calJitter').textContent = fmt2(c.madIntervalMs) + ' ms';
  $('calRange').textContent = fmt2(c.p05IntervalMs) + ' – ' + fmt2(c.p95IntervalMs) + ' ms';
  $('calDropped').textContent =
    c.droppedFrames + (c.droppedRate != null ? ' (' + (c.droppedRate * 100).toFixed(2) + '%)' : '');

  const badge = $('calGrade');
  badge.textContent = calGrade.grade;
  badge.className = 'badge badge-' + calGrade.grade;

  const list = $('calProblems');
  list.innerHTML = '';
  if (!calGrade.problems.length) {
    const li = document.createElement('li');
    li.textContent = 'Frame timing is stable. Onset times are quantised to ' +
      fmt2(c.frameIntervalMs) + ' ms, which is the physical limit of this display.';
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

let probe = null;
let probeCount = 0;

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
      'Input dispatch delay: <strong>' + fmt2(inputProbeResult.medianDelayMs) + ' ms</strong> (median). ' +
      'This delay is measured and removed — responses are timed from the OS event stamp, not from when ' +
      'this app got round to handling them.';
    $('probeStatus').className = 'ok';
  } else {
    $('probeStatus').textContent =
      'Could not verify the input clock (' + inputProbeResult.reason + '). ' +
      'Falling back to handler time; reaction times will carry a few ms of extra dispatch noise. ' +
      'This is recorded in the data.';
    $('probeStatus').className = 'warn';
  }
  $('btnBegin').hidden = false;
}

/* ---------------- setup validation ---------------- */

function validateSetup() {
  const isi = intVal('isiMs', 3000);
  const stim = intVal('stimMs', 1000);
  const crit = intVal('missCriterion', 7);
  const maxMin = intVal('maxMinutes', 40);
  if (isi < 200) return 'Stimulus interval must be at least 200 ms.';
  if (stim >= isi) return 'Stimulus duration must be shorter than the stimulus interval.';
  if (crit < 1) return 'The miss criterion must be at least 1.';
  if (maxMin < 1) return 'Maximum duration must be at least 1 minute.';
  return null;
}

/* ---------------- trial ---------------- */

function beginTrial() {
  const isiMs = intVal('isiMs', 3000);
  const stimMs = intVal('stimMs', 1000);

  framesPerEpoch = Math.max(1, Math.round(isiMs / calibration.frameIntervalMs));
  framesStimOn = Math.max(1, Math.round(stimMs / calibration.frameIntervalMs));
  if (framesStimOn >= framesPerEpoch) framesStimOn = framesPerEpoch - 1;
  achievedIsiMs = framesPerEpoch * calibration.frameIntervalMs;

  cfg = {
    isiMs,
    stimMs,
    missCriterion: intVal('missCriterion', 7),
    maxMinutes: intVal('maxMinutes', 40),
    maxMs: intVal('maxMinutes', 40) * 60000,
    presentationOffsetFrames: numVal('presOffsetFrames', 1),
    hardwareOffsetMs: numVal('hwOffsetMs', 0),
    photodiode: $('photodiode').checked,
    lapseMs: intVal('lapseMs', DEFAULT_LAPSE_MS),
    responseWindow: $('responseWindow').value,
    framesPerEpoch,
    framesStimOn,
    achievedIsiMs,
    achievedStimMs: framesStimOn * calibration.frameIntervalMs
  };
  maxEpochs = Math.floor(cfg.maxMs / achievedIsiMs);

  meta = {
    runId: 'bsrt-' + Date.now(),
    participantId: $('participantId').value.trim() || 'NA',
    sessionLabel: $('sessionLabel').value.trim() || 'NA',
    trialNumber: intVal('trialNumber', 1)
  };

  epochs = [];
  currentEpoch = null;
  consecutiveMisses = 0;
  longestMissRun = 0;
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
  trialStartTs = null;
  rafId = requestAnimationFrame(frameLoop);
}

/*
 * One loop, driven by the display. Every epoch boundary is a counted frame,
 * so presentation is periodic by construction: it cannot drift relative to the
 * screen, which is the thing the participant actually sees.
 */
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
    if (epochIdx >= maxEpochs) {
      endTask('max_duration');
      return;
    }
    startEpoch(epochIdx, ts);
  } else if (phase === framesStimOn) {
    ledOff();
  } else if (currentEpoch && currentEpoch.extinguished) {
    // A response came in since the last frame. Extinguish here rather than in
    // the event handler so the offset lands on a frame boundary like every
    // other transition — it keeps the photodiode trace interpretable.
    ledOff();
    currentEpoch.extinguished = false;
  }
}

function startEpoch(idx, frameTs) {
  ledOn();
  currentEpoch = {
    index: idx,
    frameTs,
    // The paint issued in this callback reaches the panel roughly one frame
    // later. That pipeline depth is an assumption, exposed as a setting and
    // recorded in the data, because only a photodiode can measure it.
    presentationMs: frameTs + cfg.presentationOffsetFrames * calibration.frameIntervalMs,
    onsetMs: frameTs - trialStartTs,
    responded: false,
    rtRawMs: null,
    extra: 0,
    late: 0,
    extinguished: false
  };
  epochs.push(currentEpoch);
}

function finalizeEpoch() {
  if (!currentEpoch) return;
  const ep = currentEpoch;
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

function handleResponse(evt) {
  if (!running || !currentEpoch) return;

  const t = T.eventTime(evt, useEventStamp);
  const elapsed = t - currentEpoch.presentationMs;

  // Under the 'stimulus' response window, a press after the light has gone out
  // is not a response to that stimulus, so the epoch stays a miss.
  if (cfg.responseWindow === 'stimulus' && elapsed > cfg.achievedStimMs) {
    currentEpoch.late += 1;
    return;
  }

  if (currentEpoch.responded) {
    currentEpoch.extra += 1;
    return;
  }
  currentEpoch.responded = true;
  currentEpoch.rtRawMs = elapsed;

  // Responding extinguishes the stimulus, but never shortens the epoch: the
  // next stimulus still starts a fixed number of frames after this one's onset,
  // so total task duration does not depend on how the participant responds.
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

/* ---------------- scoring ---------------- */

function buildResult(reason) {
  const rts = [];
  let hits = 0, misses = 0, extra = 0;
  for (const ep of epochs) {
    extra += ep.extra;
    if (ep.responded) { hits += 1; rts.push(ep.rtRawMs); } else { misses += 1; }
  }
  // Split anticipations out before summarising. Raw per-epoch RTs are exported
  // untouched, so any other exclusion rule can be applied downstream in R.
  const valid = rts.filter((v) => v >= ANTICIPATION_MS);
  const anticipations = rts.length - valid.length;
  const corrected = valid.map((v) => v - cfg.hardwareOffsetMs);
  const lapses = valid.filter((v) => v > cfg.lapseMs).length;
  const late = epochs.reduce((n, e) => n + e.late, 0);

  /*
   * Latency comes from the frame timestamps we actually observed, not from
   * `index x achievedIsiMs`. The calibrated frame interval carries a small
   * estimation error, and multiplying it by an epoch index would let that
   * error accumulate across a 40-minute trial. Measured onsets cannot drift.
   * The nominal arithmetic is kept only as a fallback.
   */
  let latencyFirstMissMs = null;
  let latencyCriterionMs = null;
  if (reason === 'sleep_onset' && missRunStartIndex >= 0) {
    const first = epochs[missRunStartIndex];
    const last = epochs[missRunStartIndex + cfg.missCriterion - 1];
    latencyFirstMissMs = first ? first.onsetMs : missRunStartIndex * achievedIsiMs;
    latencyCriterionMs = last
      ? last.onsetMs + achievedIsiMs
      : latencyFirstMissMs + cfg.missCriterion * achievedIsiMs;
  }

  return {
    runId: meta.runId,
    startedAt: startWall,
    participantId: meta.participantId,
    sessionLabel: meta.sessionLabel,
    trialNumber: meta.trialNumber,
    config: cfg,
    display: displayInfo,
    calibration,
    calibrationGrade: calGrade ? calGrade.grade : null,
    inputProbe: inputProbeResult,
    inputTimeSource: useEventStamp ? 'os_event_stamp' : 'handler_time',
    endReason: reason,
    sleptBeforeMax: reason === 'sleep_onset',
    latencyFirstMissMs,
    latencyCriterionMs,
    elapsedMs: epochs.length ? epochs[epochs.length - 1].onsetMs + achievedIsiMs : 0,
    nEpochs: epochs.length,
    hits,
    misses,
    hitRate: epochs.length ? hits / epochs.length : null,
    longestMissRun,
    anticipations,
    nRtValid: valid.length,
    anticipationThresholdMs: ANTICIPATION_MS,
    meanRtMs: T.stats.mean(valid),
    medianRtMs: T.stats.median(valid),
    sdRtMs: T.stats.sd(valid),
    madRtMs: T.stats.mad(valid),
    meanRtCorrectedMs: T.stats.mean(corrected),
    medianRtCorrectedMs: T.stats.median(corrected),
    lapses,
    lapseThresholdMs: cfg.lapseMs,
    lateResponses: late,
    extraResponses: extra,
    droppedFramesDuringTrial: droppedFrames,
    blurCount,
    epochs: epochs.map((e) => ({
      index: e.index,
      onsetMs: Math.round(e.onsetMs * 100) / 100,
      responded: e.responded ? 1 : 0,
      rtRawMs: e.rtRawMs == null ? null : Math.round(e.rtRawMs * 100) / 100,
      lapse: e.rtRawMs == null ? null : (e.rtRawMs > cfg.lapseMs ? 1 : 0),
      extra: e.extra,
      late: e.late
    }))
  };
}

/* ---------------- results ---------------- */

function renderResult(r) {
  $('resultWho').textContent =
    r.participantId + ' · ' + r.sessionLabel + ' · trial ' + r.trialNumber;

  if (r.sleptBeforeMax) {
    $('mLatency').textContent = fmtClock(r.latencyFirstMissMs);
    $('mLatencyFoot').textContent =
      'to the first of ' + r.config.missCriterion + ' consecutive missed stimuli (criterion confirmed at ' +
      fmtClock(r.latencyCriterionMs) + ')';
    $('mOutcome').textContent = 'Sleep onset was scored before the ' + r.config.maxMinutes + '-minute ceiling.';
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
  $('mLate').textContent = r.config.responseWindow === 'stimulus'
    ? r.lateResponses + ' (after light-off — not counted)'
    : 'n/a (full-epoch window)';
  $('mAntic').textContent = r.anticipations +
    (r.anticipations ? ' (excluded from RT statistics)' : '');
  $('mExtra').textContent = r.extraResponses;

  // A large crop of anticipations usually means the presentation-pipeline
  // assumption is wrong for this machine, not that the participant is psychic.
  const anticRate = r.hits ? r.anticipations / r.hits : 0;
  const flag = $('anticFlag');
  if (anticRate > 0.1 && r.hits >= 10) {
    flag.hidden = false;
    flag.textContent =
      (anticRate * 100).toFixed(0) + '% of responses were faster than ' + ANTICIPATION_MS +
      ' ms. That is implausible for genuine reactions and suggests the presentation ' +
      'pipeline setting (' + r.config.presentationOffsetFrames + ' frame) overestimates ' +
      'this display’s lag. Reaction times are shifted; sleep latency is unaffected.';
  } else {
    flag.hidden = true;
  }

  $('tRefresh').textContent = r.calibration.refreshHz.toFixed(2) + ' Hz';
  $('tInterval').textContent = fmt2(r.calibration.frameIntervalMs) + ' ms';
  $('tAchievedIsi').textContent = fmt2(r.config.achievedIsiMs) + ' ms (requested ' + r.config.isiMs + ')';
  $('tQuant').textContent = '± ' + fmt2(r.calibration.frameIntervalMs / 2) + ' ms';
  $('tDropped').textContent = r.droppedFramesDuringTrial;
  $('tInputSrc').textContent =
    r.inputTimeSource === 'os_event_stamp' ? 'OS event stamp' : 'handler time (fallback)';
  $('tInputDelay').textContent =
    r.inputProbe && r.inputProbe.usable ? fmt2(r.inputProbe.medianDelayMs) + ' ms (removed)' : '—';
  $('tHwOffset').textContent = r.config.hardwareOffsetMs
    ? r.config.hardwareOffsetMs + ' ms (applied to corrected RT)'
    : 'not measured';
}

/* ---------------- export ---------------- */

const EPOCH_HEADER = [
  'run_id', 'participant_id', 'session_label', 'trial_number', 'started_at',
  'epoch_index', 'onset_ms', 'responded', 'rt_raw_ms', 'rt_corrected_ms', 'lapse',
  'extra_responses', 'late_responses'
];

const epochRows = (r) => r.epochs.map((e) => [
  r.runId, r.participantId, r.sessionLabel, r.trialNumber, r.startedAt,
  e.index, e.onsetMs, e.responded,
  e.rtRawMs,
  e.rtRawMs == null ? null : Math.round((e.rtRawMs - r.config.hardwareOffsetMs) * 100) / 100,
  e.lapse, e.extra, e.late
]);

const SUMMARY_HEADER = [
  'run_id', 'participant_id', 'session_label', 'trial_number', 'started_at',
  'isi_requested_ms', 'isi_achieved_ms', 'stim_requested_ms', 'stim_achieved_ms',
  'miss_criterion', 'max_minutes', 'lapse_threshold_ms', 'response_window',
  'end_reason', 'slept_before_max', 'latency_first_miss_ms', 'latency_criterion_ms', 'elapsed_ms',
  'n_epochs', 'hits', 'misses', 'hit_rate', 'longest_miss_run',
  'mean_rt_ms', 'median_rt_ms', 'sd_rt_ms', 'mad_rt_ms',
  'mean_rt_corrected_ms', 'median_rt_corrected_ms',
  'n_rt_valid', 'anticipations', 'anticipation_threshold_ms',
  'lapses', 'extra_responses', 'late_responses',
  'refresh_hz_measured', 'refresh_hz_reported', 'frame_interval_ms', 'frame_jitter_mad_ms',
  'onset_quantisation_ms', 'dropped_frames_calibration', 'dropped_frames_trial', 'calibration_grade',
  'input_time_source', 'input_dispatch_median_ms',
  'presentation_offset_frames', 'hardware_offset_ms',
  'platform', 'os_release', 'electron_version', 'chrome_version', 'display_scale_factor'
];

const round = (v, dp) => (v == null ? null : Number(v.toFixed(dp)));

function summaryRow(r) {
  const c = r.calibration;
  const d = r.display ? r.display.current : {};
  return [
    r.runId, r.participantId, r.sessionLabel, r.trialNumber, r.startedAt,
    r.config.isiMs, round(r.config.achievedIsiMs, 3), r.config.stimMs, round(r.config.achievedStimMs, 3),
    r.config.missCriterion, r.config.maxMinutes, r.config.lapseMs, r.config.responseWindow,
    r.endReason, r.sleptBeforeMax ? 1 : 0,
    round(r.latencyFirstMissMs, 1), round(r.latencyCriterionMs, 1), round(r.elapsedMs, 1),
    r.nEpochs, r.hits, r.misses, round(r.hitRate, 4), r.longestMissRun,
    round(r.meanRtMs, 2), round(r.medianRtMs, 2), round(r.sdRtMs, 2), round(r.madRtMs, 2),
    round(r.meanRtCorrectedMs, 2), round(r.medianRtCorrectedMs, 2),
    r.nRtValid, r.anticipations, r.anticipationThresholdMs,
    r.lapses, r.extraResponses, r.lateResponses,
    round(c.refreshHz, 3), r.display ? d.displayFrequency : null,
    round(c.frameIntervalMs, 4), round(c.madIntervalMs, 4),
    round(c.frameIntervalMs / 2, 3), c.droppedFrames, r.droppedFramesDuringTrial, r.calibrationGrade,
    r.inputTimeSource, r.inputProbe && r.inputProbe.usable ? round(r.inputProbe.medianDelayMs, 3) : null,
    r.config.presentationOffsetFrames, r.config.hardwareOffsetMs,
    r.display ? r.display.platform : null, r.display ? r.display.osRelease : null,
    r.display ? r.display.electron : null, r.display ? r.display.chrome : null,
    d.scaleFactor == null ? null : d.scaleFactor
  ];
}

/* ---------------- wiring ---------------- */

$('btnStart').addEventListener('click', startCalibration);
$('btnBegin').addEventListener('click', beginTrial);
$('btnCalBack').addEventListener('click', () => show('screen-setup'));

$('btnAgain').addEventListener('click', () => {
  $('trialNumber').value = intVal('trialNumber', 1) + 1;
  show('screen-setup');
});

$('btnExportEpochs').addEventListener('click', () => {
  if (!lastResult) return;
  const r = lastResult;
  saveCsv(
    'bsrt_' + safeName(r.participantId) + '_' + safeName(r.sessionLabel) + '_t' + r.trialNumber + '_epochs.csv',
    toCsv(EPOCH_HEADER, epochRows(r))
  );
});

$('btnExportSummary').addEventListener('click', () => {
  if (!lastResult) return;
  const r = lastResult;
  saveCsv(
    'bsrt_' + safeName(r.participantId) + '_' + safeName(r.sessionLabel) + '_t' + r.trialNumber + '_summary.csv',
    toCsv(SUMMARY_HEADER, [summaryRow(r)])
  );
});

$('btnExportAll').addEventListener('click', () => {
  const all = loadSessions();
  if (!all.length) return alert('No saved trials.');
  let rows = [];
  all.forEach((r) => { rows = rows.concat(epochRows(r)); });
  saveCsv('bsrt_all_epochs.csv', toCsv(EPOCH_HEADER, rows));
});

$('btnExportAllSummary').addEventListener('click', () => {
  const all = loadSessions();
  if (!all.length) return alert('No saved trials.');
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
  if (e.code === 'Space' || e.code === 'Enter') {
    e.preventDefault();
    handleResponse(e);
  }
});

document.addEventListener('visibilitychange', () => {
  if (running && document.hidden) blurCount += 1;
});

window.addEventListener('error', (e) => {
  if (running) {
    console.error('Runtime error during trial:', e.message);
  }
});

init();
