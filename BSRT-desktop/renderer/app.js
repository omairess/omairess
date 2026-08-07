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
const L = window.BSRTi18n;
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

let framesStimOn = 0;
/* With a variable schedule each epoch has its own frame count, so the loop
 * walks a list of frame boundaries instead of dividing by a constant. */
let achievedIsis = [];
let nextEpochIdx = 0;
let currentEpochStartFrame = 0;

let consecutiveMisses = 0;
let missRunStartIndex = -1;
let blurCount = 0;
let lastResult = null;
let pendingResult = null;
let presses = [];
let kssBefore = null;
let kssAfter = null;
let kssStage = null;
let audioCtx = null;
let alarmNodes = null;

/* ---------------- helpers ---------------- */

const $ = (id) => document.getElementById(id);

function show(id) {
  document.querySelectorAll('.screen').forEach((s) => s.classList.remove('active'));
  $(id).classList.add('active');
  // Lets CSS hide the floating language control on the dark task screens.
  document.body.setAttribute('data-screen', id);
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
  // A participant record is created or updated only once a trial has actually
  // been recorded, so an abandoned setup screen never touches the roster.
  rememberParticipant(r.participant);
}

/* ---------------- participant roster ---------------- */

/*
 * A returning participant is chosen from a list rather than retyped, and the
 * app refuses to quietly overwrite anything: if the details on screen no
 * longer match the saved record, or the trial number already has data under
 * it, the conflict is shown field by field and Start is blocked until someone
 * decides which version is right.
 */

const P = window.BSRTParticipants;
const PROFILE_KEY = 'bsrt.desktop.participants.v1';

let profiles = {};
let birthEls = null;
let identityDiff = [];
let diffSig = null;          // resets the confirmation when the conflict changes
let dupSig = null;
let duplicateTrial = null;
let fillingForm = false;     // stops fillFromProfile re-entering the checker
let autoFilledFrom = null;   // roster key whose details this app filled in

/* Never write into the field someone is currently typing in: recall is
 * triggered by blurring the ID box, which happens exactly as the next field
 * gains focus, and overwriting it there would leave text half typed over. */
function assignIfIdle(id, value) {
  const el = $(id);
  if (document.activeElement === el) return;
  el.value = value;
}

function birthSelects() {
  if (!birthEls) {
    birthEls = { day: $('pBirthDay'), month: $('pBirthMonth'), year: $('pBirthYear') };
  }
  return birthEls;
}

function renderBirth() {
  P.renderDateSelects(birthSelects(), {
    lang: L.getLanguage(),
    placeholders: {
      day: L.t('participant.day'),
      month: L.t('participant.month'),
      year: L.t('participant.year')
    }
  });
}

/* Year and month first, then re-render so the day list matches that month's
 * length, then the day — otherwise 29 February is not on offer yet. */
function setBirth(iso) {
  const els = birthSelects();
  P.writeDate(els, iso);
  renderBirth();
  P.setDay(els, iso);
}

function readIdentity() {
  return {
    name: txtVal('pName'),
    birthDate: P.readDate(birthSelects()) || '',
    address: txtVal('pAddress'),
    education: txtVal('pEducation')
  };
}

function identityIsEmpty() {
  const e = readIdentity();
  return !e.name && !e.birthDate && !e.address && !e.education;
}

function refreshProfiles() {
  profiles = P.loadProfiles(localStorage, PROFILE_KEY);
  // First run after upgrading: rebuild the roster from trials already stored,
  // so an existing installation does not start with an empty list.
  if (!P.listProfiles(profiles).length) {
    P.seedFromSessions(profiles, loadSessions());
    if (P.listProfiles(profiles).length) P.saveProfiles(localStorage, PROFILE_KEY, profiles);
  }
  renderRoster();
}

function profileLabel(p) {
  const bits = [p.participantId];
  if (p.name) bits.push(p.name);
  bits.push(p.birthDate ? L.t('participant.bornOn', { date: p.birthDate })
                        : L.t('participant.noBirth'));
  return bits.join(' · ');
}

function renderRoster() {
  const sel = $('returning');
  const keep = sel.value;
  sel.innerHTML = '';
  const blank = document.createElement('option');
  blank.value = '';
  blank.setAttribute('data-i18n', 'participant.newParticipant');
  blank.textContent = L.t('participant.newParticipant');
  sel.appendChild(blank);
  P.listProfiles(profiles).forEach((p) => {
    const o = document.createElement('option');
    o.value = P.keyOf(p.participantId);
    o.textContent = profileLabel(p);
    sel.appendChild(o);
  });
  sel.value = keep;
  if (sel.value !== keep) sel.value = '';
}

function fillFromProfile(p) {
  fillingForm = true;
  $('participantId').value = p.participantId || '';
  $('pName').value = p.name || '';
  $('pAddress').value = p.address || '';
  $('pEducation').value = p.education || '';
  setBirth(p.birthDate || '');
  // Only offered, never imposed: a session label already typed is left alone.
  if (p.lastSessionLabel && !txtVal('sessionLabel')) assignIfIdle('sessionLabel', p.lastSessionLabel);
  $('returning').value = P.keyOf(p.participantId);
  autoFilledFrom = P.keyOf(p.participantId);
  fillingForm = false;
  suggestTrialNumber();
  checkParticipant();
}

/* The next number free for this participant and session, so a second visit
 * does not silently land on top of the first. */
function suggestTrialNumber() {
  const id = txtVal('participantId');
  if (!id) return;
  assignIfIdle('trialNumber',
    P.nextTrialNumber(loadSessions(), id, txtVal('sessionLabel') || 'NA'));
}

function diffLine(d) {
  const li = document.createElement('li');
  const name = document.createElement('span');
  name.className = 'fieldname';
  name.textContent = L.t('field.' + d.field) + ': ';
  const was = document.createElement('span');
  was.className = 'was';
  was.textContent = d.from || L.t('field.blank');
  const now = document.createElement('span');
  now.className = 'now';
  now.textContent = d.to || L.t('field.blank');
  li.appendChild(name);
  li.appendChild(was);
  li.appendChild(document.createTextNode(' → '));
  li.appendChild(now);
  return li;
}

function checkParticipant() {
  if (fillingForm) return;

  const id = txtVal('participantId');
  const stored = id ? (profiles[P.keyOf(id)] || null) : null;

  identityDiff = stored ? P.diffIdentity(stored, readIdentity()) : [];
  const sig = identityDiff
    .map((d) => d.field + '\u0001' + d.from + '\u0001' + d.to)
    .join('\u0002');
  // A different conflict is a different decision: never carry a tick across.
  if (sig !== diffSig) { diffSig = sig; $('confirmIdentity').checked = false; }

  const box = $('identityConflict');
  box.hidden = identityDiff.length === 0;
  if (identityDiff.length) {
    const list = $('conflictList');
    list.innerHTML = '';
    identityDiff.forEach((d) => list.appendChild(diffLine(d)));
  }

  duplicateTrial = P.findDuplicate(loadSessions(), id || 'NA',
                                   txtVal('sessionLabel') || 'NA', intVal('trialNumber', 1));
  const dsig = duplicateTrial ? duplicateTrial.participant.runId : '';
  if (dsig !== dupSig) { dupSig = dsig; $('confirmTrial').checked = false; }

  $('trialConflict').hidden = !duplicateTrial;
  if (duplicateTrial) {
    const q = duplicateTrial.participant;
    $('trialConflictDetail').textContent = L.t('conflict.dupDetail', {
      n: q.trialNumber, session: q.sessionLabel, id: q.participantId,
      date: q.date + ' ' + q.time
    });
  }

  if (!stored) $('returning').value = '';
}

function rememberParticipant(m) {
  profiles = P.loadProfiles(localStorage, PROFILE_KEY);
  P.upsertProfile(profiles, m, new Date().toISOString());
  P.saveProfiles(localStorage, PROFILE_KEY, profiles);
  renderRoster();
}

/* True when the form still holds exactly what was auto-filled from a record —
 * i.e. nobody has typed over it. Only such untouched values are ever cleared
 * again; anything a person entered by hand is left alone. */
function identityMatches(p) {
  const e = readIdentity();
  return e.name === (p.name || '') && e.birthDate === (p.birthDate || '') &&
         e.address === (p.address || '') && e.education === (p.education || '');
}

function clearIdentity() {
  $('pName').value = '';
  $('pAddress').value = '';
  $('pEducation').value = '';
  setBirth('');
}

/* Recall runs as the ID is typed rather than when the field is left. Doing it
 * on blur meant writing into the next field at the very moment it received
 * focus, which left whatever was typed there appended to the recalled value. */
function recallIfMatch() {
  const key = P.keyOf(txtVal('participantId'));
  const p = key ? profiles[key] : null;
  if (p) {
    const prev = autoFilledFrom ? profiles[autoFilledFrom] : null;
    if (identityIsEmpty() || (prev && identityMatches(prev))) fillFromProfile(p);
    return;
  }
  // Typed on past a match: take back what was auto-filled, so a new
  // participant cannot inherit someone else's name and birth date.
  const last = autoFilledFrom ? profiles[autoFilledFrom] : null;
  if (last && identityMatches(last)) clearIdentity();
  autoFilledFrom = null;
}

function initParticipants() {
  renderBirth();
  refreshProfiles();

  $('returning').addEventListener('change', () => {
    const p = profiles[$('returning').value];
    if (p) fillFromProfile(p);
  });

  $('participantId').addEventListener('input', () => {
    recallIfMatch();
    suggestTrialNumber();
    checkParticipant();
  });

  // Any hand edit ends the link to the auto-filled record, so nothing this
  // app wrote is ever cleared out from under a person's own typing.
  ['pName', 'pAddress', 'pEducation'].forEach((id) => {
    $(id).addEventListener('input', () => {
      autoFilledFrom = null;
      checkParticipant();
    });
  });
  ['pBirthDay', 'pBirthMonth', 'pBirthYear'].forEach((id) => {
    $(id).addEventListener('change', () => {
      autoFilledFrom = null;
      renderBirth();          // month or year changed: the day list may shrink
      checkParticipant();
    });
  });

  $('sessionLabel').addEventListener('input', () => {
    suggestTrialNumber();
    checkParticipant();
  });
  $('trialNumber').addEventListener('input', checkParticipant);

  $('btnRestoreProfile').addEventListener('click', () => {
    const p = profiles[P.keyOf(txtVal('participantId'))];
    if (p) fillFromProfile(p);
  });
  $('btnNextTrial').addEventListener('click', () => {
    suggestTrialNumber();
    checkParticipant();
  });

  checkParticipant();
}

/* ---------------- startup ---------------- */

async function init() {
  $('sessionCount').textContent = loadSessions().length;
  initParticipants();
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

function parseIsiSet() {
  return $('isiSetMs').value.split(/[,\s]+/)
    .map((v) => Math.round(parseFloat(v) * 1000))
    .filter((v) => isFinite(v) && v > 0);
}

function validateSetup() {
  // Re-run the identity and duplicate checks here rather than trusting the
  // last keystroke to have fired an event.
  checkParticipant();
  if (P.readDate(birthSelects()) === null) return L.t('err.birthIncomplete');
  if (identityDiff.length && !$('confirmIdentity').checked) return L.t('err.identityConflict');
  if (duplicateTrial && !$('confirmTrial').checked) return L.t('err.trialDuplicate');

  const mode = $('mode').value;
  const stim = intVal('stimMs', 1000);
  let shortestIsi;

  if (mode === 'pvt') {
    const set = parseIsiSet();
    if (set.length < 2) return L.t('err.isiSet');
    if (intVal('blockMs', 30) < 1) return L.t('err.blockMin');
    shortestIsi = Math.min.apply(null, set);
  } else {
    shortestIsi = intVal('isiMs', 3000);
    if (shortestIsi < 200) return L.t('err.isiMin');
  }

  if (stim >= shortestIsi) return L.t('err.stimShort', { ms: shortestIsi });
  if (intVal('maxMinutes', 40) < 1) return L.t('err.maxMin');
  if (intVal('lapseMs', 500) >= stim) return L.t('err.lapseBelow');
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

/* ---------------- alarm ---------------- */

/*
 * A continuous pulsing alert, not a one-shot beep: it keeps sounding until the
 * experimenter presses Stop, so a sleep-onset event cannot be missed by someone
 * who stepped away. Synthesised rather than played from a file, so the app stays
 * dependency-free and works offline.
 *
 * A square carrier gated by a square LFO gives the familiar on/off alarm pulse.
 */
function initAudio() {
  if (audioCtx) return;
  try {
    const Ctx = window.AudioContext || window.webkitAudioContext;
    if (Ctx) audioCtx = new Ctx();
  } catch (e) { audioCtx = null; }
}

function startAlarm() {
  if (!cfg || !cfg.alarm || !audioCtx || alarmNodes) return;
  try {
    if (audioCtx.state === 'suspended') audioCtx.resume();
    const carrier = audioCtx.createOscillator();
    const gate = audioCtx.createGain();
    const lfo = audioCtx.createOscillator();
    const lfoDepth = audioCtx.createGain();

    carrier.type = 'square';
    carrier.frequency.value = 880;
    lfo.type = 'square';
    lfo.frequency.value = 4;          // four pulses a second
    lfoDepth.gain.value = 0.11;
    gate.gain.value = 0.11;           // swings between silence and ~0.22

    lfo.connect(lfoDepth);
    lfoDepth.connect(gate.gain);
    carrier.connect(gate);
    gate.connect(audioCtx.destination);
    carrier.start();
    lfo.start();

    alarmNodes = { carrier: carrier, lfo: lfo, gate: gate };
    $('alarmBar').hidden = false;
  } catch (e) { alarmNodes = null; }
}

function stopAlarm() {
  if (alarmNodes) {
    try {
      alarmNodes.carrier.stop();
      alarmNodes.lfo.stop();
      alarmNodes.gate.disconnect();
    } catch (e) { /* already stopped */ }
    alarmNodes = null;
  }
  const bar = $('alarmBar');
  if (bar) bar.hidden = true;
}

/* ---------------- trial ---------------- */

function beginTrial() {
  const isiMs = intVal('isiMs', 3000);
  const stimMs = intVal('stimMs', 1000);
  const corr = $('correction').value;
  const mode = $('mode').value;
  const now = new Date();
  const seedTyped = parseInt($('seed').value, 10);
  const seed = isFinite(seedTyped) ? seedTyped : Math.floor(Math.random() * 2147483647);

  framesStimOn = Math.max(1, Math.round(stimMs / calibration.frameIntervalMs));

  cfg = {
    mode,
    isiMs,
    isiSetMs: mode === 'pvt' ? parseIsiSet() : null,
    blockMs: intVal('blockMs', 30) * 1000,
    seed,
    stimMs,
    /* The hit window is the NOMINAL stimulus duration, not the frame-quantised
     * one, so classification is identical across machines with different
     * refresh rates. Only the light's actual duration is quantised. */
    hitWindowMs: stimMs,
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
    language: L.getLanguage(),
    presentationOffsetFrames: numVal('presOffsetFrames', 1),
    hardwareOffsetMs: numVal('hwOffsetMs', 0),
    photodiode: $('photodiode').checked,
    framesStimOn,
    achievedStimMs: framesStimOn * calibration.frameIntervalMs
  };

  cfg.schedule = S.buildSchedule({
    mode: cfg.mode,
    isiMs: cfg.isiMs,
    isiSetMs: cfg.isiSetMs,
    blockMs: cfg.blockMs,
    maxMs: cfg.maxMs,
    seed: cfg.seed
  });

  /*
   * Stimuli are fired on the first frame at or after their intended time,
   * compared against the REAL frame clock rather than converted into frame
   * indices with the calibrated interval. Converting would multiply any error
   * in that estimate by the elapsed time: a 16.70 ms estimate against a true
   * 16.667 ms drifts 0.2%, which is ~5 s across a 40-minute trial. Comparing
   * against measured time keeps presentation on a frame boundary while bounding
   * the error at half a frame, with nothing accumulating.
   */
  achievedIsis = cfg.schedule.epochIsi.slice();
  cfg.achievedIsiMs = mode === 'pvt' ? null : (achievedIsis[0] == null ? null : achievedIsis[0]);
  cfg.leadInMs = cfg.schedule.leadInMs;

  meta = {
    runId: 'bsrt-' + Date.now(),
    date: now.toISOString().slice(0, 10),
    time: now.toTimeString().slice(0, 8),
    participantId: txtVal('participantId') || 'NA',
    name: txtVal('pName'),
    address: txtVal('pAddress'),
    birthDate: P.readDate(birthSelects()) || '',
    education: txtVal('pEducation'),
    sessionLabel: txtVal('sessionLabel') || 'NA',
    trialNumber: intVal('trialNumber', 1)
  };

  epochs = [];
  presses = [];
  kssBefore = null;
  kssAfter = null;
  currentEpoch = null;
  consecutiveMisses = 0;
  missRunStartIndex = -1;
  blurCount = 0;
  droppedFrames = 0;
  frameIndex = -1;
  lastFrameTs = null;

  nextEpochIdx = 0;
  currentEpochStartFrame = 0;

  $('patch').hidden = !cfg.photodiode;
  initAudio();
  stopAlarm();
  window.bsrt.preventDisplaySleep(true);
  if ($('useFullscreen').checked) window.bsrt.setFullscreen(true);

  if (cfg.kssWhen === 'before' || cfg.kssWhen === 'both') showKss('before');
  else countdown();
}

/* ---------------- Karolinska Sleepiness Scale ---------------- */

function showKss(stage) {
  kssStage = stage;
  $('kssTitle').textContent = L.t(stage === 'before' ? 'kss.beforeTitle' : 'kss.afterTitle');
  $('kssQuestion').textContent = L.t('kss.question');
  $('kssInstruction').textContent = L.t('kss.instruction');

  const box = $('kssOptions');
  box.innerHTML = '';
  L.kssAnchors().forEach((label, i) => {
    const b = document.createElement('button');
    b.className = 'kss-option';
    b.type = 'button';
    const n = document.createElement('span');
    n.className = 'kss-num';
    n.textContent = i + 1;
    const tx = document.createElement('span');
    tx.textContent = label;
    b.appendChild(n);
    b.appendChild(tx);
    b.addEventListener('click', () => answerKss(i + 1));
    box.appendChild(b);
  });
  show('screen-kss');
}

function answerKss(value) {
  if (kssStage === 'before') { kssBefore = value; countdown(); }
  else { kssAfter = value; finishResult(); }
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

  const elapsed = ts - trialStartTs;
  const halfFrame = calibration.frameIntervalMs / 2;

  if (elapsed >= cfg.schedule.plannedDurationMs - halfFrame) {
    finalizeEpoch();
    if (!running) return;
    endTask('max_duration');
    return;
  }

  if (nextEpochIdx < cfg.schedule.nStimuli &&
      elapsed >= cfg.schedule.onsets[nextEpochIdx] - halfFrame) {
    finalizeEpoch();
    if (!running) return;
    startEpoch(nextEpochIdx, ts);
    currentEpochStartFrame = frameIndex;
    nextEpochIdx += 1;
    return;
  }

  if (currentEpoch && cfg.mode === 'pvt' && !currentEpoch.frozen && !$('clock').hidden) {
    // Display only — the counter never feeds the reaction time, which is
    // measured from stimulus presentation to the keypress.
    let shown = ts - currentEpoch.frameTs;
    if (shown < 0) shown = 0;
    if (shown > cfg.hitWindowMs) shown = cfg.hitWindowMs;
    $('clock').textContent = Math.round(shown);
  }

  if (currentEpoch) {
    if (frameIndex - currentEpochStartFrame === framesStimOn) {
      ledOff();
    } else if (currentEpoch.extinguished) {
      // Extinguish here rather than in the event handler so the offset lands on
      // a frame boundary like every other transition — keeps a photodiode trace
      // interpretable.
      ledOff();
      currentEpoch.extinguished = false;
    }
  }
}

function startEpoch(idx, frameTs) {
  ledOn();
  currentEpoch = {
    index: idx,
    frameTs,
    presentationMs: frameTs + cfg.presentationOffsetFrames * calibration.frameIntervalMs,
    onsetMs: frameTs - trialStartTs,
    block: cfg.schedule.blocks[idx],
    // Taken from the schedule, not re-derived from the measured onset: frame
    // quantisation alone can move an onset several ms off its intended time.
    minute: cfg.schedule.minutes[idx],
    // epochIsiMs is this stimulus's response window; isiBeforeMs is the wait
    // that preceded it, defined for the first stimulus too.
    epochIsiMs: cfg.schedule.epochIsi[idx],
    isiBeforeMs: cfg.schedule.isiBefore[idx],
    achievedIsiMs: cfg.schedule.epochIsi[idx],
    responded: false,
    rtRawMs: null,
    extra: 0,
    extinguished: false,
    frozen: false
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
  if (cfg.missCriterion > 0 && consecutiveMisses >= cfg.missCriterion) endTask('sleep_onset');
}

/*
 * Every response is recorded with its raw RT, however late. Hit/miss
 * classification happens at scoring time, not here.
 */
function handleResponse(evt) {
  if (!running) return;
  const t = T.eventTime(evt, useEventStamp);

  // Every press is logged, including presses outside any epoch, because the
  // integrity check needs the full pattern — not just the scored responses.
  presses.push({
    tMs: t - (trialStartTs === null ? t : trialStartTs),
    epochIndex: currentEpoch ? currentEpoch.index : null
  });

  if (!currentEpoch) return;
  if (currentEpoch.responded) { currentEpoch.extra += 1; return; }

  currentEpoch.responded = true;
  currentEpoch.rtRawMs = t - currentEpoch.presentationMs;

  if (cfg.mode === 'pvt') {
    // PVT convention: the counter stops and shows the achieved reaction time.
    currentEpoch.frozen = true;
    $('clock').textContent = Math.max(0, Math.round(currentEpoch.rtRawMs - cfg.hardwareOffsetMs));
    $('clock').classList.add('frozen');
  } else {
    // Extinguishes the light, but never shortens the epoch.
    currentEpoch.extinguished = true;
  }
}

function ledOn() {
  if (cfg.mode === 'pvt') {
    const el = $('clock');
    el.classList.remove('frozen');
    el.textContent = '0';
    el.hidden = false;
  } else {
    $('led').classList.add('on');
  }
  if (cfg.photodiode) $('patch').classList.add('on');
}

function ledOff() {
  $('led').classList.remove('on');
  const el = $('clock');
  el.hidden = true;
  el.classList.remove('frozen');
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

  if (reason === 'sleep_onset') startAlarm();

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
  /*
   * The measured hardware offset is applied BEFORE scoring, because it is a
   * measurement correction rather than a statistical exclusion: with it, rtMs
   * is the participant's actual reaction time, so hit/miss classification uses
   * the right number. Both raw and offset-applied values are exported.
   */
  const scored = S.score(epochs.map((e) => ({
    index: e.index,
    onsetMs: e.onsetMs,
    rtMs: e.rtRawMs === null ? null : e.rtRawMs - cfg.hardwareOffsetMs,
    block: e.block,
    minute: e.minute,
    epochIsiMs: e.epochIsiMs,
    isiBeforeMs: e.isiBeforeMs
  })), cfg, presses);

  let sleepOnsetMs = null, sleepOnsetCriterionMs = null;
  if (reason === 'sleep_onset' && missRunStartIndex >= 0) {
    const first = epochs[missRunStartIndex];
    const last = epochs[missRunStartIndex + cfg.missCriterion - 1];
    sleepOnsetMs = first ? first.onsetMs : cfg.schedule.onsets[missRunStartIndex];
    // With a variable schedule the criterion is confirmed one INTERVAL after the
    // last missed stimulus, and that interval differs per epoch.
    sleepOnsetCriterionMs = last ? last.onsetMs + last.achievedIsiMs
                                 : sleepOnsetMs + cfg.missCriterion * (achievedIsis[0] || cfg.isiMs);
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
    elapsedMs: epochs.length ? epochs[epochs.length - 1].onsetMs + epochs[epochs.length - 1].achievedIsiMs : 0,
    norms: S.normativeReport({
      trials: scored.trials,
      elapsedMs: epochs.length ? epochs[epochs.length - 1].onsetMs + epochs[epochs.length - 1].achievedIsiMs : 0,
      minuteBuckets: scored.perMinute.length,
      mode: cfg.mode,
      isiMs: cfg.isiMs,
      hour: parseInt(meta.time.slice(0, 2), 10),
      lapseMs: cfg.lapseMs,
      falseStartMs: cfg.falseStartMs,
      missCriterion: cfg.missCriterion
    }, window.BSRTNorms),
    extraResponses: epochs.reduce((n, e) => n + e.extra, 0),
    kssBefore: null,
    kssAfter: null,
    language: cfg.language,
    extras: epochs.map((e) => e.extra),
    rawRts: epochs.map((e) => e.rtRawMs),
    droppedFramesDuringTrial: droppedFrames,
    framesDuringTrial: frameIndex + 1,
    blurCount,
    scored
  };
}

/* ---------------- rendering ---------------- */

/* A dropped frame delays a stimulus onset by up to one frame interval. The
 * onset is measured rather than assumed, so the delay lands in the recorded
 * ISI and never in the reaction time — which is why the rate, not the raw
 * count, is what matters. Anything under about 1% is ordinary background
 * scheduling noise. */
function formatDropped(r) {
  const dropped = r.droppedFramesDuringTrial;
  const shown = r.framesDuringTrial;
  if (!shown) return String(dropped);
  const expected = shown + dropped;
  const pct = (dropped / expected) * 100;
  const verdict = pct < 1 ? 'negligible' : pct < 5 ? 'noticeable' : 'high — close other applications';
  return dropped.toLocaleString() + ' of ' + expected.toLocaleString() +
         ' frames (' + pct.toFixed(2) + '% — ' + verdict + ')';
}

function renderResult(r) {
  const p = r.participant, t = r.scored.totals, e = r.scored.errorProfiles;

  $('resultWho').textContent =
    (p.name || p.participantId) + ' · ' + p.sessionLabel + ' · trial ' + p.trialNumber +
    ' · ' + p.date + ' ' + p.time;

  if (r.sleptBeforeMax) {
    $('mLatency').textContent = fmtClock(r.sleepOnsetMs);
    $('mLatencyFoot').textContent =
      L.t('results.latencyFoot', { n: r.config.missCriterion, t: fmtClock(r.sleepOnsetCriterionMs) });
    $('mOutcome').textContent = L.t('results.outcomeSlept', { n: r.config.maxMinutes });
  } else {
    $('mLatency').textContent = '> ' + fmtClock(r.config.maxMs);
    $('mLatencyFoot').textContent = L.t('results.censored');
    $('mOutcome').textContent = r.endReason === 'aborted'
      ? L.t('results.outcomeAborted')
      : L.t('results.outcomeFull', { n: r.config.maxMinutes });
  }

  $('oTrials').textContent = t.trials;
  $('oHitRatio').textContent = t.hitRatio == null ? '—' : (t.hitRatio * 100).toFixed(1) + '%';
  $('oHits').textContent = t.hits;
  $('oMisses').textContent = t.lateResponses
    ? L.t('results.lateSuffix', { n: t.misses, late: t.lateResponses }) : t.misses;
  $('oLapses').textContent = L.t('results.lapseSuffix', { n: t.lapses, ms: r.config.lapseMs });
  $('oFalseStarts').textContent = L.t('results.fsSuffix', { n: t.falseStarts, ms: r.config.falseStartMs });
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
  $('rtIpr').textContent = n1(t.iprRt);
  $('rtcIpr').textContent = n1(t.corrIprRt);

  $('rsAvg').textContent = n3(t.avgRs);
  $('rsMedian').textContent = n3(t.medianRs);
  $('rsSd').textContent = n3(t.sdRs);
  $('rsFast').textContent = n3(t.fastest10Rs);
  $('rsSlow').textContent = n3(t.slowest10Rs);
  $('rscAvg').textContent = n3(t.corrAvgRs);
  $('rscMedian').textContent = n3(t.corrMedianRs);
  $('rscFast').textContent = n3(t.corrFastest10Rs);
  $('rscSlow').textContent = n3(t.corrSlowest10Rs);
  $('rsIpr').textContent = n3(t.iprRs);
  $('rscIpr').textContent = n3(t.corrIprRs);

  $('corrNote').textContent = r.config.correction === 'none'
    ? L.t('corr.none')
    : L.t('corr.applied', { desc: describeCorrection(r.config),
                            a: t.nFalseStartsRemoved, b: t.nOutliersRemoved, n: t.n });

  renderNorms(r);
  renderIntegrity(r);
  renderKss(r);
  renderPerMinute(r);

  $('tRefresh').textContent = r.calibration.refreshHz.toFixed(2) + ' Hz';
  $('tInterval').textContent = n2(r.calibration.frameIntervalMs) + ' ms';
  $('tAchievedIsi').textContent = r.config.mode === 'pvt'
    ? 'variable — ' + r.config.isiSetMs.map((v) => v / 1000).join('/') + ' s, ' +
      r.config.schedule.method.replace(/_/g, ' ')
    : n2(r.config.achievedIsiMs) + ' ms (requested ' + r.config.isiMs + ')';
  $('tQuant').textContent = '± ' + n2(r.calibration.frameIntervalMs / 2) + ' ms';
  $('tDropped').textContent = formatDropped(r);
  $('tInputSrc').textContent = r.inputTimeSource === 'os_event_stamp' ? 'OS event stamp' : 'handler time (fallback)';
  $('tInputDelay').textContent = r.inputProbe && r.inputProbe.usable
    ? n2(r.inputProbe.medianDelayMs) + ' ms (removed)' : '—';
  $('tHwOffset').textContent = r.config.hardwareOffsetMs
    ? r.config.hardwareOffsetMs + ' ms (applied before scoring)' : 'not measured';
}

function describeCorrection(c) {
  if (c.correction === 'falseStarts') return L.t('corr.descFs', { ms: c.falseStartMs });
  if (c.correction === 'outliers') return L.t('corr.descOut', { sd: c.sdMultiplier });
  if (c.correction === 'both') return L.t('corr.descBoth', { ms: c.falseStartMs, sd: c.sdMultiplier });
  return L.t('corrections.none');
}

/* ---------------- normative comparison ---------------- */

/*
 * Labels reuse the keys the rest of the results screen already uses, so the
 * panel speaks the interface language without another 84 strings to translate.
 * The reference table's own English labels stay in norms.js as the data's
 * documentation.
 */
const NORM_LABEL = {
  trials: 'results.totalTrials', hitRatio: 'results.hitRatio', hits: 'results.hits',
  misses: 'results.misses', lapses: 'results.lapses', falseStarts: 'results.falseStarts',
  ep12: 'results.ep12', ep36: 'results.ep36', ep7: 'results.ep7',
  rtMean: 'results.average', rtMedian: 'results.median', rtSd: 'results.sd',
  rtFast10: 'results.fastest', rtSlow10: 'results.slowest', rtIpr: 'results.ipr',
  rsMean: 'results.average', rsMedian: 'results.median', rsSd: 'results.sd',
  rsFast10: 'results.fastest', rsSlow10: 'results.slowest', rsIpr: 'results.ipr'
};
const NORM_SECTION = {
  'Test overview': 'results.overview',
  'Reaction time': 'norms.sectRt',
  'Reaction speed': 'norms.sectRs'
};

function normFmt(key, v) {
  if (v == null) return '—';
  if (key.indexOf('rs') === 0) return v.toFixed(3);
  if (key === 'hitRatio') return v.toFixed(1) + '%';
  if (key.indexOf('rt') === 0) return v.toFixed(1);
  return String(Math.round(v * 100) / 100);
}

function normReasonText(nm) {
  if (!nm) return L.t('norms.rNoNorms');
  switch (nm.reason) {
    case 'not_bsrt': return L.t('norms.rNotBsrt');
    case 'isi_mismatch': return L.t('norms.rIsi', { isi: nm.isiMs });
    case 'no_hour': return L.t('norms.rNoHour');
    case 'too_short': return L.t('norms.rTooShort');
    default: return L.t('norms.rNoNorms');
  }
}

function renderNorms(r) {
  var nm = r.norms;
  var body = $('normsBody');
  var un = $('normsUnavailable');

  if (!nm || !nm.available) {
    body.hidden = true;
    un.hidden = false;
    un.textContent = L.t('norms.unavailable') + ' ' + normReasonText(nm);
    $('normsContext').textContent = '';
    return;
  }

  un.hidden = true;
  body.hidden = false;

  var hh = (nm.hour < 10 ? '0' : '') + nm.hour;
  var summary = (nm.nOrange + nm.nRed) === 0
    ? L.t('norms.summaryClean')
    : L.t('norms.summaryFlagged', { orange: nm.nOrange, red: nm.nRed });
  $('normsContext').textContent =
    L.t('norms.context', { hour: hh, min: nm.windowMinutes, n: nm.n }) + ' ' + summary;

  var tbody = $('normsTable').getElementsByTagName('tbody')[0];
  tbody.innerHTML = '';
  var lastSection = null;

  nm.rows.forEach(function (row) {
    var c = row.comparison;
    var tr = document.createElement('tr');
    if (row.section !== lastSection) { tr.className = 'section-start'; lastSection = row.section; }

    var th = document.createElement('th');
    if (tr.className === 'section-start') {
      var sect = document.createElement('span');
      sect.className = 'sect';
      sect.textContent = L.t(NORM_SECTION[row.section] || row.section);
      th.appendChild(sect);
    }
    th.appendChild(document.createTextNode(L.t(NORM_LABEL[row.key] || row.label)));
    tr.appendChild(th);

    var tdVal = document.createElement('td');
    tdVal.textContent = normFmt(row.key, c ? c.value : null);
    // The value carries the band as well as the z cell, so the thing being
    // judged is the thing that looks red.
    if (c) tdVal.className = 'band-' + c.band;
    tr.appendChild(tdVal);

    var tdRef = document.createElement('td');
    tdRef.textContent = c ? normFmt(row.key, c.mean) + ' ± ' + normFmt(row.key, c.sd) : '—';
    tr.appendChild(tdRef);

    var tdZ = document.createElement('td');
    if (!c) {
      tdZ.className = 'band-none';
      tdZ.textContent = '—';
    } else if (c.dir === 0) {
      tdZ.className = 'band-neutral';
      tdZ.textContent = L.t('norms.notJudged');
    } else if (c.degenerate) {
      tdZ.className = 'band-' + c.band;
      // No spread in the reference, so no z exists — say that rather than
      // dividing by zero and printing an infinity.
      tdZ.textContent = c.band === 'red' ? L.t('norms.noSpread') : '—';
    } else {
      tdZ.className = 'band-' + c.band;
      // Signed so that positive always reads "worse", whichever direction is
      // bad for this variable, and labelled so the sign is never ambiguous.
      var word = c.z >= 0 ? L.t('norms.worse') : L.t('norms.better');
      tdZ.textContent = Math.abs(c.z).toFixed(2) + ' SD ' + word;
    }
    tr.appendChild(tdZ);
    tbody.appendChild(tr);
  });

  var items = [L.t('norms.dPreliminary', { sessions: nm.sessions, participants: nm.participants })];
  if (nm.belowProtocol) {
    items.push(L.t('norms.dShort', { protocol: nm.protocolMinutes, min: nm.windowMinutes }));
  }
  if (nm.truncated) {
    items.push(L.t('norms.dTruncated', { test: nm.testMinutes, protocol: nm.protocolMinutes }));
  }
  // A reference SD of well under one count turns a two-trial difference into a
  // double-digit z. Say so when it actually happens rather than always.
  var extreme = nm.rows.some(function (x) {
    return x.comparison && x.comparison.z != null && Math.abs(x.comparison.z) > 5;
  });
  if (extreme) items.push(L.t('norms.dExtreme'));
  items.push(L.t('norms.dConvention'));
  items.push(L.t('norms.dHour', { hour: hh }));

  var html = '<strong>' + L.t('norms.disclaimerHeading') + '</strong><ul>';
  for (var i = 0; i < items.length; i++) html += '<li>' + items[i] + '</li>';
  $('normsDisclaimer').innerHTML = html + '</ul>';
}

function renderIntegrity(r) {
  const g = r.scored.integrity;
  $('inPresses').textContent = g.totalPresses;
  $('inExtra').textContent = g.extraRate
    ? L.t('results.extraVal', { n: g.extraPresses, pct: (g.extraRate * 100).toFixed(0) })
    : g.extraPresses;
  $('inBurst').textContent = L.t('results.burstVal', { n: g.burstMax, ms: g.thresholds.burstWindowMs });
  $('inRapid').textContent = L.t('results.rapidVal', { n: g.rapidPairs, ms: g.thresholds.rapidGapMs });

  const box = $('integrityVerdict');
  if (g.suspected) {
    box.className = 'note warnbox';
    box.textContent = L.t('results.integrityFlag') + ' ' + g.reasons.join('; ') + '.';
  } else {
    box.className = 'note tight';
    box.textContent = L.t('results.integrityOk');
  }
}

function renderKss(r) {
  const card = $('kssResult');
  if (r.kssBefore == null && r.kssAfter == null) { card.hidden = true; return; }
  card.hidden = false;
  const anchors = L.kssAnchors(r.language);
  $('kssBeforeVal').textContent = r.kssBefore == null ? '—' : r.kssBefore + ' — ' + anchors[r.kssBefore - 1];
  $('kssAfterVal').textContent = r.kssAfter == null ? '—' : r.kssAfter + ' — ' + anchors[r.kssAfter - 1];
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
    ? L.t('results.noTrend')
    : L.t('results.trend', { v: n3(r.scored.dynamicsRs.slope) });
}

/* ---------------- exports ---------------- */

const PARTICIPANT_COLS = ['run_id', 'date', 'time', 'language', 'participant_id', 'name', 'address',
                          'birth_date', 'educational_level', 'session_label', 'trial_number'];

function participantVals(r) {
  const p = r.participant;
  return [r.runId, p.date, p.time, r.language || 'en', p.participantId, p.name, p.address,
          p.birthDate, p.education, p.sessionLabel, p.trialNumber];
}

const RAW_HEADER = PARTICIPANT_COLS.concat([
  'epoch_index', 'block', 'minute', 'onset_ms', 'epoch_isi_ms', 'isi_before_ms', 'responded',
  'rt_raw_ms', 'rt_ms', 'rs_per_sec',
  'outcome', 'lapse', 'late_response', 'false_start', 'extra_responses'
]);

function rawRows(r) {
  const pv = participantVals(r);
  const extras = r.extras || [];
  const raws = r.rawRts || [];
  return r.scored.trials.map((t, i) => pv.concat([
    // 1-based on the way out, like block and minute: epoch_index was the only
    // column that still started at zero.
    t.index + 1, t.block === null ? '' : t.block + 1, t.minute + 1, round(t.onsetMs, 2),
    t.epochIsiMs, t.isiBeforeMs,
    t.rtMs === null ? 0 : 1,
    round(raws[i] == null ? null : raws[i], 3),   // before the hardware offset
    round(t.rtMs, 3),                              // after it — what was scored
    round(t.rsPerSec, 5),
    t.outcome, t.lapse, t.lateResponse, t.falseStart,
    extras[i] == null ? '' : extras[i]
  ]));
}

const PM_HEADER = PARTICIPANT_COLS.concat([
  'minute', 'trials', 'hits', 'misses', 'lapses', 'late_responses', 'hit_ratio',
  'avg_rt', 'median_rt', 'stdev_rt', 'fastest10_rt', 'slowest10_rt',
  'corr_avg_rt', 'corr_median_rt', 'corr_stdev_rt', 'corr_fastest10_rt', 'corr_slowest10_rt',
  'avg_rs', 'median_rs', 'stdev_rs', 'fastest10_rs', 'slowest10_rs',
  'corr_avg_rs', 'corr_median_rs', 'corr_stdev_rs', 'corr_fastest10_rs', 'corr_slowest10_rs',
  'n_rt', 'n_rt_corrected', 'false_starts_removed', 'outliers_removed',
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
    m.n, m.nCorrected, m.nFalseStartsRemoved, m.nOutliersRemoved,
    round(vRs.velocity[i], 5), round(vRs.acceleration[i], 5),
    round(vRt.velocity[i], 3), round(vRt.acceleration[i], 3)
  ]));
}

/*
 * Normative comparison, long format: one row per variable so files rbind
 * cleanly, the same shape as the other exports. Values here follow the
 * REFERENCE workbook's definitions, not the app's own scoring — see
 * normativeSummary in scoring.js — so norm_value will not always equal the
 * matching column in the summary export. Both are kept on purpose.
 */
const NORMS_HEADER = PARTICIPANT_COLS.concat([
  'norm_available', 'norm_reason', 'norm_source', 'norm_hour', 'norm_window_min',
  'norm_test_min', 'norm_protocol_min', 'norm_below_protocol', 'norm_truncated',
  'norm_ref_sessions', 'norm_ref_n_hour_bin',
  'section', 'variable', 'unit', 'higher_is',
  'norm_value', 'norm_ref_mean', 'norm_ref_sd', 'z_worse', 'band', 'zero_variance_reference'
]);

function normsRows(r) {
  var pv = participantVals(r);
  var nm = r.norms;
  if (!nm || !nm.available) {
    return [pv.concat([0, nm ? nm.reason : 'no_norms', '', '', '', '', '', '', '', '', '',
                       '', '', '', '', '', '', '', '', '', ''])];
  }
  var ctx = [1, '', nm.source, nm.hour, nm.windowMinutes, nm.testMinutes, nm.protocolMinutes,
             nm.belowProtocol ? 1 : 0, nm.truncated ? 1 : 0, nm.sessions, nm.n];
  return nm.rows.map(function (row) {
    var c = row.comparison;
    return pv.concat(ctx).concat([
      row.section, row.label, row.unit,
      row.dir === 0 ? 'not judged' : (row.dir > 0 ? 'better' : 'worse'),
      c ? round(c.value, 5) : null,
      c ? round(c.mean, 5) : null,
      c ? round(c.sd, 5) : null,
      c && c.z != null ? round(c.z, 4) : null,
      c ? c.band : '',
      c && c.degenerate ? 1 : 0
    ]);
  });
}

const SUMMARY_HEADER = PARTICIPANT_COLS.concat([
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
  'mode', 'isi_requested_ms', 'isi_achieved_ms', 'isi_set_s', 'block_s',
  'schedule_method', 'schedule_seed', 'scheduled_stimuli',
  'stim_requested_ms', 'stim_achieved_ms',
  'hit_window_ms', 'lapse_threshold_ms', 'miss_criterion', 'max_minutes',
  'correction', 'false_start_threshold_ms', 'sd_multiplier',
  'kss_when', 'kss_before', 'kss_after',
  'total_presses', 'extra_presses', 'burst_max', 'rapid_pairs', 'cheating_suspected', 'cheating_reasons',
  'alarm_enabled',
  'refresh_hz_measured', 'refresh_hz_reported', 'frame_interval_ms', 'frame_jitter_mad_ms',
  'onset_quantisation_ms', 'dropped_frames_calibration', 'dropped_frames_trial',
  'frames_trial', 'dropped_rate_trial', 'calibration_grade',
  'input_time_source', 'input_dispatch_median_ms',
  'presentation_offset_frames', 'hardware_offset_ms',
  'platform', 'os_release', 'electron_version', 'chrome_version', 'display_scale_factor',
  'extra_responses', 'page_blur_count'
]);

function summaryRow(r) {
  const t = r.scored.totals, e = r.scored.errorProfiles, c = r.config;
  const g = r.scored.integrity;
  const cal = r.calibration, d = r.display ? r.display.current : {};
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
    c.mode, c.mode === 'pvt' ? '' : c.isiMs, c.mode === 'pvt' ? '' : round(c.achievedIsiMs, 3),
    c.isiSetMs ? c.isiSetMs.map((v) => v / 1000).join(' ') : '',
    c.blockMs / 1000, c.schedule.method, c.schedule.seed, c.schedule.nStimuli,
    c.stimMs, round(c.achievedStimMs, 3),
    c.hitWindowMs, c.lapseMs, c.missCriterion, c.maxMinutes,
    c.correction, c.falseStartMs, c.sdMultiplier,
    c.kssWhen, r.kssBefore, r.kssAfter,
    g.totalPresses, g.extraPresses, g.burstMax, g.rapidPairs,
    g.suspected ? 1 : 0, g.reasons.join('; '),
    c.alarm ? 1 : 0,
    round(cal.refreshHz, 3), d.displayFrequency == null ? null : d.displayFrequency,
    round(cal.frameIntervalMs, 4), round(cal.madIntervalMs, 4),
    round(cal.frameIntervalMs / 2, 3), cal.droppedFrames, r.droppedFramesDuringTrial,
    r.framesDuringTrial,
    r.framesDuringTrial
      ? round(r.droppedFramesDuringTrial / (r.framesDuringTrial + r.droppedFramesDuringTrial), 5)
      : null,
    r.calibrationGrade,
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
  if (!all.length) { alert(L.t('export.none')); return; }
  let rows = [];
  all.forEach((r) => { rows = rows.concat(builder(r)); });
  saveCsv(name, toCsv(header, rows));
}

/* ---------------- wiring ---------------- */

$('btnStart').addEventListener('click', startCalibration);
$('btnBegin').addEventListener('click', beginTrial);
$('btnCalBack').addEventListener('click', () => show('screen-setup'));

$('btnStopAlarm').addEventListener('click', stopAlarm);

$('btnAgain').addEventListener('click', () => {
  stopAlarm();
  // Straight to the next free number for this participant and session, so a
  // second run cannot land on a trial number that already holds data.
  suggestTrialNumber();
  checkParticipant();
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

$('btnExportNorms').addEventListener('click', () => {
  if (lastResult) saveCsv(fileStem(lastResult) + '_norms.csv', toCsv(NORMS_HEADER, normsRows(lastResult)));
});
$('btnExportAllNorms').addEventListener('click', () => bulk(normsRows, NORMS_HEADER, 'bsrt_all_norms.csv'));

$('btnExportAllRaw').addEventListener('click', () => bulk(rawRows, RAW_HEADER, 'bsrt_all_raw.csv'));
$('btnExportAllPm').addEventListener('click', () => bulk(pmRows, PM_HEADER, 'bsrt_all_perminute.csv'));
$('btnExportAllSummary').addEventListener('click', () => {
  const all = loadSessions();
  if (!all.length) { alert(L.t('export.none')); return; }
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
  $('lapseHint').textContent = L.t(short ? 'lapse.short' : 'lapse.regular');
});

function applyLanguage() {
  L.setLanguage($('language').value);
  L.applyTranslations(document);
  $('lapseHint').textContent = L.t(intVal('maxMinutes', 40) < 5 ? 'lapse.short' : 'lapse.regular');
  applyMode();
  // Month names and the roster labels are built in code, so they need
  // rebuilding by hand when the language changes.
  renderBirth();
  renderRoster();
  if (lastResult) renderResult(lastResult);
}

function applyMode() {
  const mode = $('mode').value;
  document.body.setAttribute('data-mode', mode);
  // A PVT runs for a fixed duration; the sleep-onset criterion is an OSLER idea.
  $('criterionOn').checked = mode !== 'pvt';
  $('modeNote').textContent = mode === 'pvt'
    ? 'PVT: intervals vary within each block, drawn so that every block contains one of each. The sleep-onset criterion is off by default — a PVT runs to time.'
    : 'BSRT / OSLER: a fixed interval between stimuli, with sleep onset scored from consecutive misses.';
  $('instructionsText').textContent = L.t(mode === 'pvt' ? 'instructions.pvt' : 'instructions.bsrt');
  $('taskHint').textContent = L.t(mode === 'pvt' ? 'task.hintPvt' : 'task.hintBsrt');
}

$('mode').addEventListener('change', applyMode);
$('language').addEventListener('change', applyLanguage);
applyLanguage();
init();
