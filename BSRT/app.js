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
var T = window.BSRTTiming;
var D = window.BSRTDiagnostics;

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
var alarmNodes = null;

/* ---- device diagnostics ---- */
var env = null;            // passive environment, read once at load
var deviceCal = null;      // display calibration measured during the countdown
var inputProbe = null;     // fed by every real response during the trial
var monitorRaf = null;
var frameCount = 0;
var droppedInTrial = 0;
var lastFrameTs = null;
var frameIntervals = [];   // bounded sample, enough for a jitter estimate
var touchUsed = false;
var wasFullscreen = false;
var rtStampRejected = 0;   // responses whose event stamp failed its sanity check
var deviceCheckStop = null;

/* ---------------- helpers ---------------- */

function $(id) { return document.getElementById(id); }

function show(id) {
  var all = document.querySelectorAll('.screen');
  for (var i = 0; i < all.length; i++) all[i].classList.remove('active');
  $(id).classList.add('active');
  // Lets CSS hide the floating language control on the dark task screens.
  document.body.setAttribute('data-screen', id);
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
  // A participant record is created or updated only once a trial has actually
  // been recorded, so an abandoned setup screen never touches the roster.
  rememberParticipant(r.participant);
}

function refreshCount() { $('sessionCount').textContent = loadSessions().length; }

/* ---------------- participant roster ---------------- */

/*
 * A returning participant is chosen from a list rather than retyped, and the
 * app refuses to quietly overwrite anything: if the details on screen no
 * longer match the saved record, or the trial number already has data under
 * it, the conflict is shown field by field and Start is blocked until someone
 * decides which version is right.
 */

var P = window.BSRTParticipants;
var PROFILE_KEY = 'bsrt.participants.v1';

var profiles = {};
var birthEls = null;
var identityDiff = [];
var diffSig = null;          // resets the confirmation when the conflict changes
var dupSig = null;
var duplicateTrial = null;
var fillingForm = false;     // stops fillFromProfile re-entering the checker
var autoFilledFrom = null;   // roster key whose details this app filled in

/* Never write into the field someone is currently typing in: recall is
 * triggered by blurring the ID box, which happens exactly as the next field
 * gains focus, and overwriting it there would leave text half typed over. */
function assignIfIdle(id, value) {
  var el = $(id);
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
  var els = birthSelects();
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
  var e = readIdentity();
  return !e.name && !e.birthDate && !e.address && !e.education;
}

function refreshProfiles() {
  // Seeds from stored trials on the first run only. An empty roster is a
  // choice once someone has cleared it, not a gap to be refilled.
  profiles = P.loadOrSeed(localStorage, PROFILE_KEY, loadSessions());
  renderRoster();
}

function profileLabel(p) {
  var bits = [p.participantId];
  if (p.name) bits.push(p.name);
  bits.push(p.birthDate ? L.t('participant.bornOn', { date: p.birthDate })
                        : L.t('participant.noBirth'));
  return bits.join(' · ');
}

function renderRoster() {
  var sel = $('returning');
  var keep = sel.value;
  sel.innerHTML = '';
  var blank = document.createElement('option');
  blank.value = '';
  blank.setAttribute('data-i18n', 'participant.newParticipant');
  blank.textContent = L.t('participant.newParticipant');
  sel.appendChild(blank);
  P.listProfiles(profiles).forEach(function (p) {
    var o = document.createElement('option');
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
  var id = txtVal('participantId');
  if (!id) return;
  assignIfIdle('trialNumber',
    P.nextTrialNumber(loadSessions(), id, txtVal('sessionLabel') || 'NA'));
}

function diffLine(d) {
  var li = document.createElement('li');
  var name = document.createElement('span');
  name.className = 'fieldname';
  name.textContent = L.t('field.' + d.field) + ': ';
  var was = document.createElement('span');
  was.className = 'was';
  was.textContent = d.from || L.t('field.blank');
  var now = document.createElement('span');
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

  var id = txtVal('participantId');
  var stored = id ? (profiles[P.keyOf(id)] || null) : null;

  identityDiff = stored ? P.diffIdentity(stored, readIdentity()) : [];
  var sig = identityDiff.map(function (d) {
    return d.field + '\u0001' + d.from + '\u0001' + d.to;
  }).join('\u0002');
  // A different conflict is a different decision: never carry a tick across.
  if (sig !== diffSig) { diffSig = sig; $('confirmIdentity').checked = false; }

  var box = $('identityConflict');
  box.hidden = identityDiff.length === 0;
  if (identityDiff.length) {
    var list = $('conflictList');
    list.innerHTML = '';
    identityDiff.forEach(function (d) { list.appendChild(diffLine(d)); });
  }

  duplicateTrial = P.findDuplicate(loadSessions(), id || 'NA',
                                   txtVal('sessionLabel') || 'NA', intVal('trialNumber', 1));
  var dsig = duplicateTrial ? duplicateTrial.participant.runId : '';
  if (dsig !== dupSig) { dupSig = dsig; $('confirmTrial').checked = false; }

  $('trialConflict').hidden = !duplicateTrial;
  if (duplicateTrial) {
    var q = duplicateTrial.participant;
    $('trialConflictDetail').textContent = L.t('conflict.dupDetail', {
      n: q.trialNumber, session: q.sessionLabel, id: q.participantId,
      date: q.date + ' ' + q.time
    });
  }

  if (!stored) $('returning').value = '';
}

function rememberParticipant(m) {
  refreshProfilesQuietly();
  P.upsertProfile(profiles, m, new Date().toISOString());
  P.saveProfiles(localStorage, PROFILE_KEY, profiles);
  renderRoster();
}

function refreshProfilesQuietly() { profiles = P.loadProfiles(localStorage, PROFILE_KEY); }

/* True when the form still holds exactly what was auto-filled from a record —
 * i.e. nobody has typed over it. Only such untouched values are ever cleared
 * again; anything a person entered by hand is left alone. */
function identityMatches(p) {
  var e = readIdentity();
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
  var key = P.keyOf(txtVal('participantId'));
  var p = key ? profiles[key] : null;
  if (p) {
    var prev = autoFilledFrom ? profiles[autoFilledFrom] : null;
    if (identityIsEmpty() || (prev && identityMatches(prev))) fillFromProfile(p);
    return;
  }
  // Typed on past a match: take back what was auto-filled, so a new
  // participant cannot inherit someone else's name and birth date.
  var last = autoFilledFrom ? profiles[autoFilledFrom] : null;
  if (last && identityMatches(last)) clearIdentity();
  autoFilledFrom = null;
}

/*
 * Removing someone from the list never removes their data. The roster is a
 * convenience for filling the form in; recorded trials live separately and are
 * left exactly as they were.
 */
function forgetParticipant() {
  var key = $('returning').value || P.keyOf(txtVal('participantId'));
  var p = key ? profiles[key] : null;
  if (!p) { alert(L.t('participant.forgetNone')); return; }
  var who = p.participantId + (p.name ? ' (' + p.name + ')' : '');
  if (!confirm(L.t('participant.forgetConfirm', { who: who }))) return;
  P.forgetProfile(localStorage, PROFILE_KEY, profiles, p.participantId);
  if (autoFilledFrom === P.keyOf(p.participantId)) {
    clearIdentity();
    $('participantId').value = '';
    autoFilledFrom = null;
  }
  $('returning').value = '';
  renderRoster();
  checkParticipant();
}

function forgetAllParticipants() {
  var n = P.listProfiles(profiles).length;
  if (!n) { alert(L.t('participant.forgetNone')); return; }
  if (!confirm(L.t('participant.forgetAllConfirm', { n: n }))) return;
  profiles = P.forgetAllProfiles(localStorage, PROFILE_KEY);
  autoFilledFrom = null;
  $('returning').value = '';
  renderRoster();
  checkParticipant();
}

function initParticipants() {
  renderBirth();
  refreshProfiles();

  $('returning').addEventListener('change', function () {
    var p = profiles[$('returning').value];
    if (p) fillFromProfile(p);
  });

  $('participantId').addEventListener('input', function () {
    recallIfMatch();
    suggestTrialNumber();
    checkParticipant();
  });

  // Any hand edit ends the link to the auto-filled record, so nothing this
  // app wrote is ever cleared out from under a person's own typing.
  ['pName', 'pAddress', 'pEducation'].forEach(function (id) {
    $(id).addEventListener('input', function () {
      autoFilledFrom = null;
      checkParticipant();
    });
  });
  ['pBirthDay', 'pBirthMonth', 'pBirthYear'].forEach(function (id) {
    $(id).addEventListener('change', function () {
      autoFilledFrom = null;
      renderBirth();          // month or year changed: the day list may shrink
      checkParticipant();
    });
  });

  $('sessionLabel').addEventListener('input', function () {
    suggestTrialNumber();
    checkParticipant();
  });
  $('trialNumber').addEventListener('input', checkParticipant);

  $('btnRestoreProfile').addEventListener('click', function () {
    var p = profiles[P.keyOf(txtVal('participantId'))];
    if (p) fillFromProfile(p);
  });
  $('btnNextTrial').addEventListener('click', function () {
    suggestTrialNumber();
    checkParticipant();
  });

  $('btnForgetParticipant').addEventListener('click', forgetParticipant);
  $('btnForgetAll').addEventListener('click', forgetAllParticipants);
  checkParticipant();
}

/* ---------------- device diagnostics ---------------- */

/*
 * A browser cannot ask the OS what monitor it is on, so everything here is
 * either read from the page or measured. See diagnostics.js for why network
 * latency is recorded as provenance and never as a timing figure.
 */

function initDiagnostics() {
  env = D.snapshot();
  renderEnvLine();
}

function renderEnvLine() {
  if (!env) return;
  var bits = [D.describeBrowser(env.userAgent)];
  if (env.screenW) bits.push(env.screenW + '×' + env.screenH);
  if (env.devicePixelRatio) bits.push(L.t('diag.scale', { x: env.devicePixelRatio }));
  if (deviceCal && deviceCal.refreshHz) {
    bits.push(deviceCal.refreshHz.toFixed(1) + ' Hz');
  } else {
    bits.push(L.t('diag.hzUnmeasured'));
  }
  if (env.timerResolutionMs != null) {
    bits.push(L.t('diag.timer', { res: env.timerResolutionMs.toFixed(2) }));
  }
  $('envInfo').textContent = bits.join(' · ');

  // Only the things that would actually spoil a measurement get a warning.
  var warn = [];
  if (env.timerResolutionMs != null && env.timerResolutionMs >= 1) {
    warn.push(L.t('diag.warnTimer', { res: env.timerResolutionMs.toFixed(2) }));
  }
  if (env.maxTouchPoints > 0 && /Android|iPhone|iPad/.test(env.userAgent || '')) {
    warn.push(L.t('diag.warnTouch'));
  }
  $('envWarn').hidden = !warn.length;
  $('envWarn').textContent = warn.join(' ');
}

/* Frame monitoring during the trial: the same idea as the desktop build's
 * frame loop, but purely observational — this build still schedules with
 * setTimeout, and the monitor never touches the epoch clock. */
function startFrameMonitor() {
  frameCount = 0;
  droppedInTrial = 0;
  lastFrameTs = null;
  frameIntervals = [];
  function step(ts) {
    if (!running) { monitorRaf = null; return; }
    if (lastFrameTs !== null) {
      var dt = ts - lastFrameTs;
      var nominal = deviceCal && deviceCal.frameIntervalMs ? deviceCal.frameIntervalMs : null;
      if (nominal) {
        var mult = Math.round(dt / nominal);
        if (mult >= 2) droppedInTrial += mult - 1;
      }
      // Bounded: a 40-minute trial would otherwise hold 144,000 numbers, and
      // a few thousand is already more than a MAD estimate needs.
      if (frameIntervals.length < 5000) frameIntervals.push(dt);
    }
    lastFrameTs = ts;
    frameCount += 1;
    monitorRaf = requestAnimationFrame(step);
  }
  monitorRaf = requestAnimationFrame(step);
}

function stopFrameMonitor() {
  if (monitorRaf) { cancelAnimationFrame(monitorRaf); monitorRaf = null; }
}

/* What the trial actually ran under, assembled at the end. */
function buildConditions() {
  var trialFrames = frameIntervals.length > 1 ? T.analyseFrames(
    frameIntervals.reduce(function (acc, dt) {
      acc.push(acc[acc.length - 1] + dt); return acc;
    }, [0])
  ) : null;

  return {
    environment: env,
    display: deviceCal,
    displayGrade: deviceCal ? T.gradeCalibration(deviceCal, null) : null,
    inputDelay: inputProbe ? inputProbe.result() : null,
    framesDuringTrial: frameCount,
    droppedFramesDuringTrial: droppedInTrial,
    trialFrameMadMs: trialFrames ? trialFrames.madIntervalMs : null,
    touchUsed: touchUsed,
    fullscreen: wasFullscreen
  };
}

function diagRows(c, opts) {
  var e = c.environment || {};
  var cal = c.display;
  var inp = c.inputDelay;
  var rows = [];

  rows.push([L.t('diag.rBrowser'), D.describeBrowser(e.userAgent)]);
  rows.push([L.t('diag.rScreen'),
    (e.screenW ? e.screenW + '×' + e.screenH : '—') +
    (e.devicePixelRatio ? '  ·  ' + L.t('diag.scale', { x: e.devicePixelRatio }) : '')]);
  rows.push([L.t('diag.rRefresh'), cal && cal.refreshHz
    ? cal.refreshHz.toFixed(2) + ' Hz  (' + cal.frameIntervalMs.toFixed(2) + ' ms)'
    : L.t('diag.notMeasured')]);
  rows.push([L.t('diag.rQuant'), cal && cal.frameIntervalMs
    ? '± ' + (cal.frameIntervalMs / 2).toFixed(2) + ' ms' : '—']);
  rows.push([L.t('diag.rJitter'), cal && cal.madIntervalMs != null
    ? cal.madIntervalMs.toFixed(2) + ' ms' : '—']);
  rows.push([L.t('diag.rTimer'), e.timerResolutionMs != null
    ? e.timerResolutionMs.toFixed(3) + ' ms' +
      (e.crossOriginIsolated ? '' : '  ·  ' + L.t('diag.notIsolated'))
    : '—']);
  // Three outcomes, and conflating them would mislead: measured, nothing to
  // measure yet, or a browser whose event clock cannot be trusted.
  var inputText;
  if (inp && inp.usable) {
    inputText = inp.medianDelayMs.toFixed(1) + ' ms  ·  ' + L.t('diag.nSamples', { n: inp.n }) +
                '  ·  ' + L.t('diag.notRemoved');
  } else if (inp && inp.reason === 'too few samples') {
    inputText = L.t('diag.inputWaiting', { n: inp.n });
  } else if (inp) {
    inputText = L.t('diag.inputUnusable');
  } else {
    inputText = L.t('diag.notMeasured');
  }
  rows.push([L.t('diag.rInput'), inputText]);

  if (opts && opts.trial) {
    // The same responses read off both clocks, so the size of the correction
    // is visible rather than assumed.
    var ev = (opts.rtEventMs || []).filter(function (v) { return v != null; });
    var hd = (opts.rtHandlerMs || []).filter(function (v) { return v != null; });
    var pairs = [];
    (opts.rtEventMs || []).forEach(function (v, i) {
      var h = (opts.rtHandlerMs || [])[i];
      if (v != null && h != null) pairs.push(h - v);
    });
    var meanOf = function (a) {
      return a.length ? a.reduce(function (x, y) { return x + y; }, 0) / a.length : null;
    };
    var scoredIsEvent = (opts.rtSource || 'event') === 'event';
    var mark = '  ← ' + L.t('diag.scoredMark');
    rows.push([L.t('diag.rRtEvent'), ev.length
      ? n1(meanOf(ev)) + ' ms' + (scoredIsEvent ? mark : '') : L.t('diag.notMeasured')]);
    rows.push([L.t('diag.rRtHandler'), hd.length
      ? n1(meanOf(hd)) + ' ms' + (scoredIsEvent ? '' : mark) : '—']);
    rows.push([L.t('diag.rRtDiff'), pairs.length
      ? n2(meanOf(pairs)) + ' ms  ·  ' + L.t('diag.nSamples', { n: pairs.length })
      : '—']);

    rows.push([L.t('diag.rDropped'), formatDroppedWeb(c)]);
    rows.push([L.t('diag.rFullscreen'), c.fullscreen ? L.t('diag.yes') : L.t('diag.no')]);
    rows.push([L.t('diag.rFocus'), opts.blurCount
      ? L.t('diag.focusLost', { n: opts.blurCount }) : L.t('diag.focusKept')]);
  }
  rows.push([L.t('diag.rDelivery'), describeDelivery(e.delivery)]);
  if (e.cores) rows.push([L.t('diag.rMachine'),
    L.t('diag.cores', { n: e.cores }) + (e.memoryGb ? '  ·  ' + e.memoryGb + ' GB' : '')]);
  return rows;
}

/* Dropped frames need their denominator to be readable — the same lesson as
 * the desktop build's counter. */
function formatDroppedWeb(c) {
  if (!c.framesDuringTrial) return '—';
  var expected = c.framesDuringTrial + c.droppedFramesDuringTrial;
  var pct = (c.droppedFramesDuringTrial / expected) * 100;
  var verdict = pct < 1 ? L.t('diag.dropOk')
              : pct < 5 ? L.t('diag.dropSome') : L.t('diag.dropBad');
  return c.droppedFramesDuringTrial.toLocaleString() + ' / ' + expected.toLocaleString() +
         '  (' + pct.toFixed(2) + '% — ' + verdict + ')';
}

function describeDelivery(d) {
  if (!d) return '—';
  if (d.local) return L.t('diag.deliveredLocal', { p: d.protocol });
  var bits = [d.protocol + '://' + (d.host || '?')];
  if (d.networkMs != null) bits.push(L.t('diag.loadedIn', { ms: Math.round(d.networkMs) }));
  if (d.connectionType) bits.push(d.connectionType);
  return bits.join('  ·  ') + '  ·  ' + L.t('diag.provenanceOnly');
}

function renderDiagTable(tbodyId, rows) {
  var tb = $(tbodyId);
  tb.innerHTML = '';
  rows.forEach(function (r) {
    var tr = document.createElement('tr');
    var th = document.createElement('th');
    th.textContent = r[0];
    var td = document.createElement('td');
    td.textContent = r[1];
    tr.appendChild(th); tr.appendChild(td);
    tb.appendChild(tr);
  });
}

function renderVerdict(elId, report) {
  var el = $(elId);
  el.className = 'note verdict ' + report.grade;
  var html = '<strong>' + L.t('diag.grade.' + report.grade) + '</strong>';
  if (report.notes.length) {
    html += '<ul>';
    report.notes.forEach(function (n) {
      html += '<li class="' + n.level + '">' + L.t(n.key, n.vars) + '</li>';
    });
    html += '</ul>';
  } else {
    html += ' ' + L.t('diag.allClear');
  }
  el.innerHTML = html;
}

/* The on-demand check from the setup screen: same measurement the trial makes,
 * run now so a machine can be judged before anyone is asked to sit down. */
function runDeviceCheck() {
  var btn = $('btnDeviceCheck');
  btn.disabled = true;
  var was = btn.textContent;
  btn.textContent = L.t('diag.running');

  var probe = T.makeInputProbe();
  function grab(e) { probe.record(e); }
  document.addEventListener('keydown', grab);
  document.addEventListener('pointerdown', grab);

  T.calibrateDisplay(80, 10).then(function (cal) {
    deviceCal = cal;
    env = D.snapshot();

    function draw() {
      var c = {
        environment: env, display: cal,
        inputDelay: probe.result(),
        framesDuringTrial: 0, droppedFramesDuringTrial: 0,
        fullscreen: !!document.fullscreenElement
      };
      renderDiagTable('deviceCheckTable', diagRows(c, { trial: false }));
      renderVerdict('deviceCheckVerdict',
        D.screeningReport(env, cal, c.inputDelay, { touchUsed: false }));
    }

    draw();
    $('deviceCheckResult').hidden = false;
    renderEnvLine();
    btn.disabled = false;
    btn.textContent = was;

    // Input delay needs real presses. Rather than blocking the check behind a
    // key-tapping ritual, the display results appear at once and the input row
    // fills itself in as soon as someone taps — and stops listening when the
    // trial starts, so it cannot compete with the trial's own probe.
    var redraw = function (e) { probe.record(e); draw(); };
    document.removeEventListener('keydown', grab);
    document.removeEventListener('pointerdown', grab);
    document.addEventListener('keydown', redraw);
    document.addEventListener('pointerdown', redraw);
    deviceCheckStop = function () {
      document.removeEventListener('keydown', redraw);
      document.removeEventListener('pointerdown', redraw);
      deviceCheckStop = null;
    };
  });
}

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
    var Ctx = window.AudioContext || window.webkitAudioContext;
    if (Ctx) audioCtx = new Ctx();
  } catch (e) { audioCtx = null; }
}

function startAlarm() {
  if (!cfg || !cfg.alarm || !audioCtx || alarmNodes) return;
  try {
    if (audioCtx.state === 'suspended') audioCtx.resume();
    var carrier = audioCtx.createOscillator();
    var gate = audioCtx.createGain();
    var lfo = audioCtx.createOscillator();
    var lfoDepth = audioCtx.createGain();

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
  var bar = $('alarmBar');
  if (bar) bar.hidden = true;
}

/* ---------------- trial ---------------- */

function parseIsiSet() {
  return $('isiSetMs').value.split(/[,\s]+/)
    .map(function (v) { return Math.round(parseFloat(v) * 1000); })
    .filter(function (v) { return isFinite(v) && v > 0; });
}

function validate() {
  // Re-run the identity and duplicate checks here rather than trusting the
  // last keystroke to have fired an event.
  checkParticipant();
  if (P.readDate(birthSelects()) === null) return L.t('err.birthIncomplete');
  if (identityDiff.length && !$('confirmIdentity').checked) return L.t('err.identityConflict');
  if (duplicateTrial && !$('confirmTrial').checked) return L.t('err.trialDuplicate');

  var mode = $('mode').value;
  var stim = intVal('stimMs', 1000);
  var shortestIsi;

  if (mode === 'pvt') {
    var set = parseIsiSet();
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
    rtSource: $('rtSource').value,          // 'event' | 'handler'
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
    birthDate: P.readDate(birthSelects()) || '',
    education: txtVal('pEducation'),
    sessionLabel: txtVal('sessionLabel') || 'NA',
    trialNumber: intVal('trialNumber', 1)
  };

  epochs = [];
  presses = [];
  // A fresh input probe per trial: dispatch delay is measured from the real
  // responses rather than a separate key-tapping ritual.
  inputProbe = T.makeInputProbe();
  touchUsed = false;
  rtStampRejected = 0;
  currentEpoch = null;
  consecutiveMisses = 0;
  missRunStartIndex = -1;
  blurCount = 0;
  running = false;
  kssBefore = null;
  kssAfter = null;

  if (deviceCheckStop) deviceCheckStop();
  initAudio();
  stopAlarm();
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

  /*
   * Measure the display while the countdown runs. Three seconds of otherwise
   * dead time is enough for ~90 frames at 60 Hz, and measuring here rather
   * than on the setup screen means the numbers describe the display in the
   * state the task will actually use — same window, same fullscreen mode.
   * A slow display can take longer than the countdown, so the task waits for
   * the measurement rather than starting without it.
   */
  var calDone = T.calibrateDisplay(80, 10).then(function (c) {
    deviceCal = c;
  }).catch(function () { deviceCal = null; });

  var iv = setInterval(function () {
    n -= 1;
    if (n <= 0) {
      clearInterval(iv);
      $('countdownNum').textContent = '';
      calDone.then(startTask);
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
  startFrameMonitor();
  startClock();
  // Schedule the first stimulus rather than firing it: onsets[0] is the lead-in,
  // so the screen stays blank for one interval after the countdown.
  scheduleEpoch(0);
}

/*
 * Absolute targets taken from the pre-built schedule, so presentation does not
 * drift across a long trial and a variable-interval (PVT) schedule is delivered
 * exactly as generated.
 */
function scheduleEpoch(n) {
  // Onsets already include the lead-in, so the first stimulus is one interval
  // after the countdown rather than the instant it ends. Past the final
  // stimulus there is no onset to aim at, but the trial still has to be ended,
  // so the last call is scheduled at the planned end of the trial.
  var target = n < cfg.schedule.nStimuli
    ? cfg.schedule.onsets[n]
    : cfg.schedule.plannedDurationMs;
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
    // Taken from the schedule, not re-derived from the measured onset, so a
    // stimulus landing a millisecond either side of a boundary cannot slide
    // into the neighbouring minute.
    minute: cfg.schedule.minutes[n],
    // epochIsiMs is this stimulus's response window (time to the next one);
    // isiBeforeMs is the wait that preceded it, defined for the first one too.
    epochIsiMs: cfg.schedule.epochIsi[n],
    isiBeforeMs: cfg.schedule.isiBefore[n],
    responded: false,
    rtMs: null,
    rtEventMs: null,
    rtHandlerMs: null,
    rtUsed: null,
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
function handleResponse(evt) {
  if (!running) return;
  // Measured first, so the delta reflects handler entry rather than whatever
  // work happens below.
  if (evt && inputProbe) inputProbe.record(evt);
  if (evt && evt.pointerType && evt.pointerType !== 'mouse') touchUsed = true;

  /*
   * Both clock readings, for every response, exactly as the desktop build
   * does — see its README section on input dispatch.
   *
   * evt.timeStamp is when the platform delivered the event; performance.now()
   * here is when this code got to look at it. The gap is the dispatch delay.
   * Which one is scored is a setting; both are always recorded so the
   * correction can be checked, reproduced or undone rather than trusted.
   *
   * The stamp is validated per response rather than once per trial: a value
   * that is missing, or that implies a delay outside a plausible range, means
   * this browser is not putting event.timeStamp on the performance clock, and
   * that response silently falls back to handler time — recorded, not hidden.
   */
  var tHandler = performance.now();
  var tEvent = (evt && typeof evt.timeStamp === 'number' && evt.timeStamp > 0)
    ? evt.timeStamp : null;
  var delta = tEvent === null ? null : tHandler - tEvent;
  var stampUsable = delta !== null && delta >= -1 && delta < 150;
  var useEvent = cfg.rtSource === 'event' && stampUsable;
  var now = useEvent ? tEvent : tHandler;
  if (!stampUsable && tEvent !== null) rtStampRejected += 1;

  // Every press is logged, including presses outside any epoch, because the
  // integrity check needs the full pattern — not just the scored responses.
  presses.push({ tMs: now - startTime, epochIndex: currentEpoch ? currentEpoch.index : null });

  if (!currentEpoch) return;
  if (currentEpoch.responded) { currentEpoch.extra += 1; return; }

  currentEpoch.responded = true;
  currentEpoch.rtMs = now - currentEpoch.onsetPerf;
  currentEpoch.rtEventMs = stampUsable ? tEvent - currentEpoch.onsetPerf : null;
  currentEpoch.rtHandlerMs = tHandler - currentEpoch.onsetPerf;
  currentEpoch.rtUsed = useEvent ? 'event' : 'handler';

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
  stopFrameMonitor();
  // Read the fullscreen state BEFORE leaving it, or the record would always
  // say the task ran in a window.
  wasFullscreen = !!document.fullscreenElement;
  ledOff();
  releaseWakeLock();
  if (document.fullscreenElement && document.exitFullscreen) document.exitFullscreen().catch(function () {});

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
  var scored = S.score(epochs.map(function (e) {
    return {
      index: e.index, onsetMs: e.onsetMs, rtMs: e.rtMs,
      block: e.block, minute: e.minute, epochIsiMs: e.epochIsiMs, isiBeforeMs: e.isiBeforeMs
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
    rtEventMs: epochs.map(function (e) { return e.rtEventMs; }),
    rtHandlerMs: epochs.map(function (e) { return e.rtHandlerMs; }),
    rtStampRejected: rtStampRejected,
    conditions: buildConditions(),
    norms: S.normativeReport({
      trials: scored.trials,
      elapsedMs: epochs.length ? epochs[epochs.length - 1].onsetMs + epochs[epochs.length - 1].epochIsiMs : 0,
      minuteBuckets: scored.perMinute.length,
      mode: cfg.mode,
      isiMs: cfg.isiMs,
      hour: parseInt(meta.time.slice(0, 2), 10),
      lapseMs: cfg.lapseMs,
      falseStartMs: cfg.falseStartMs,
      missCriterion: cfg.missCriterion
    }, window.BSRTNorms),
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

  renderConditions(r);
  renderNorms(r);
  renderIntegrity(r);
  renderKss(r);
  renderPerMinute(r);
}

/* ---------------- normative comparison ---------------- */

/*
 * Labels reuse the keys the rest of the results screen already uses, so the
 * panel speaks the interface language without another 84 strings to translate.
 * The reference table's own English labels stay in norms.js as the data's
 * documentation.
 */
var NORM_LABEL = {
  trials: 'results.totalTrials', hitRatio: 'results.hitRatio', hits: 'results.hits',
  misses: 'results.misses', lapses: 'results.lapses', falseStarts: 'results.falseStarts',
  ep12: 'results.ep12', ep36: 'results.ep36', ep7: 'results.ep7',
  rtMean: 'results.average', rtMedian: 'results.median', rtSd: 'results.sd',
  rtFast10: 'results.fastest', rtSlow10: 'results.slowest', rtIpr: 'results.ipr',
  rsMean: 'results.average', rsMedian: 'results.median', rsSd: 'results.sd',
  rsFast10: 'results.fastest', rsSlow10: 'results.slowest', rsIpr: 'results.ipr'
};
var NORM_SECTION = {
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

/*
 * The norms panel is off by default: the results screen shows the plain,
 * uncoloured tables unless someone asks for the comparison. The choice is
 * remembered, so a researcher who wants norms does not re-tick it every trial.
 * The report itself is always computed and always exported — only the display
 * is optional.
 */
var SHOW_NORMS_KEY = 'bsrt.showNorms.v1';

function normsWanted() {
  try { return localStorage.getItem(SHOW_NORMS_KEY) === '1'; } catch (e) { return false; }
}

function applyNormsVisibility() {
  var on = $('showNorms').checked;
  $('normsPanel').hidden = !on;
  try { localStorage.setItem(SHOW_NORMS_KEY, on ? '1' : '0'); } catch (e) { /* private mode */ }
}

function initNormsToggle() {
  $('showNorms').checked = normsWanted();
  $('normsPanel').hidden = !$('showNorms').checked;
  $('showNorms').addEventListener('change', applyNormsVisibility);
}

function renderConditions(r) {
  var c = r.conditions;
  if (!c) { $('diagTable').innerHTML = ''; $('diagVerdict').textContent = ''; return; }
  renderDiagTable('diagTable', diagRows(c, {
    trial: true, blurCount: r.blurCount,
    rtEventMs: r.rtEventMs, rtHandlerMs: r.rtHandlerMs,
    rtSource: r.config && r.config.rtSource
  }));
  renderVerdict('diagVerdict', D.screeningReport(c.environment, c.display, c.inputDelay, {
    touchUsed: c.touchUsed, blurCount: r.blurCount, fullscreen: c.fullscreen
  }));
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
  var g = r.scored.integrity;
  $('inPresses').textContent = g.totalPresses;
  $('inExtra').textContent = g.extraRate
    ? L.t('results.extraVal', { n: g.extraPresses, pct: (g.extraRate * 100).toFixed(0) })
    : g.extraPresses;
  $('inBurst').textContent = L.t('results.burstVal', { n: g.burstMax, ms: g.thresholds.burstWindowMs });
  $('inRapid').textContent = L.t('results.rapidVal', { n: g.rapidPairs, ms: g.thresholds.rapidGapMs });

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
  if (c.correction === 'falseStarts') return L.t('corr.descFs', { ms: c.falseStartMs });
  if (c.correction === 'outliers') return L.t('corr.descOut', { sd: c.sdMultiplier });
  if (c.correction === 'both') return L.t('corr.descBoth', { ms: c.falseStartMs, sd: c.sdMultiplier });
  return L.t('corrections.none');
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
    ? L.t('results.noTrend')
    : L.t('results.trend', { v: n3(r.scored.dynamicsRs.slope) });
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
  // rt_event_ms and rt_handler_ms are the same response read off two clocks;
  // input_delay_ms is their difference and rt_source_used says which was scored.
  'responded', 'rt_event_ms', 'rt_handler_ms', 'input_delay_ms', 'rt_source_used',
  'rt_ms', 'rs_per_sec',
  'outcome', 'lapse', 'late_response', 'false_start', 'extra_responses'
]);

function rawRows(r) {
  var pv = participantVals(r);
  var extras = r.extras || [];
  var evRts = r.rtEventMs || [];
  var hdRts = r.rtHandlerMs || [];
  var used = (r.scored.trials || []).map(function (t, i) {
    if (evRts[i] == null && hdRts[i] == null) return '';
    return (r.config && r.config.rtSource === 'handler') || evRts[i] == null ? 'handler' : 'event';
  });
  return r.scored.trials.map(function (t, i) {
    return pv.concat([
      // 1-based on the way out, like block and minute: epoch_index was the
      // only column that still started at zero.
      t.index + 1, t.block === null ? '' : t.block + 1, t.minute + 1, round(t.onsetMs, 2),
      t.epochIsiMs, t.isiBeforeMs,
      t.rtMs === null ? 0 : 1,
      round(evRts[i], 3), round(hdRts[i], 3),
      (evRts[i] == null || hdRts[i] == null) ? null : round(hdRts[i] - evRts[i], 3),
      used[i] || '',
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

/*
 * Normative comparison, long format: one row per variable so files rbind
 * cleanly, the same shape as the other exports. Values here follow the
 * REFERENCE workbook's definitions, not the app's own scoring — see
 * normativeSummary in scoring.js — so norm_value will not always equal the
 * matching column in the summary export. Both are kept on purpose.
 */
var NORMS_HEADER = PARTICIPANT_COLS.concat([
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
  'alarm_enabled', 'extra_responses', 'page_blur_count',

  'device_browser', 'device_platform', 'screen_w', 'screen_h', 'device_pixel_ratio',
  'refresh_hz_measured', 'frame_interval_ms', 'onset_quantisation_ms', 'frame_mad_ms',
  'timer_resolution_ms', 'cross_origin_isolated',
  'input_dispatch_median_ms', 'input_dispatch_mad_ms', 'input_dispatch_n', 'input_stamp_usable',
  'rt_source_setting', 'mean_rt_event_ms', 'mean_rt_handler_ms', 'mean_input_delay_ms',
  'median_input_delay_ms', 'n_rt_compared', 'n_stamp_rejected',
  'frames_trial', 'dropped_frames_trial', 'dropped_rate_trial',
  'ran_fullscreen', 'touch_used', 'cpu_cores', 'device_memory_gb',
  'page_protocol', 'page_host', 'page_local', 'page_network_ms', 'connection_type',
  'device_grade',]);

/* Device conditions for the summary row. Kept in one place so the header and
 * the values cannot drift apart. */
/* The two clock readings summarised, so the size of the input-delay
 * correction sits on the same row as the metrics it affected. */
function rtSourceVals(r) {
  var ev = (r.rtEventMs || []).filter(function (v) { return v != null; });
  var hd = (r.rtHandlerMs || []).filter(function (v) { return v != null; });
  var pairs = [];
  (r.rtEventMs || []).forEach(function (v, i) {
    var h = (r.rtHandlerMs || [])[i];
    if (v != null && h != null) pairs.push(h - v);
  });
  var mean = function (a) {
    return a.length ? a.reduce(function (x, y) { return x + y; }, 0) / a.length : null;
  };
  return [
    (r.config && r.config.rtSource) || 'event',
    round(mean(ev), 3), round(mean(hd), 3), round(mean(pairs), 4),
    pairs.length ? round(T.stats.median(pairs), 4) : null,
    pairs.length, r.rtStampRejected || 0
  ];
}

function conditionVals(r) {
  var c = r.conditions || {};
  var e = c.environment || {};
  var d = e.delivery || {};
  var cal = c.display || {};
  var inp = c.inputDelay || {};
  var frames = c.framesDuringTrial || 0;
  var dropped = c.droppedFramesDuringTrial || 0;
  var expected = frames + dropped;
  return [
    D.describeBrowser(e.userAgent), e.platform, e.screenW, e.screenH, e.devicePixelRatio,
    round(cal.refreshHz, 3), round(cal.frameIntervalMs, 4),
    cal.frameIntervalMs ? round(cal.frameIntervalMs / 2, 4) : null,
    round(cal.madIntervalMs, 4),
    round(e.timerResolutionMs, 5), e.crossOriginIsolated == null ? '' : (e.crossOriginIsolated ? 1 : 0),
    inp.usable ? round(inp.medianDelayMs, 3) : null,
    inp.usable ? round(inp.madDelayMs, 3) : null,
    inp.n == null ? null : inp.n,
    inp.usable == null ? '' : (inp.usable ? 1 : 0),
    ...rtSourceVals(r),
    frames || null, dropped,
    expected ? round(dropped / expected, 5) : null,
    c.fullscreen ? 1 : 0, c.touchUsed ? 1 : 0,
    e.cores, e.memoryGb,
    d.protocol, d.host, d.local == null ? '' : (d.local ? 1 : 0),
    round(d.networkMs, 1), d.connectionType,
    c.environment ? D.screeningReport(c.environment, c.display, c.inputDelay,
      { touchUsed: c.touchUsed, blurCount: r.blurCount, fullscreen: c.fullscreen }).grade : null
  ];
}

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
  ]).concat(conditionVals(r));
}

function allRows(builder, header) {
  var all = loadSessions();
  if (!all.length) { alert(L.t('export.none')); return null; }
  var rows = [];
  for (var i = 0; i < all.length; i++) rows = rows.concat(builder(all[i]));
  return toCsv(header, rows);
}

/* ---------------- wiring ---------------- */

$('btnStart').addEventListener('click', beginTrial);

$('btnStopAlarm').addEventListener('click', stopAlarm);

$('btnAgain').addEventListener('click', function () {
  stopAlarm();
  // Straight to the next free number for this participant and session, so a
  // second run cannot land on a trial number that already holds data.
  suggestTrialNumber();
  checkParticipant();
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

$('btnExportNorms').addEventListener('click', function () {
  if (lastResult) download(fileStem(lastResult) + '_norms.csv', toCsv(NORMS_HEADER, normsRows(lastResult)));
});
$('btnExportAllNorms').addEventListener('click', function () {
  var csv = allRows(normsRows, NORMS_HEADER);
  if (csv) download('bsrt_all_norms.csv', csv);
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
  if (!all.length) { alert(L.t('export.none')); return; }
  download('bsrt_all_summaries.csv', toCsv(SUMMARY_HEADER, all.map(summaryRow)));
});

$('btnAbort').addEventListener('click', function (e) {
  e.stopPropagation();
  if (running) endTask('aborted');
});

$('screen-task').addEventListener('pointerdown', function (e) {
  if (e.target && e.target.id === 'btnAbort') return;
  handleResponse(e);
});

document.addEventListener('keydown', function (e) {
  if (!running) return;
  if (e.code === 'Space' || e.key === ' ' || e.code === 'Enter' || e.key === 'Enter') {
    e.preventDefault();
    handleResponse(e);
  }
});

document.addEventListener('visibilitychange', function () {
  if (running && document.hidden) blurCount += 1;
  if (!document.hidden && running) requestWakeLock();
});

/* Suggest the short-protocol lapse threshold when the ceiling drops below 5 min. */
$('maxMinutes').addEventListener('change', function () {
  var short = intVal('maxMinutes', 40) < 5;
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
  var mode = $('mode').value;
  document.body.setAttribute('data-mode', mode);
  // A PVT runs for a fixed duration; the sleep-onset criterion is an OSLER idea.
  $('criterionOn').checked = mode !== 'pvt';
  $('modeNote').textContent = L.t(mode === 'pvt' ? 'mode.notePvt' : 'mode.noteBsrt');
  $('instructionsText').textContent = L.t(mode === 'pvt' ? 'instructions.pvt' : 'instructions.bsrt');
  $('taskHint').textContent = L.t(mode === 'pvt' ? 'task.hintPvt' : 'task.hintBsrt');
}

$('mode').addEventListener('change', applyMode);
$('language').addEventListener('change', applyLanguage);
applyLanguage();
refreshCount();
initParticipants();
initNormsToggle();
initDiagnostics();
$('btnDeviceCheck').addEventListener('click', runDeviceCheck);
