'use strict';

/* Participant identity — shared byte-for-byte between the browser build
 * (BSRT/) and the desktop build (BSRT-desktop/renderer/). `npm run
 * check:shared` in BSRT-desktop refuses to package if the copies drift.
 *
 * Two jobs.
 *
 * 1. BIRTH DATES IN ONE FORMAT. A free-text field collects 1985-03-07,
 *    07/03/1985, 3/7/85 and "March 1985" from four different testers, and the
 *    resulting column cannot be parsed. Three dropdowns make an out-of-format
 *    answer impossible to produce: there is nothing to type. They are used in
 *    preference to a native calendar popup because a birth date sits decades
 *    back — picking 1985 from a list is one click, where a calendar means
 *    paging through months. The stored value is always ISO YYYY-MM-DD, and
 *    month names are rendered in the interface language via Intl, so nothing
 *    has to be translated by hand.
 *
 * 2. RETURNING PARTICIPANTS. Session two should not mean retyping an ID and a
 *    birth date — retyping is exactly where a P001/POO1 or a transposed year
 *    creeps in and silently splits one participant into two. Details are kept
 *    in a small local roster and recalled from a dropdown.
 *
 *    The roster is also what makes accidental overwriting detectable: once a
 *    participant has a saved record, any change to their identity is a
 *    conflict that the app reports field by field and refuses to act on until
 *    someone says which version is right. The same applies to re-running a
 *    trial number that already has data filed under it.
 */

(function (global) {

  var MAX_AGE_YEARS = 120;

  /* Fields that identify a person and are therefore protected. The session
   * label and trial number are deliberately absent: those are expected to
   * change from one visit to the next. */
  var IDENTITY_FIELDS = ['name', 'birthDate', 'address', 'education'];

  /* ---------------- birth date ---------------- */

  function daysInMonth(year, month) {
    if (!month) return 31;
    // With no year chosen yet, offer 29 February rather than hide a valid day.
    if (!year) return month === 2 ? 29 : new Date(2001, month, 0).getDate();
    return new Date(year, month, 0).getDate();
  }

  function monthNames(lang) {
    var out = [];
    var fmt = null;
    try {
      fmt = new Intl.DateTimeFormat(lang || 'en', { month: 'long' });
    } catch (e) { fmt = null; }
    for (var m = 1; m <= 12; m++) {
      var label = null;
      if (fmt) {
        try { label = fmt.format(new Date(2001, m - 1, 1)); } catch (e) { label = null; }
      }
      out.push(label || String(m));
    }
    return out;
  }

  function option(value, label) {
    var o = document.createElement('option');
    o.value = value;
    o.textContent = label;
    return o;
  }

  function fill(sel, placeholder, values, labels) {
    var keep = sel.value;
    sel.innerHTML = '';
    sel.appendChild(option('', placeholder));
    for (var i = 0; i < values.length; i++) {
      sel.appendChild(option(String(values[i]), labels[i]));
    }
    // Restore the previous choice when it still exists (31 January stays
    // selected through a re-render; 31 February cannot and is dropped).
    sel.value = keep;
    if (sel.value !== keep) sel.value = '';
  }

  /* Populate the three selects. Called again on every language change and
   * whenever month or year moves, so the day list matches the month length. */
  function renderDateSelects(els, opts) {
    opts = opts || {};
    var lang = opts.lang || 'en';
    var ph = opts.placeholders || { day: 'Day', month: 'Month', year: 'Year' };
    var thisYear = (opts.today || new Date()).getFullYear();

    var years = [], yl = [];
    for (var y = thisYear; y >= thisYear - MAX_AGE_YEARS; y--) { years.push(y); yl.push(String(y)); }
    fill(els.year, ph.year, years, yl);

    var months = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
    fill(els.month, ph.month, months, monthNames(lang));

    var n = daysInMonth(parseInt(els.year.value, 10) || null,
                        parseInt(els.month.value, 10) || null);
    var days = [], dl = [];
    for (var d = 1; d <= n; d++) { days.push(d); dl.push(String(d)); }
    fill(els.day, ph.day, days, dl);
  }

  function pad(n) { return (n < 10 ? '0' : '') + n; }

  /* '' when nothing has been chosen, null when the answer is half finished,
   * an ISO date when it is complete. The three states are distinct: blank is
   * a legitimate answer for an anonymised participant, half finished is not. */
  function readDate(els) {
    var d = parseInt(els.day.value, 10);
    var m = parseInt(els.month.value, 10);
    var y = parseInt(els.year.value, 10);
    var chosen = (els.day.value ? 1 : 0) + (els.month.value ? 1 : 0) + (els.year.value ? 1 : 0);
    if (chosen === 0) return '';
    if (chosen < 3) return null;
    if (d > daysInMonth(y, m)) return null;
    return y + '-' + pad(m) + '-' + pad(d);
  }

  function writeDate(els, iso) {
    var parts = /^(\d{4})-(\d{2})-(\d{2})$/.exec(String(iso || ''));
    if (!parts) {
      els.day.value = ''; els.month.value = ''; els.year.value = '';
      return false;
    }
    els.year.value = String(parseInt(parts[1], 10));
    els.month.value = String(parseInt(parts[2], 10));
    return true;   // the caller re-renders, then sets the day
  }

  function setDay(els, iso) {
    var parts = /^(\d{4})-(\d{2})-(\d{2})$/.exec(String(iso || ''));
    els.day.value = parts ? String(parseInt(parts[3], 10)) : '';
  }

  /* ---------------- roster ---------------- */

  function keyOf(id) { return String(id == null ? '' : id).trim().toLowerCase(); }

  function loadProfiles(storage, key) {
    try {
      var raw = JSON.parse(storage.getItem(key) || '{}');
      return raw && typeof raw === 'object' && !Array.isArray(raw) ? raw : {};
    } catch (e) { return {}; }
  }

  function saveProfiles(storage, key, map) {
    try { storage.setItem(key, JSON.stringify(map)); return true; }
    catch (e) { return false; }
  }

  function upsertProfile(map, meta, now) {
    var k = keyOf(meta.participantId);
    if (!k) return map;
    var prev = map[k] || {};
    map[k] = {
      participantId: meta.participantId,
      name: meta.name || '',
      birthDate: meta.birthDate || '',
      address: meta.address || '',
      education: meta.education || '',
      lastSessionLabel: meta.sessionLabel || '',
      sessions: (prev.sessions || 0) + 1,
      firstSeen: prev.firstSeen || (now || new Date().toISOString()),
      updatedAt: now || new Date().toISOString()
    };
    return map;
  }

  /* Existing installations already have trials on disk. Derive a roster from
   * them once, so upgrading does not present an empty participant list. */
  function seedFromSessions(map, sessions) {
    for (var i = 0; i < sessions.length; i++) {
      var p = sessions[i] && sessions[i].participant;
      if (!p) continue;
      var k = keyOf(p.participantId);
      if (!k || k === 'na' || map[k]) continue;
      map[k] = {
        participantId: p.participantId,
        name: p.name || '',
        birthDate: p.birthDate || '',
        address: p.address || '',
        education: p.education || '',
        lastSessionLabel: p.sessionLabel || '',
        sessions: 0,
        firstSeen: p.date || '',
        updatedAt: p.date || ''
      };
    }
    return map;
  }

  function listProfiles(map) {
    var out = [];
    for (var k in map) if (Object.prototype.hasOwnProperty.call(map, k)) out.push(map[k]);
    out.sort(function (a, b) {
      return String(a.participantId).localeCompare(String(b.participantId), undefined, { numeric: true });
    });
    return out;
  }

  /* Every identity field on which the form and the saved record disagree.
   * Clearing a field counts: losing a stored birth date by accident is as
   * damaging as replacing it with the wrong one. */
  function diffIdentity(stored, entered) {
    var out = [];
    if (!stored) return out;
    for (var i = 0; i < IDENTITY_FIELDS.length; i++) {
      var f = IDENTITY_FIELDS[i];
      var from = String(stored[f] == null ? '' : stored[f]).trim();
      var to = String(entered[f] == null ? '' : entered[f]).trim();
      if (from !== to) out.push({ field: f, from: from, to: to });
    }
    return out;
  }

  /* ---------------- trial numbering ---------------- */

  function sameTrial(p, id, label, n) {
    return keyOf(p.participantId) === keyOf(id) &&
           keyOf(p.sessionLabel) === keyOf(label) &&
           Number(p.trialNumber) === Number(n);
  }

  function findDuplicate(sessions, id, label, n) {
    for (var i = sessions.length - 1; i >= 0; i--) {
      var p = sessions[i] && sessions[i].participant;
      if (p && sameTrial(p, id, label, n)) return sessions[i];
    }
    return null;
  }

  function nextTrialNumber(sessions, id, label) {
    var highest = 0;
    for (var i = 0; i < sessions.length; i++) {
      var p = sessions[i] && sessions[i].participant;
      if (!p) continue;
      if (keyOf(p.participantId) !== keyOf(id)) continue;
      if (keyOf(p.sessionLabel) !== keyOf(label)) continue;
      var n = Number(p.trialNumber);
      if (isFinite(n) && n > highest) highest = n;
    }
    return highest + 1;
  }

  global.BSRTParticipants = {
    MAX_AGE_YEARS: MAX_AGE_YEARS,
    IDENTITY_FIELDS: IDENTITY_FIELDS,
    daysInMonth: daysInMonth,
    monthNames: monthNames,
    renderDateSelects: renderDateSelects,
    readDate: readDate,
    writeDate: writeDate,
    setDay: setDay,
    keyOf: keyOf,
    loadProfiles: loadProfiles,
    saveProfiles: saveProfiles,
    upsertProfile: upsertProfile,
    seedFromSessions: seedFromSessions,
    listProfiles: listProfiles,
    diffIdentity: diffIdentity,
    findDuplicate: findDuplicate,
    nextTrialNumber: nextTrialNumber
  };

})(typeof window !== 'undefined' ? window : globalThis);
