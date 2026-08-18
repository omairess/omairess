#!/usr/bin/env python
"""BSRT / PVT for PsychoPy.

    python bsrt_psychopy.py                 run it (a dialog collects everything)
    python bsrt_psychopy.py --simulate      run the whole pipeline with no screen

Open this in the PsychoPy **Coder** view and press Run, or run it from a shell
in an environment that has PsychoPy installed. There is no Builder .psyexp: the
task is a timed loop with a stopping rule, which Builder expresses awkwardly and
which has to stay identical to the browser and desktop builds.

WHY THIS BUILD EXISTS. It measures a stimulus onset as the vsync timestamp of
the frame that carried it, and a reaction time against the keyboard's own
hardware stamp. The browser can only approximate both. If you have PsychoPy set
up already, this is the most accurate of the three.

WHAT IS SHARED. Scheduling, scoring and the normative comparison come from
bsrt_core.py, which is a checked port of scoring.js — tests/test_equivalence.py
drives both and asserts they agree to 1e-9. The CSV columns are the browser
build's columns, asserted by tests/test_headers.py. The on-screen text and the
KSS anchors are extracted from i18n.js by tools/gen_strings.py. Nothing that
matters is retyped here.

CONTROLS
    space / any key   respond
    escape            abort the trial (the data collected so far is still saved)
    s                 silence the sleep-onset alarm
"""

import argparse
import datetime
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import bsrt_io as bio          # noqa: E402
import bsrt_ui as ui           # noqa: E402
import bsrt_norms as bnorms    # noqa: E402
import bsrt_task as btask      # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
STRINGS = json.load(open(os.path.join(HERE, 'strings.json')))


def T(lang, key, **kw):
    """One on-screen string, in the chosen language.

    strings.json is generated from i18n.js by tools/gen_strings.py, so the
    wording is the browser build's wording and the four languages stay in step.
    Falls back to English rather than showing a bare key, matching i18n.js.
    """
    table = STRINGS.get(lang) or STRINGS['en']
    text = table.get(key) or STRINGS['en'].get(key) or key
    for k, v in kw.items():
        text = text.replace('{%s}' % k, str(v))
    return text


DEFAULTS = {
    'participant_id': '',
    'session_label': 'baseline',
    'trial_number': 1,
    'name': '',
    'birth_date': '',
    'education': '',
    'language': 'en',
    'mode': 'bsrt',
    'isi_ms': 3000,
    'isi_set_s': '2/4/6/8/10',
    'stim_ms': 1000,
    'hit_window_ms': 1000,
    'lapse_ms': 500,
    'miss_criterion': 7,
    'max_minutes': 40,
    'kss_when': 'none',
    'alarm': True,
    'fullscreen': True,
    'remove_false_starts': False,
    'remove_outliers': False,
    'out_dir': 'data',
    'seed': 0,
}

# What changes when the paradigm changes. The browser build ticks the
# sleep-onset criterion off for a PVT (app.js: criterionOn.checked = mode !==
# 'pvt'); without this the PsychoPy build would run a PVT that stops early at
# seven consecutive misses, which is not the same test and not what the other
# two builds do.
MODE_DEFAULTS = {
    'bsrt': {'miss_criterion': 7, 'isi_ms': 3000, 'max_minutes': 40},
    'pvt': {'miss_criterion': 0, 'isi_ms': 3000, 'max_minutes': 10},
}


def parse_isi_set(text):
    """'2/4/6/8/10' -> [2000, 4000, 6000, 8000, 10000] milliseconds.

    Accepts slashes, commas or spaces. Falls back to the standard PVT set
    rather than raising, because an unparseable box should not lose a
    participant who is already sitting down.
    """
    out = []
    for part in str(text).replace(',', ' ').replace('/', ' ').split():
        try:
            v = float(part)
        except ValueError:
            continue
        if v > 0:
            out.append(int(round(v * 1000)))
    return out or [2000, 4000, 6000, 8000, 10000]


RESPONSE_KEYS = ['space', 'return', 'lshift', 'rshift', 'a', 'l', 'num_1']
QUIT_KEYS = ['escape']
SILENCE_KEYS = ['s']


def cfg_from(settings):
    """Turn the dialog's flat answers into the config the engine takes."""
    return {
        'mode': settings['mode'],
        'isiMs': int(settings['isi_ms']),
        'stimMs': int(settings['stim_ms']),
        'hitWindowMs': int(settings['hit_window_ms']),
        'lapseMs': int(settings['lapse_ms']),
        'missCriterion': int(settings['miss_criterion']),
        'maxMinutes': int(settings['max_minutes']),
        'isiSetMs': parse_isi_set(settings.get('isi_set_s', '2/4/6/8/10')),
        'blockMs': 30000,
        'seed': int(settings['seed']) or random.randint(1, 2 ** 31 - 1),
        'removeFalseStarts': bool(settings['remove_false_starts']),
        'removeOutliers': bool(settings['remove_outliers']),
        'falseStartMs': 100,
        'sdMultiplier': 2,
        'alarm': bool(settings['alarm']),
        'fullscreen': bool(settings['fullscreen']),
        'feedbackMs': 500,
    }


# --------------------------------------------------------------------------
# The real display


class PsychoPyDisplay(object):
    """The Display the engine drives, backed by a PsychoPy window.

    Onsets are the times `win.flip()` reports, which is when the frame reached
    the screen. Responses are the keyboard's own timestamps, on the same clock.
    Neither is re-derived from a software timer, which is the whole point.
    """

    def __init__(self, win, kb, mode, frame_ms):
        from psychopy import visual

        self.win = win
        self.kb = kb
        self.frame = frame_ms
        self.quit = False
        self.silence = False

        # The same red as the browser build's LED (styles.css --led), and the
        # same neutral white for the PVT counter.
        self.dot = visual.Circle(win, radius=0.03, units='height',
                                 fillColor=ui.LED, colorSpace='rgb255',
                                 lineColor=None)
        self.counter = visual.TextStim(win, text='', height=0.12, color=ui.TEXT,
                                       colorSpace='rgb255', units='height')
        self.banner = None

    def frame_interval_ms(self):
        return self.frame

    def flip(self, state):
        kind = state[0]
        if kind == 'dot':
            self.dot.draw()
        elif kind == 'counter':
            # Whole milliseconds only: the counter is feedback, not a clock, and
            # a flickering final digit is a distraction.
            self.counter.text = '%d' % int(state[1])
            self.counter.draw()
        elif kind == 'feedback':
            self.counter.text = '%d ms' % int(round(state[1]))
            self.counter.draw()
        if self.banner is not None:
            self.banner.draw()
        t = self.win.flip()
        return _ms(t)

    def poll(self):
        out = []
        for k in self.kb.getKeys(keyList=None, waitRelease=False, clear=True):
            name = getattr(k, 'name', None)
            if name in QUIT_KEYS:
                self.quit = True
                continue
            if name in SILENCE_KEYS:
                self.silence = True
                continue
            out.append((_ms(k.tDown), name))
        return out

    def quit_requested(self):
        return self.quit


def _ms(t):
    """PsychoPy reports seconds; everything here is milliseconds."""
    return None if t is None else float(t) * 1000.0


# --------------------------------------------------------------------------
# A display with no screen, so the whole program can be exercised


class SimulatedDisplay(object):
    """Drives the same engine on a virtual clock, for --simulate.

    Not a test double for its own sake: it is what lets the dialog-to-CSV path
    be run on a machine with no PsychoPy and no display, which is how this file
    was checked.
    """

    def __init__(self, rts, frame_ms=1000.0 / 60.0):
        self.t = 0.0
        self.frame = frame_ms
        self.rts = list(rts)
        self._epoch = -1
        self._onset = None
        self._pending = []

    def frame_interval_ms(self):
        return self.frame

    def flip(self, state):
        self.t += self.frame
        showing = state[0] in ('dot', 'counter')
        if showing and self._onset is None:
            self._epoch += 1
            self._onset = self.t
            rt = self.rts[self._epoch] if self._epoch < len(self.rts) else None
            if rt is not None:
                self._pending.append(self._onset + rt)
        elif not showing and state[0] != 'feedback':
            self._onset = None
        return self.t

    def poll(self):
        out = [(t, 'space') for t in self._pending if t <= self.t]
        self._pending = [t for t in self._pending if t > self.t]
        return out

    def quit_requested(self):
        return False


# --------------------------------------------------------------------------
# Screens


def _wait_key(win, kb, draw, accept=None):
    """Draw a screen every frame and return the first accepted key.

    Redrawing and flipping each pass rather than spinning on getKeys keeps the
    window responsive to the OS and stops the loop pegging a core.
    """
    kb.clearEvents()
    while True:
        draw()
        win.flip()
        for k in kb.getKeys(waitRelease=False):
            if k.name in QUIT_KEYS:
                return 'escape'
            if accept is None or k.name in accept:
                return k.name


def show_instructions(win, ui_, kb, lang, mode):
    hint = T(lang, 'task.hintPvt' if mode == 'pvt' else 'task.hintBsrt')
    title = T(lang, 'pp.titlePvt' if mode == 'pvt' else 'pp.titleBsrt')

    def draw():
        ui.draw_instructions(ui_, title, [hint], T(lang, 'pp.anyKey'))
    return _wait_key(win, kb, draw)


def ask_kss(win, ui_, kb, lang, title):
    """The 9-point Karolinska Sleepiness Scale as a scale, not a list.

    Anchors come from strings.json, extracted from i18n.js, so the wording is
    the browser build's wording. Highlighting follows the number under the
    cursor keys as well as the digits, because a row of boxes invites arrows.
    """
    anchors = (STRINGS.get(lang) or STRINGS['en'])['_anchors']
    state = {'sel': 0}

    def draw():
        ui.draw_kss(ui_, title, T(lang, 'kss.question'), anchors, state['sel'],
                    T(lang, 'pp.kssKeys'))

    kb.clearEvents()
    while True:
        draw()
        win.flip()
        for k in kb.getKeys(waitRelease=False):
            name = k.name or ''
            if name in QUIT_KEYS:
                return None
            digit = name[4:] if name.startswith('num_') else name
            if digit.isdigit() and 1 <= int(digit) <= 9:
                return int(digit)
            if name in ('left', 'down'):
                state['sel'] = max(1, (state['sel'] or 1) - 1)
            elif name in ('right', 'up'):
                state['sel'] = min(9, (state['sel'] or 0) + 1)
            elif name in ('return', 'space') and state['sel']:
                return state['sel']


def countdown(win, ui_, kb, lang, seconds=5):
    from psychopy import core

    clock = core.Clock()
    while True:
        left = seconds - clock.getTime()
        if left <= 0:
            return True
        ui.draw_countdown(ui_, int(left) + 1, T(lang, 'countdown.hint'))
        win.flip()
        for k in kb.getKeys(waitRelease=False):
            if k.name in QUIT_KEYS:
                return False


def show_results(win, ui_, kb, rec):
    """The end-of-trial screen: how the trial went, and how it compares."""
    t = rec['scored']['totals']
    raw = rec['raw']
    cfg = rec['config']
    lang = rec.get('language', 'en')

    end_label = T(lang, {'max_duration': 'pp.endMax',
                         'sleep_onset': 'pp.endSleep',
                         'aborted': 'pp.endAborted'}.get(raw['endReason'], 'pp.endMax'))
    rows = [
        (T(lang, 'pp.outcome'), end_label),
        (T(lang, 'pp.stimuli'), '%d' % t['trials']),
        (T(lang, 'pp.hitsMisses'), '%d / %d' % (t['hits'], t['misses'])),
        (T(lang, 'pp.lapsesOver', ms=cfg['lapseMs']), '%d' % t['lapses']),
        (T(lang, 'pp.meanRt'), '—' if t['avgRt'] is None else '%.0f ms' % t['avgRt']),
        (T(lang, 'pp.framesDropped'),
         '%d (%d)' % (raw['frames'], raw['droppedFrames'])),
    ]
    if raw.get('sleepOnsetMs'):
        rows.insert(1, (T(lang, 'pp.sleepOnsetRow'), _mmss(raw['sleepOnsetMs'])))

    # A few headline norm comparisons, if one was available. The full table is
    # in the CSV; this is the experimenter's glance, not the analysis.
    norm_rows = []
    nm = rec.get('norms')
    if nm and nm.get('available'):
        by_key = {r['key']: r for r in nm['rows']}
        for key in ('rtMean', 'lapses', 'misses', 'rtIpr'):
            r = by_key.get(key)
            if not r or not r['comparison']:
                continue
            c = r['comparison']
            if c['z'] is None:
                value = T(lang, 'pp.noSpread')
            else:
                value = T(lang, 'pp.sdWorse' if c['z'] >= 0 else 'pp.sdBetter',
                          n='%.1f' % abs(c['z']))
            norm_rows.append((r['label'], value, c['band']))

    where = os.path.basename(os.path.dirname(rec['paths'][0]) or '.')

    def draw():
        ui.draw_results(ui_, T(lang, 'pp.finished'), rows, norm_rows,
                        T(lang, 'pp.savedTo', dir=where),
                        heading=T(lang, 'pp.normsHeading'))
    _wait_key(win, kb, draw)


def _mmss(ms):
    total = int(round(ms / 1000.0))
    return '%d:%02d' % (total // 60, total % 60)


def measure_frame_interval(win):
    """The real frame interval, measured, with the nominal rate as a fallback.

    Reported in the CSV either way, so an analysis can tell what the timing was
    actually quantised to rather than assuming 60 Hz.
    """
    try:
        hz = win.getActualFrameRate(nIdentical=20, nMaxFrames=240,
                                    nWarmUpFrames=20, threshold=1.0)
    except Exception:
        hz = None
    if hz and hz > 20:
        return 1000.0 / hz, hz
    period = getattr(win, 'monitorFramePeriod', None)
    if period:
        return period * 1000.0, 1.0 / period
    return 1000.0 / 60.0, 60.0


def run(settings, display_factory=None, rts=None):
    """Collect one trial and write its four CSVs. Returns the record."""
    cfg = cfg_from(settings)
    schedule = btask.build_schedule_for(cfg)
    lang = settings['language']
    now = datetime.datetime.now()

    simulated = display_factory is None and rts is not None
    win = kb = None
    alarm = None

    if simulated:
        display = SimulatedDisplay(rts)
        frame_ms, hz = display.frame, 1000.0 / display.frame
        kss_before = 4 if settings['kss_when'] in ('before', 'both') else ''
        kss_after = 7 if settings['kss_when'] in ('after', 'both') else ''
    else:
        from psychopy import visual, core                      # noqa: F401
        from psychopy.hardware import keyboard

        win = visual.Window(fullscr=bool(settings['fullscreen']), color=ui.BG,
                            colorSpace='rgb255', units='height', allowGUI=False,
                            waitBlanking=True)
        kb = keyboard.Keyboard()
        ui_ = ui.UI(win)
        frame_ms, hz = measure_frame_interval(win)

        s = STRINGS.get(lang, STRINGS['en'])
        if show_instructions(win, ui_, kb, lang, cfg['mode']) == 'escape':
            win.close()
            return None

        kss_before = ''
        kss_after = ''      # asked after the trial, but the record is built first
        if settings['kss_when'] in ('before', 'both'):
            kss_before = ask_kss(win, ui_, kb, lang, s['kss.beforeTitle']) or ''

        if not countdown(win, ui_, kb, lang):
            win.close()
            return None

        display = PsychoPyDisplay(win, kb, cfg['mode'], frame_ms)
        # The keyboard clock is deliberately NOT reset here. Key timestamps and
        # flip timestamps have to share an origin — a reaction time is one minus
        # the other — and PsychoPy puts both on the monotonic clock already.
        # Resetting one of them mid-run would silently offset every RT.
        kb.clearEvents()
        alarm = _make_alarm(cfg)

    alarm_failed = []

    def on_sleep_onset(index, t_ms):
        """Start the alarm — best effort, and never at the cost of the trial.

        This runs INSIDE the trial loop. An exception here propagates out of
        run_trial and out of run(), so before it was guarded a sound device
        that refused to play took the whole session with it: the participant
        reached sleep onset and no data was written at all. A missing alarm is
        an inconvenience; missing data is the experiment.
        """
        if alarm is None:
            return
        try:
            alarm.play(loops=-1)
        except Exception as e:
            alarm_failed.append('%s: %s' % (type(e).__name__, e))

    raw = btask.run_trial(display, schedule, cfg, on_sleep_onset=on_sleep_onset)
    raw['refreshHz'] = hz
    raw['frameMadMs'] = None
    scored = btask.score_trial(raw, cfg)

    norms = None
    try:
        norms = bnorms.Norms.load()
    except Exception:
        norms = None

    report = bnorms.normative_report({
        'mode': cfg['mode'], 'isiMs': cfg['isiMs'], 'hour': now.hour,
        'elapsedMs': raw['elapsedMs'], 'minuteBuckets': len(scored['perMinute']),
        'trials': scored['trials'], 'falseStartMs': cfg['falseStartMs'],
        'lapseMs': cfg['lapseMs'], 'missCriterion': cfg['missCriterion'],
    }, norms) if norms else None

    def build(kss_after_value):
        return {
            'runId': '%s-%s' % (now.strftime('%Y%m%dT%H%M%S'), os.getpid()),
            'language': lang,
            'participant': {
                'participantId': settings['participant_id'],
                'name': settings['name'], 'address': '',
                'birthDate': settings['birth_date'],
                'education': settings['education'],
                'sessionLabel': settings['session_label'],
                'trialNumber': settings['trial_number'],
                'date': now.strftime('%Y-%m-%d'), 'time': now.strftime('%H:%M:%S'),
            },
            'config': cfg, 'schedule': schedule, 'raw': raw, 'scored': scored,
            'norms': report,
            'kssWhen': settings['kss_when'],
            'kssBefore': kss_before, 'kssAfter': kss_after_value,
        }

    # SAVE BEFORE ANYTHING ELSE CAN GO WRONG.
    #
    # What remains after this point is screens: the sleep-onset banner, the
    # closing KSS, the results. Every one of them can fail on some machine, and
    # none of them is worth a participant's trial. The files are written here,
    # and rewritten below if the closing KSS adds an answer — the filenames are
    # the same, so the second write replaces the first rather than accumulating.
    rec = build(kss_after)
    paths = bio.write_all(rec, settings['out_dir'])
    rec['paths'] = paths

    if not simulated:
        try:
            if raw['endReason'] == 'sleep_onset':
                # The banner is shown whether or not the sound worked: it is the
                # visual half of the alert, and on a machine with no working
                # audio it is the only half there is.
                _wait_for_silence(win, ui_, kb, alarm, raw.get('sleepOnsetMs'),
                                  bool(alarm_failed), lang)
            if settings['kss_when'] in ('after', 'both'):
                answer = ask_kss(win, ui_, kb, lang,
                                 STRINGS.get(lang, STRINGS['en'])['kss.afterTitle'])
                if answer:
                    rec = build(answer)
                    rec['paths'] = bio.write_all(rec, settings['out_dir'])
        except Exception as e:
            _report_but_continue(e, paths)

        try:
            show_results(win, ui_, kb, rec)
        except Exception as e:
            _report_but_continue(e, rec['paths'])
        try:
            win.close()
        except Exception:
            pass
    return rec


def _report_but_continue(exc, paths):
    """Say what broke, and where the data already is.

    Printed rather than raised: by the time anything here runs the trial is
    over and saved, so the useful thing is to tell the experimenter both facts
    instead of ending with a traceback that looks like the session was lost.
    """
    import traceback
    sys.stderr.write('\n--- a screen after the trial failed ---\n')
    traceback.print_exc()
    sys.stderr.write('THE DATA IS SAFE. Already written:\n')
    for p in paths or []:
        sys.stderr.write('  %s\n' % p)
    sys.stderr.write('---------------------------------------\n')


def _make_alarm(cfg):
    if not cfg.get('alarm'):
        return None
    try:
        from psychopy import sound
        # A pulsing alert rather than a single beep; it loops until silenced.
        return sound.Sound(value=880, secs=0.4, stereo=True)
    except Exception:
        return None


def _wait_for_silence(win, ui_, kb, alarm, sleep_onset_ms, sound_failed=False,
                      lang='en'):
    """Hold the sleep-onset banner until the experimenter silences it.

    Shown whether or not the sound is working — on a machine with no usable
    audio device the banner is the whole alert, so it must not depend on the
    alarm having played.
    """
    when = _mmss(sleep_onset_ms) if sleep_onset_ms else '—'
    hint = T(lang, 'pp.silenceQuiet' if (alarm is None or sound_failed) else 'pp.silence')
    kb.clearEvents()
    while True:
        ui.draw_sleep_onset(ui_, T(lang, 'pp.sleepOnset'),
                            T(lang, 'pp.sleepOnsetAt', t=when), hint)
        win.flip()
        for k in kb.getKeys(waitRelease=False):
            if k.name in SILENCE_KEYS or k.name in QUIT_KEYS:
                if alarm is not None:
                    try:
                        alarm.stop()
                    except Exception:
                        pass
                return


def ask_settings():
    """Collect the settings, in two steps.

    The paradigm is asked first, on its own, because it changes what the rest of
    the defaults should be. A PVT with the sleep-onset criterion left on is not
    a PVT — it would stop at seven consecutive misses — and expecting the
    experimenter to notice and zero that box is exactly the kind of silent
    mis-configuration this instrument should not have. The browser build solves
    it by re-ticking the box when the mode changes; a modal dialog cannot, so it
    is asked in two passes instead.
    """
    from psychopy import gui

    first = {'mode': ['bsrt', 'pvt'], 'language': ['en', 'fr', 'nl', 'de']}
    dlg = gui.DlgFromDict(first, title='BSRT — paradigm',
                          order=['mode', 'language'],
                          tip={'mode': 'bsrt = flashing dot every 3 s; '
                                       'pvt = counter at 2-10 s intervals'})
    if not dlg.OK:
        return None

    mode = _first(first['mode'], 'bsrt')
    settings = dict(DEFAULTS)
    settings.update(MODE_DEFAULTS.get(mode, {}))
    settings['mode'] = mode
    settings['language'] = _first(first['language'], 'en')

    order = ['participant_id', 'session_label', 'trial_number', 'name',
             'birth_date', 'education', 'kss_when', 'max_minutes',
             'miss_criterion', 'stim_ms', 'hit_window_ms', 'lapse_ms',
             'alarm', 'fullscreen', 'remove_false_starts', 'remove_outliers',
             'out_dir', 'seed']
    # Only show the interval control that applies: a fixed ISI is meaningless
    # for a PVT, and an interval set is meaningless for a BSRT.
    if mode == 'pvt':
        settings.pop('isi_ms', None)
        order.insert(8, 'isi_set_s')
    else:
        settings.pop('isi_set_s', None)
        order.insert(8, 'isi_ms')
    settings.pop('mode', None)
    settings.pop('language', None)

    tips = {
        'seed': '0 picks a random seed and records it, so the schedule is reproducible',
        'kss_when': 'none, before, after or both',
        'birth_date': 'YYYY-MM-DD',
        'max_minutes': 'the ceiling; a BSRT can still stop earlier at sleep onset',
        'miss_criterion': 'consecutive misses scored as sleep onset; 0 turns it off',
        'isi_set_s': 'PVT intervals in seconds, e.g. 2/4/6/8/10',
        'out_dir': 'where the four CSVs are written',
    }
    dlg2 = gui.DlgFromDict(settings, title='BSRT — %s' % mode.upper(),
                           order=[k for k in order if k in settings], tip=tips)
    if not dlg2.OK:
        return None

    settings['mode'] = mode
    settings['language'] = _first(first['language'], 'en')
    settings.setdefault('isi_ms', DEFAULTS['isi_ms'])
    settings.setdefault('isi_set_s', DEFAULTS['isi_set_s'])
    return settings


def _first(value, fallback):
    """A DlgFromDict list field comes back as the chosen string, but an
    untouched one can still be the list. Take the first entry either way."""
    if isinstance(value, (list, tuple)):
        return value[0] if value else fallback
    return value or fallback


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--simulate', action='store_true',
                    help='run the whole pipeline with no screen and no PsychoPy')
    ap.add_argument('--minutes', type=int, default=3, help='with --simulate')
    ap.add_argument('--mode', default='bsrt', choices=['bsrt', 'pvt'])
    ap.add_argument('--out', default='data', help='output directory')
    args = ap.parse_args(argv)

    if args.simulate:
        settings = dict(DEFAULTS, participant_id='SIM001', mode=args.mode,
                        max_minutes=args.minutes, out_dir=args.out, seed=20260818,
                        kss_when='both', fullscreen=False)
        n = args.minutes * (20 if args.mode == 'bsrt' else 10) + 10
        r = random.Random(1)
        rts = [None if i % 23 == 22 else 250.0 + r.random() * 300.0 for i in range(n)]
        rec = run(settings, rts=rts)
        t = rec['scored']['totals']
        print('simulated %s, %d minutes' % (args.mode, args.minutes))
        print('  stimuli %d  hits %d  misses %d  lapses %d  mean RT %.1f ms'
              % (t['trials'], t['hits'], t['misses'], t['lapses'], t['avgRt']))
        print('  ended: %s   frames %d, dropped %d'
              % (rec['raw']['endReason'], rec['raw']['frames'], rec['raw']['droppedFrames']))
        for p in rec['paths']:
            print('  wrote', p)
        return 0

    settings = ask_settings()
    if settings is None:
        return 1
    rec = run(settings)
    return 0 if rec else 1


if __name__ == '__main__':
    sys.exit(main())
