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
    s = STRINGS.get(lang, STRINGS['en'])
    hint = s['task.hintPvt'] if mode == 'pvt' else s['task.hintBsrt']
    title = 'Psychomotor Vigilance Task' if mode == 'pvt' else 'Sleep Resistance Task'

    def draw():
        ui.draw_instructions(ui_, title, [hint], 'press any key to continue')
    return _wait_key(win, kb, draw)


def ask_kss(win, ui_, kb, lang, title):
    """The 9-point Karolinska Sleepiness Scale as a scale, not a list.

    Anchors come from strings.json, extracted from i18n.js, so the wording is
    the browser build's wording. Highlighting follows the number under the
    cursor keys as well as the digits, because a row of boxes invites arrows.
    """
    s = STRINGS.get(lang, STRINGS['en'])
    anchors = s['_anchors']
    state = {'sel': 0}

    def draw():
        ui.draw_kss(ui_, title, s['kss.question'], anchors, state['sel'],
                    'press 1-9, or use the arrow keys and Enter')

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

    s = STRINGS.get(lang, STRINGS['en'])
    clock = core.Clock()
    while True:
        left = seconds - clock.getTime()
        if left <= 0:
            return True
        ui.draw_countdown(ui_, int(left) + 1, s['countdown.hint'])
        win.flip()
        for k in kb.getKeys(waitRelease=False):
            if k.name in QUIT_KEYS:
                return False


def show_results(win, ui_, kb, rec):
    """The end-of-trial screen: how the trial went, and how it compares."""
    t = rec['scored']['totals']
    raw = rec['raw']
    cfg = rec['config']

    end_label = {'max_duration': 'ran to the ceiling',
                 'sleep_onset': 'stopped at sleep onset',
                 'aborted': 'stopped by the experimenter'}.get(raw['endReason'],
                                                               raw['endReason'])
    rows = [
        ('outcome', end_label),
        ('stimuli', '%d' % t['trials']),
        ('hits / misses', '%d / %d' % (t['hits'], t['misses'])),
        ('lapses (> %d ms)' % cfg['lapseMs'], '%d' % t['lapses']),
        ('mean reaction time', '—' if t['avgRt'] is None else '%.0f ms' % t['avgRt']),
        ('frames (dropped)', '%d (%d)' % (raw['frames'], raw['droppedFrames'])),
    ]
    if raw.get('sleepOnsetMs'):
        rows.insert(1, ('sleep onset at', _mmss(raw['sleepOnsetMs'])))

    # A few headline norm comparisons, if one was available. The full table is
    # in the CSV; this is the experimenter's glance, not the analysis.
    norm_rows = []
    nm = rec.get('norms')
    if nm and nm.get('available'):
        wanted = ('rtMean', 'lapses', 'misses', 'rtIpr')
        by_key = {r['key']: r for r in nm['rows']}
        for key in wanted:
            r = by_key.get(key)
            if not r or not r['comparison']:
                continue
            c = r['comparison']
            if c['z'] is None:
                value = 'no spread in the reference'
            else:
                value = '%.1f SD %s' % (abs(c['z']), 'worse' if c['z'] >= 0 else 'better')
            norm_rows.append((r['label'], value, c['band']))

    def draw():
        ui.draw_results(ui_, 'Trial finished', rows, norm_rows,
                        'saved to %s   —   press any key to close'
                        % os.path.basename(os.path.dirname(rec['paths'][0]) or '.'))
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

    fired = {'at': None}

    def on_sleep_onset(index, t_ms):
        fired['at'] = t_ms
        if alarm is not None:
            alarm.play(loops=-1)

    raw = btask.run_trial(display, schedule, cfg, on_sleep_onset=on_sleep_onset)
    raw['refreshHz'] = hz
    raw['frameMadMs'] = None
    scored = btask.score_trial(raw, cfg)

    if not simulated:
        # The alarm keeps sounding until the experimenter silences it, so a
        # sleep-onset event cannot be missed by someone who stepped away.
        if raw['endReason'] == 'sleep_onset' and alarm is not None:
            _wait_for_silence(win, ui_, kb, alarm, raw.get('sleepOnsetMs'))
        if settings['kss_when'] in ('after', 'both'):
            kss_after = ask_kss(win, ui_, kb, lang,
                                STRINGS.get(lang, STRINGS['en'])['kss.afterTitle']) or ''
        else:
            kss_after = ''

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

    rec = {
        'runId': '%s-%s' % (now.strftime('%Y%m%dT%H%M%S'), os.getpid()),
        'language': lang,
        'participant': {
            'participantId': settings['participant_id'],
            'name': settings['name'], 'address': '',
            'birthDate': settings['birth_date'], 'education': settings['education'],
            'sessionLabel': settings['session_label'],
            'trialNumber': settings['trial_number'],
            'date': now.strftime('%Y-%m-%d'), 'time': now.strftime('%H:%M:%S'),
        },
        'config': cfg, 'schedule': schedule, 'raw': raw, 'scored': scored,
        'norms': report,
        'kssWhen': settings['kss_when'], 'kssBefore': kss_before, 'kssAfter': kss_after,
    }

    paths = bio.write_all(rec, settings['out_dir'])
    rec['paths'] = paths

    if not simulated:
        show_results(win, ui_, kb, rec)
        win.close()
    return rec


def _make_alarm(cfg):
    if not cfg.get('alarm'):
        return None
    try:
        from psychopy import sound
        # A pulsing alert rather than a single beep; it loops until silenced.
        return sound.Sound(value=880, secs=0.4, stereo=True)
    except Exception:
        return None


def _wait_for_silence(win, ui_, kb, alarm, sleep_onset_ms):
    """Hold the alarm until the experimenter silences it.

    The alarm keeps sounding until someone presses S, so a sleep-onset event
    cannot be missed by an experimenter who stepped away — the same behaviour
    as the browser build's banner.
    """
    when = _mmss(sleep_onset_ms) if sleep_onset_ms else '—'
    kb.clearEvents()
    while True:
        ui.draw_sleep_onset(ui_, when, 'press S to silence the alarm')
        win.flip()
        for k in kb.getKeys(waitRelease=False):
            if k.name in SILENCE_KEYS or k.name in QUIT_KEYS:
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
