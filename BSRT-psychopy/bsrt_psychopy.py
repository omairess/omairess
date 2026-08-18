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
        'isiSetMs': [2000, 4000, 6000, 8000, 10000],
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

        self.dot = visual.Circle(win, radius=0.03, units='height',
                                 fillColor='red', lineColor=None)
        self.counter = visual.TextStim(win, text='', height=0.12, color='white',
                                       units='height')
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


def show_text(win, kb, text, wait_keys=True, height=0.045):
    from psychopy import visual

    stim = visual.TextStim(win, text=text, height=height, color='white',
                           units='height', wrapWidth=1.4)
    if not wait_keys:
        stim.draw()
        win.flip()
        return None
    kb.clearEvents()
    # Redraw and flip every pass rather than spinning on getKeys: a bare polling
    # loop pegs a core and leaves the window unresponsive to the OS.
    while True:
        stim.draw()
        win.flip()
        for k in kb.getKeys(waitRelease=False):
            if k.name in QUIT_KEYS:
                return 'quit'
            return k.name


def ask_kss(win, kb, lang, title):
    """The 9-point Karolinska Sleepiness Scale, one screen, keys 1-9.

    Anchors come from strings.json, which is extracted from i18n.js — the same
    wording the browser build shows.
    """
    s = STRINGS.get(lang, STRINGS['en'])
    anchors = s['_anchors']
    lines = [title, '', s['kss.question'], '']
    for i, a in enumerate(anchors):
        lines.append('%d   %s' % (i + 1, a))
    lines += ['', s['kss.instruction']]
    text = '\n'.join(lines)

    from psychopy import visual
    stim = visual.TextStim(win, text=text, height=0.038, color='white',
                           units='height', wrapWidth=1.6)
    kb.clearEvents()
    while True:
        stim.draw()
        win.flip()
        for k in kb.getKeys(waitRelease=False):
            if k.name in QUIT_KEYS:
                return None
            if k.name and k.name.isdigit() and 1 <= int(k.name) <= 9:
                return int(k.name)
            if k.name and k.name.startswith('num_'):
                d = k.name[4:]
                if d.isdigit() and 1 <= int(d) <= 9:
                    return int(d)


def countdown(win, kb, lang, seconds=5):
    """A visible countdown, so nobody is surprised by the first stimulus."""
    from psychopy import visual

    s = STRINGS.get(lang, STRINGS['en'])
    stim = visual.TextStim(win, text='', height=0.15, color='white', units='height')
    hint = visual.TextStim(win, text=s['countdown.hint'], height=0.035, color='grey',
                           units='height', pos=(0, -0.3), wrapWidth=1.4)
    from psychopy import core
    clock = core.Clock()
    while True:
        left = seconds - clock.getTime()
        if left <= 0:
            break
        stim.text = '%d' % int(left + 1)
        stim.draw()
        hint.draw()
        win.flip()
        for k in kb.getKeys(waitRelease=False):
            if k.name in QUIT_KEYS:
                return False
    return True


# --------------------------------------------------------------------------


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

        win = visual.Window(fullscr=bool(settings['fullscreen']), color='black',
                            units='height', allowGUI=False, waitBlanking=True)
        kb = keyboard.Keyboard()
        frame_ms, hz = measure_frame_interval(win)

        s = STRINGS.get(lang, STRINGS['en'])
        hint = s['task.hintPvt'] if cfg['mode'] == 'pvt' else s['task.hintBsrt']
        if show_text(win, kb, hint + '\n\n\n[ press any key to continue ]') == 'quit':
            win.close()
            return None

        kss_before = ''
        if settings['kss_when'] in ('before', 'both'):
            kss_before = ask_kss(win, kb, lang, s['kss.beforeTitle']) or ''

        if not countdown(win, kb, lang):
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
            _wait_for_silence(win, kb, alarm)
        if settings['kss_when'] in ('after', 'both'):
            kss_after = ask_kss(win, kb, lang, STRINGS.get(lang, STRINGS['en'])['kss.afterTitle']) or ''
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
        show_text(win, kb, _end_text(rec, paths))
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


def _wait_for_silence(win, kb, alarm):
    from psychopy import visual

    msg = visual.TextStim(
        win, text='SLEEP ONSET\n\nthe criterion was reached\n\n[ press S to silence ]',
        height=0.06, color='red', units='height', wrapWidth=1.4)
    kb.clearEvents()
    while True:
        msg.draw()
        win.flip()
        for k in kb.getKeys(waitRelease=False):
            if k.name in SILENCE_KEYS or k.name in QUIT_KEYS:
                try:
                    alarm.stop()
                except Exception:
                    pass
                return


def _end_text(rec, paths):
    t = rec['scored']['totals']
    raw = rec['raw']
    lines = [
        'Trial finished.', '',
        'ended: %s' % raw['endReason'],
        'stimuli: %d    hits: %d    misses: %d    lapses: %d'
        % (t['trials'], t['hits'], t['misses'], t['lapses']),
        'mean RT: %s ms' % ('—' if t['avgRt'] is None else '%.0f' % t['avgRt']),
        'frames: %d    dropped: %d' % (raw['frames'], raw['droppedFrames']),
        '', 'saved:',
    ]
    lines += ['   ' + os.path.basename(p) for p in paths]
    lines += ['', '[ press any key to close ]']
    return '\n'.join(lines)


def ask_settings():
    """The setup dialog. Everything the browser build's setup screen asks."""
    from psychopy import gui

    settings = dict(DEFAULTS)
    dlg = gui.DlgFromDict(
        settings, title='BSRT',
        order=['participant_id', 'session_label', 'trial_number', 'name',
               'birth_date', 'education', 'language', 'mode', 'isi_ms', 'stim_ms',
               'hit_window_ms', 'lapse_ms', 'miss_criterion', 'max_minutes',
               'kss_when', 'alarm', 'fullscreen', 'remove_false_starts',
               'remove_outliers', 'out_dir', 'seed'],
        tip={'seed': '0 picks a random seed and records it, so the schedule is reproducible',
             'mode': 'bsrt or pvt',
             'kss_when': 'none, before, after or both',
             'birth_date': 'YYYY-MM-DD',
             'max_minutes': 'the ceiling; a BSRT can still stop earlier at sleep onset'})
    if not dlg.OK:
        return None
    return settings


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
