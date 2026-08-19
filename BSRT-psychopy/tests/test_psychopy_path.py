"""Run bsrt_psychopy.py's real PsychoPy branch against a stub.

This is the code that only executes on a machine with PsychoPy: window setup,
the instruction and KSS screens, the countdown, the trial itself, the
alarm-and-silence flow, and the export. Without this it would ship completely
unexecuted.

Read tests/stub_psychopy/psychopy/__init__.py for what this does NOT prove: the
stub is written from the same reading of the PsychoPy docs as the code it tests,
so it catches our mistakes but not a misread API.
"""

import csv
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'stub_psychopy'))     # stub before anything else
sys.path.insert(0, os.path.dirname(HERE))

import psychopy.core as pcore                                # noqa: E402
import psychopy.gui as pgui                                  # noqa: E402
from psychopy.hardware import keyboard as pkb                # noqa: E402

import psychopy as ppkg                                      # noqa: E402
import psychopy.sound as psound                              # noqa: E402
import psychopy.visual as pvisual                            # noqa: E402

import bsrt_psychopy as app                                  # noqa: E402


def _last_window():
    return pvisual.LAST[0]


def _last_sound():
    return psound.LAST[0]

failures = []


def check(cond, msg):
    print(('  ok    ' if cond else '  FAIL  ') + msg)
    if not cond:
        failures.append(msg)


FRAME = 1.0 / 60.0


def schedule_keys(events):
    pkb.Keyboard.SCHEDULE = sorted(events)


def base_settings(**kw):
    s = dict(app.DEFAULTS)
    s.update({'participant_id': 'STUB01', 'session_label': 'baseline',
              'trial_number': 1, 'max_minutes': 1, 'seed': 4242,
              'fullscreen': False, 'out_dir': None})
    s.update(kw)
    return s


def responses_for(schedule, rt_s=0.30, skip=()):
    """A press `rt_s` after each stimulus onset, on the virtual clock.

    Onsets are the intended ones plus the instruction/countdown time that runs
    before the trial, which the caller passes in as `offset`.
    """
    return [(o / 1000.0 + rt_s, 'space')
            for i, o in enumerate(schedule['onsets']) if i not in skip]


def run_once(settings, key_events, pre_seconds):
    """Reset the virtual world, queue answers and keys, run the program."""
    pcore.reset()
    pgui.QUEUE = [dict(settings)]
    # The trial does not start at t=0: the instruction screen, the KSS and the
    # countdown all consume virtual time first, so responses are offset by it.
    schedule_keys([(t + pre_seconds, k) for (t, k) in key_events])
    return app.run(settings)


out_dir = tempfile.mkdtemp()
try:
    print('=== a whole trial through the PsychoPy path ===')
    settings = base_settings(out_dir=out_dir, kss_when='none')
    cfg = app.cfg_from(settings)
    sched = app.btask.build_schedule_for(cfg)

    # One key to dismiss the instructions, then 5 s of countdown, then the trial.
    # The instruction screen ends on the first key at or after t=0.
    pre = FRAME + 5.0 + FRAME * 2
    events = [(FRAME * 0.5, 'space')] + [(t + pre, k) for (t, k) in responses_for(sched)]
    # Dismiss the end screen just after the trial finishes. Scheduling it far in
    # the future would make the virtual clock grind through every frame between.
    end_s = pre + sched['plannedDurationMs'] / 1000.0 + 1.0
    events.append((end_s, 'space'))
    pcore.reset()
    pgui.QUEUE = [dict(settings)]
    schedule_keys(events)
    rec = app.run(settings)

    check(rec is not None, 'the program ran and returned a record')
    check(len(rec['scored']['trials']) == 20,
          'a 1-minute BSRT delivered 20 stimuli (got %d)' % len(rec['scored']['trials']))
    hits = rec['scored']['totals']['hits']
    check(hits == 20, 'every response was caught as a hit (got %d)' % hits)
    rts = [t['rtMs'] for t in rec['scored']['trials'] if t['rtMs'] is not None]
    check(rts and all(abs(r - 300.0) < FRAME * 1000 + 1 for r in rts),
          'RTs land within a frame of the 300 ms pressed (mean %.1f ms)'
          % (sum(rts) / len(rts)))
    check(rec['raw']['endReason'] == 'max_duration', 'it ended at the ceiling')
    check(rec['raw']['droppedFrames'] == 0, 'no frames dropped on a virtual display')

    print('\n=== it wrote the four files ===')
    check(len(rec['paths']) == 4, 'four CSVs')
    for p in rec['paths']:
        with open(p, newline='') as f:
            rows = list(csv.reader(f))
        check(len(rows) >= 2 and len({len(r) for r in rows}) == 1,
              '%-38s %d rows, all the same width' % (os.path.basename(p), len(rows) - 1))

    print('\n=== the stimulus actually reached the screen ===')
    # Every frame the engine asked for a dot should have gone through the window.
    win_frames = [f for f in _last_window().frames if f]
    check(len(win_frames) > 0, 'frames were drawn (%d non-empty of %d)'
          % (len(win_frames), len(_last_window().frames)))

    print('\n=== the KSS is asked and recorded ===')
    settings2 = base_settings(out_dir=out_dir, kss_when='both')
    cfg2 = app.cfg_from(settings2)
    sched2 = app.btask.build_schedule_for(cfg2)
    # The instruction screen polls at t=FRAME and consumes every key ready by
    # then, so the KSS answer has to be scheduled strictly after that poll or it
    # is swallowed along with the key that dismissed the instructions.
    pre2 = FRAME * 3 + 5.0 + FRAME * 2         # instructions, KSS answer, countdown
    ev = [(FRAME * 0.5, 'space'), (FRAME * 1.5, '3')]
    ev += [(t + pre2, k) for (t, k) in responses_for(sched2)]
    end2 = pre2 + sched2['plannedDurationMs'] / 1000.0 + 1.0
    ev += [(end2, '8'), (end2 + 1.0, 'space')]   # KSS after, then the end screen
    pcore.reset()
    pgui.QUEUE = [dict(settings2)]
    schedule_keys(ev)
    rec2 = app.run(settings2)
    check(rec2['kssBefore'] == 3, 'the before answer was recorded (got %r)' % rec2['kssBefore'])
    check(rec2['kssAfter'] == 8, 'the after answer was recorded (got %r)' % rec2['kssAfter'])

    print('\n=== sleep onset sounds the alarm and waits to be silenced ===')
    settings3 = base_settings(out_dir=out_dir, kss_when='none', max_minutes=2,
                              miss_criterion=7)
    cfg3 = app.cfg_from(settings3)
    sched3 = app.btask.build_schedule_for(cfg3)
    pre3 = FRAME + 5.0 + FRAME * 2
    # Respond to the first 5, then stop entirely: seven consecutive misses.
    ev3 = [(FRAME * 0.5, 'space')]
    ev3 += [(o / 1000.0 + 0.30 + pre3, 'space') for o in sched3['onsets'][:5]]
    # Sleep onset fires when the twelfth epoch CLOSES, not when it starts. The
    # silence key has to come after that: during the trial the display swallows
    # 's' as an experimenter key, so an early one would simply vanish.
    end3 = pre3 + (sched3['onsets'][11] + sched3['epochIsi'][11]) / 1000.0
    ev3 += [(end3 + 1.0, 's'), (end3 + 2.0, 'space')]
    pcore.reset()
    pgui.QUEUE = [dict(settings3)]
    schedule_keys(ev3)
    rec3 = app.run(settings3)
    check(rec3['raw']['endReason'] == 'sleep_onset', 'the trial stopped at sleep onset')
    check(len(rec3['scored']['trials']) == 12,
          '5 hits then 7 misses = 12 epochs (got %d)' % len(rec3['scored']['trials']))
    check(rec3['raw']['sleepOnsetMs'] is not None, 'a sleep-onset time was recorded')
    snd = _last_sound()
    check(snd is not None and snd.plays == 1, 'the alarm played once')
    check(snd is not None and snd.stops == 1, 'and was stopped when S was pressed')
    check(snd is not None and snd.loops == -1, 'it was looping, not a single beep')

    print('\n=== escape aborts but still saves ===')
    settings4 = base_settings(out_dir=out_dir, kss_when='none', max_minutes=2)
    cfg4 = app.cfg_from(settings4)
    sched4 = app.btask.build_schedule_for(cfg4)
    pre4 = FRAME + 5.0 + FRAME * 2
    ev4 = [(FRAME * 0.5, 'space')]
    ev4 += [(o / 1000.0 + 0.30 + pre4, 'space') for o in sched4['onsets'][:4]]
    ev4 += [(pre4 + 20.0, 'escape'), (pre4 + 22.0, 'space')]
    pcore.reset()
    pgui.QUEUE = [dict(settings4)]
    schedule_keys(ev4)
    rec4 = app.run(settings4)
    check(rec4['raw']['endReason'] == 'aborted', 'the abort is recorded as such')
    check(len(rec4['paths']) == 4, 'and the data collected so far was still written')
    check(len(rec4['scored']['trials']) < sched4['nStimuli'],
          'with fewer epochs than the full schedule (%d of %d)'
          % (len(rec4['scored']['trials']), sched4['nStimuli']))

    print('\n=== a screen failing after the trial cannot cost the data ===')
    # The reported symptom: reach sleep onset, get no data folder. The cause was
    # alarm.play() raising inside the trial loop on a machine whose audio device
    # would not open, which took run_trial and the whole save with it.
    settings6 = base_settings(out_dir=out_dir, kss_when='none', max_minutes=2,
                              miss_criterion=7, participant_id='NOSOUND')
    cfg6 = app.cfg_from(settings6)
    sched6 = app.btask.build_schedule_for(cfg6)
    pre6 = FRAME + 5.0 + FRAME * 2
    ev6 = [(FRAME * 0.5, 'space')]
    ev6 += [(o / 1000.0 + 0.30 + pre6, 'space') for o in sched6['onsets'][:5]]
    end6 = pre6 + (sched6['onsets'][11] + sched6['epochIsi'][11]) / 1000.0
    ev6 += [(end6 + 1.0, 's'), (end6 + 2.0, 'space')]
    ppkg.FAIL_SOUND[0] = True
    try:
        pcore.reset()
        pgui.QUEUE = [dict(settings6)]
        schedule_keys(ev6)
        rec6 = app.run(settings6)
    finally:
        ppkg.FAIL_SOUND[0] = False
    check(rec6 is not None, 'the trial survived an audio device that will not play')
    check(rec6['raw']['endReason'] == 'sleep_onset', 'it still ended at sleep onset')
    written = [p for p in os.listdir(out_dir) if 'NOSOUND' in p]
    check(len(written) == 4, 'and all four files were written anyway (got %d)' % len(written))
    check(len(rec6['scored']['trials']) == 12, '12 epochs recorded')

    print('\n=== the data is on disk before the closing screens run ===')
    # Any screen after the trial is best-effort; the files must already exist.
    settings7 = base_settings(out_dir=out_dir, kss_when='none', max_minutes=1,
                              participant_id='CRASHER')
    cfg7 = app.cfg_from(settings7)
    sched7 = app.btask.build_schedule_for(cfg7)
    pre7 = FRAME + 5.0 + FRAME * 2
    ev7 = [(FRAME * 0.5, 'space')]
    ev7 += [(o / 1000.0 + 0.30 + pre7, 'space') for o in sched7['onsets']]
    pcore.reset()
    pgui.QUEUE = [dict(settings7)]
    schedule_keys(ev7)          # no key to dismiss the results screen
    real_results = app.show_results
    app.show_results = lambda *a, **k: (_ for _ in ()).throw(RuntimeError('screen blew up'))
    try:
        rec7 = app.run(settings7)
    finally:
        app.show_results = real_results
    check(rec7 is not None, 'a crashing results screen does not end the run')
    got7 = [p for p in os.listdir(out_dir) if 'CRASHER' in p]
    check(len(got7) == 4, 'the four files are already on disk (got %d)' % len(got7))

    print('\n=== a PVT runs through the same path ===')
    settings5 = base_settings(out_dir=out_dir, kss_when='none', mode='pvt',
                              max_minutes=2, miss_criterion=0)
    cfg5 = app.cfg_from(settings5)
    sched5 = app.btask.build_schedule_for(cfg5)
    pre5 = FRAME + 5.0 + FRAME * 2
    ev5 = [(FRAME * 0.5, 'space')]
    ev5 += [(o / 1000.0 + 0.28 + pre5, 'space') for o in sched5['onsets']]
    ev5.append((pre5 + sched5['plannedDurationMs'] / 1000.0 + 1.0, 'space'))
    pcore.reset()
    pgui.QUEUE = [dict(settings5)]
    schedule_keys(ev5)
    rec5 = app.run(settings5)
    check(rec5 is not None, 'the PVT ran')
    check(len(rec5['scored']['trials']) == 20,
          'a 2-minute PVT delivered 20 stimuli (got %d)' % len(rec5['scored']['trials']))
    check(rec5['schedule']['method'] == 'block_permutation',
          'on a block-permuted schedule')
    isis = sorted(set(rec5['schedule']['isis']))
    check(isis == [2000, 4000, 6000, 8000, 10000],
          'with genuinely variable intervals: %s' % isis)
    # Each 30 s block must hold exactly one of each interval — that is what
    # 'balanced' means, and it is the property a naive random draw would lose.
    blocks = {}
    for b, v in zip(rec5['schedule']['blocks'], rec5['schedule']['isis']):
        blocks.setdefault(b, []).append(v)
    balanced = all(sorted(v) == isis for v in blocks.values() if len(v) == 5)
    check(balanced, 'and every full block holds one of each (%d blocks)' % len(blocks))
    check(rec5['raw']['endReason'] == 'max_duration',
          'it ran to the ceiling rather than stopping at sleep onset')
    check(len(rec5['paths']) == 4, 'and wrote its four CSVs')
    import csv as _csv
    with open([p for p in rec5['paths'] if p.endswith('_summary.csv')][0], newline='') as f:
        srow = list(_csv.DictReader(f))[0]
    check(srow['mode'] == 'pvt', 'the summary records mode=pvt')
    check(srow['isi_set_s'] == '2.0/4.0/6.0/8.0/10.0',
          'and the interval set used: %s' % srow['isi_set_s'])
    check(srow['schedule_method'] == 'block_permutation', 'and the schedule rule')

    print('\n=== the screens are translated ===')
    # The symptom was a half-translated interface: the strings that came from
    # i18n.js changed language, the ones hardcoded here did not. This runs a
    # whole trial in French and reads back every piece of text that was drawn.
    def texts_from(win):
        out = []
        for frame in win.frames:
            for stim in frame:
                t = getattr(stim, 'text', None)
                if isinstance(t, str) and t.strip():
                    out.append(t)
        return out

    settings8 = base_settings(out_dir=out_dir, kss_when='both', max_minutes=1,
                              language='fr', participant_id='FR01')
    cfg8 = app.cfg_from(settings8)
    sched8 = app.btask.build_schedule_for(cfg8)
    pre8 = FRAME * 3 + 5.0 + FRAME * 2
    ev8 = [(FRAME * 0.5, 'space'), (FRAME * 1.5, '4')]
    ev8 += [(o / 1000.0 + 0.30 + pre8, 'space') for o in sched8['onsets']]
    end8 = pre8 + sched8['plannedDurationMs'] / 1000.0 + 1.0
    ev8 += [(end8, '7'), (end8 + 1.0, 'space')]
    pcore.reset()
    pgui.QUEUE = [dict(settings8)]
    schedule_keys(ev8)
    rec8 = app.run(settings8)
    drawn = ' | '.join(texts_from(_last_window()))

    check('Test de résistance au sommeil' in drawn,
          'the instruction title is in French')
    check('appuyez sur une touche' in drawn, 'and its key hint')
    check('Extrêmement alerte' in drawn, 'the KSS anchors are French')
    check('appuyez sur 1' in drawn, 'and the KSS key hint')
    check('Essai terminé' in drawn, 'the results screen title is French')
    check('temps de réaction moyen' in drawn, 'and its row labels')
    check('appuyez sur une touche pour fermer' in drawn, 'and its footer')

    # 'stimuli' is deliberately not probed: it is the same word in French.
    english = [phrase for phrase in
               ('press any key', 'Trial finished', 'mean reaction time',
                'Sleep Resistance Task', 'press 1-9', 'press 1–9',
                'against the preliminary norms', 'hits / misses',
                'ran to the ceiling')
               if phrase in drawn]
    check(not english, 'no English text leaked through: %s' % (english or 'none'))

    # And the other two languages reach the screen at all.
    for lang, needle in (('nl', 'Slaapweerstandstaak'), ('de', 'Schlafresistenztest')):
        st = base_settings(out_dir=out_dir, kss_when='none', max_minutes=1,
                           language=lang, participant_id=lang.upper() + '01')
        c = app.cfg_from(st)
        sc = app.btask.build_schedule_for(c)
        pre = FRAME + 5.0 + FRAME * 2
        evs = [(FRAME * 0.5, 'space')]
        evs += [(o / 1000.0 + 0.30 + pre, 'space') for o in sc['onsets']]
        evs.append((pre + sc['plannedDurationMs'] / 1000.0 + 1.0, 'space'))
        pcore.reset()
        pgui.QUEUE = [dict(st)]
        schedule_keys(evs)
        app.run(st)
        check(needle in ' | '.join(texts_from(_last_window())),
              '%s renders in %s' % (needle, lang))

    print('\n=== the KSS shows an anchor on every step ===')
    # Labelling only the two ends turns the KSS into a bare 1-9 rating, which
    # is a different instrument. The browser build shows all nine, so this must.
    import bsrt_ui as ui_mod
    probe_win = pvisual.Window()
    probe = ui_mod.UI(probe_win)
    for lang in ('en', 'fr', 'nl', 'de'):
        anchors = app.STRINGS[lang]['_anchors']
        ui_mod.draw_kss(probe, 'KSS', app.T(lang, 'kss.question'), anchors, 5,
                        app.T(lang, 'pp.kssKeys'))
        probe_win.flip()
        shown = [st.text for st in probe_win.frames[-1] if st.text]
        missing = [a for a in anchors if a not in shown]
        check(not missing, 'all nine %s anchors are on screen%s'
              % (lang, '' if not missing else ' — missing %s' % missing))
        numbers = [t for t in shown if t.isdigit()]
        check(sorted(numbers, key=int) == [str(i) for i in range(1, 10)],
              'and every step is numbered 1-9')
    probe_win.close()

    print('\n=== the setup dialog ===')

    def queue(mode, lang, second):
        """Queue the two dialogs the way the user sees them: the first in
        English, the second labelled in the chosen language."""
        first = {app.T('en', 'ppf.mode'): mode, app.T('en', 'ppf.language'): lang}
        labels = app._labels(lang, second.keys())
        pgui.QUEUE = [first, dict((labels[k], v) for k, v in second.items())]

    pgui.QUEUE = []
    check(app.ask_settings() is None, 'cancelling the first dialog returns nothing')

    queue('bsrt', 'en', {'participant_id': 'P9'})
    got = app.ask_settings()
    check(got is not None and got['participant_id'] == 'P9',
          'an accepted dialog returns the answers by their internal names')
    check(got['mode'] == 'bsrt' and 'max_minutes' in got, 'with the defaults filled in')
    check(got['miss_criterion'] == 7, 'a BSRT keeps the sleep-onset criterion (7)')
    check(got['kss_when'] == 'none',
          'kss_when comes back as a plain value, not a list (%r)' % got['kss_when'])

    # The bug this two-stage dialog exists to prevent: a PVT must NOT stop at
    # seven consecutive misses, matching the browser build.
    queue('pvt', 'en', {'participant_id': 'P9'})
    got_pvt = app.ask_settings()
    check(got_pvt['mode'] == 'pvt', 'choosing pvt gives a PVT')
    check(got_pvt['miss_criterion'] == 0,
          'and the sleep-onset criterion is OFF (got %r)' % got_pvt['miss_criterion'])
    cfg_pvt = app.cfg_from(got_pvt)
    check(cfg_pvt['isiSetMs'] == [2000, 4000, 6000, 8000, 10000],
          'with the standard PVT interval set: %s' % cfg_pvt['isiSetMs'])
    sched_pvt = app.btask.build_schedule_for(dict(cfg_pvt, maxMinutes=3))
    check(sched_pvt['method'] == 'block_permutation',
          'and a block-permuted variable schedule')
    check(sched_pvt['nStimuli'] == 30,
          'a 3-minute PVT presents 30 stimuli (got %d)' % sched_pvt['nStimuli'])

    queue('pvt', 'en', {'participant_id': 'P9', 'isi_set_s': '1/2/3'})
    got_custom = app.ask_settings()
    check(app.cfg_from(got_custom)['isiSetMs'] == [1000, 2000, 3000],
          'a custom PVT interval set is parsed: %s'
          % app.cfg_from(got_custom)['isiSetMs'])

    print('\n=== the dialog is translated too ===')
    queue('bsrt', 'fr', {'participant_id': 'FR9'})
    got_fr = app.ask_settings()
    check(got_fr is not None and got_fr['participant_id'] == 'FR9',
          'a French dialog still answers by internal name')
    check(got_fr['language'] == 'fr', 'and records the language')
    check(app.T('fr', 'ppf.participant_id') == 'identifiant du participant',
          'the labels themselves are French')

    # Labels become dictionary keys, so a collision would silently drop a field.
    for lang in ('en', 'fr', 'nl', 'de'):
        labels = list(app._labels(lang, app.DEFAULTS.keys()).values())
        check(len(labels) == len(set(labels)),
              'every %s dialog label is unique (%d fields)' % (lang, len(labels)))

    # Values stay English so the exported columns do not depend on the
    # experimenter's interface language.
    check(app.CHOICES['kss_when'] == ['none', 'before', 'after', 'both'],
          'kss_when is a dropdown with untranslated values')
    check(app.CHOICES['mode'] == ['bsrt', 'pvt'], 'and so is the task')

finally:
    shutil.rmtree(out_dir, ignore_errors=True)

print()
if failures:
    print('%d PSYCHOPY-PATH TEST(S) FAILED' % len(failures))
    sys.exit(1)
print('THE PSYCHOPY CODE PATH RUNS END TO END')
