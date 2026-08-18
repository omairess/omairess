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

    print('\n=== the setup dialog ===')
    pgui.QUEUE = []
    check(app.ask_settings() is None, 'cancelling the first dialog returns nothing')

    # Two stages: the paradigm, then everything else.
    pgui.QUEUE = [{'mode': 'bsrt', 'language': 'en'}, {'participant_id': 'P9'}]
    got = app.ask_settings()
    check(got is not None and got['participant_id'] == 'P9',
          'an accepted dialog returns the answers')
    check(got['mode'] == 'bsrt' and 'max_minutes' in got, 'with the defaults filled in')
    check(got['miss_criterion'] == 7, 'a BSRT keeps the sleep-onset criterion (7)')

    # The bug this two-stage dialog exists to prevent: a PVT must NOT stop at
    # seven consecutive misses, matching the browser build.
    pgui.QUEUE = [{'mode': 'pvt', 'language': 'en'}, {'participant_id': 'P9'}]
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
    check(len(set(sched_pvt['isis'])) == 5,
          'using all five intervals: %s' % sorted(set(sched_pvt['isis'])))

    # A custom interval set typed into the box is honoured.
    pgui.QUEUE = [{'mode': 'pvt', 'language': 'en'},
                  {'participant_id': 'P9', 'isi_set_s': '1/2/3'}]
    got_custom = app.ask_settings()
    check(app.cfg_from(got_custom)['isiSetMs'] == [1000, 2000, 3000],
          'a custom PVT interval set is parsed: %s'
          % app.cfg_from(got_custom)['isiSetMs'])
finally:
    shutil.rmtree(out_dir, ignore_errors=True)

print()
if failures:
    print('%d PSYCHOPY-PATH TEST(S) FAILED' % len(failures))
    sys.exit(1)
print('THE PSYCHOPY CODE PATH RUNS END TO END')
