"""The trial engine, driven by a fake display on a virtual clock.

No screen, no PsychoPy: the engine takes a display object, so the whole task can
be run deterministically and its timing asserted exactly. What this checks is
the part that is new in the PsychoPy build — the frame loop, epoch boundaries,
response capture and the sleep-onset stop. Scoring is checked separately, and
against the JavaScript, in test_equivalence.py.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import bsrt_core as core      # noqa: E402
import bsrt_task as task      # noqa: E402

failures = []


def check(cond, msg):
    if cond:
        print('  ok    ' + msg)
    else:
        print('  FAIL  ' + msg)
        failures.append(msg)


class FakeDisplay:
    """A display on a virtual clock.

    `rts` gives the reaction time to produce for each epoch in turn; None means
    no response at all. Presses are stamped at exactly onset + rt, the way a
    hardware-timestamped keyboard reports a press that happened between two
    flips, so the recorded RT can be asserted to the millisecond.
    """

    def __init__(self, rts=(), frame_ms=1000.0 / 60.0, drop_frames=(), abort_at=None,
                 extras=()):
        self.t = 1000.0                 # start away from zero to catch t0 bugs
        self.frame = frame_ms
        self.rts = list(rts)
        self.states = []
        self.flips = 0
        self.drop_frames = set(drop_frames)
        # extras: (epoch index, delay after onset) for presses beyond the first.
        self.extras = list(extras)
        self.abort_at = abort_at
        self._epoch = -1
        self._onset = None
        self._pending = []
        self._fired = False

    def frame_interval_ms(self):
        return self.frame

    def flip(self, state):
        # A dropped frame is one the display missed: the next flip lands two
        # (or more) intervals later instead of one.
        step = self.frame * (2 if self.flips in self.drop_frames else 1)
        self.t += step
        self.flips += 1

        showing = state[0] in ('dot', 'counter')
        if showing and self._onset is None:
            self._epoch += 1
            self._onset = self.t
            self._fired = False
            rt = self.rts[self._epoch] if self._epoch < len(self.rts) else None
            if rt is not None:
                self._pending.append(self._onset + rt)
            for ep, delay in self.extras:
                if ep == self._epoch:
                    self._pending.append(self._onset + delay)
        elif not showing and state[0] != 'feedback':
            self._onset = None

        self.states.append(state[0])
        return self.t

    def poll(self):
        out = [(t, 'space') for t in self._pending if t <= self.t]
        self._pending = [t for t in self._pending if t > self.t]
        return out

    def quit_requested(self):
        return self.abort_at is not None and self.flips >= self.abort_at


CFG = {'mode': 'bsrt', 'isiMs': 3000, 'stimMs': 1000, 'hitWindowMs': 1000,
       'lapseMs': 500, 'missCriterion': 7, 'maxMinutes': 3, 'seed': 42,
       'removeFalseStarts': False, 'removeOutliers': False}


def run(cfg, rts, **kw):
    sched = task.build_schedule_for(cfg)
    d = FakeDisplay(rts=rts, **kw)
    raw = task.run_trial(d, sched, cfg)
    return sched, d, raw


print('=== a 3-minute BSRT presents 60 stimuli ===')
sched, d, raw = run(CFG, [300.0] * 60)
check(sched['nStimuli'] == 60, 'the schedule holds 60 stimuli')
check(len(raw['epochs']) == 60, 'and 60 epochs were run (got %d)' % len(raw['epochs']))
check(raw['endReason'] == 'max_duration', 'it ended at the ceiling')
check(all(e['rtMs'] is not None for e in raw['epochs']), 'every epoch caught its response')

print('\n=== reaction times are exact ===')
want = [250.0, 480.0, 999.0, 1200.0, 120.0, 45.0]
sched, d, raw = run(dict(CFG, maxMinutes=1), want + [300.0] * 20)
got = [e['rtMs'] for e in raw['epochs'][:len(want)]]
check(all(abs(a - b) < 1e-9 for a, b in zip(got, want)),
      'recorded RTs match what was pressed: %s' % [round(g, 1) for g in got])

scored = task.score_trial(raw, dict(CFG, maxMinutes=1))
outcomes = [t['outcome'] for t in scored['trials'][:len(want)]]
check(outcomes == ['hit', 'hit', 'hit', 'miss', 'hit', 'hit'],
      'a 1200 ms response is a miss, 999 ms is still a hit: %s' % outcomes)
check(scored['trials'][3]['lateResponse'] == 1, 'and it is marked as a late response')
check(scored['trials'][1]['lapse'] == 0, '480 ms is under the lapse threshold')
check(scored['trials'][5]['falseStart'] == 1, '45 ms is a false start')

print('\n=== onsets follow the intended schedule ===')
sched, d, raw = run(CFG, [300.0] * 60)
errs = [e['onsetErrorMs'] for e in raw['epochs']]
worst = max(abs(e) for e in errs)
check(worst <= d.frame + 1e-6,
      'no stimulus is more than one frame late (worst %.2f ms, frame %.2f ms)' % (worst, d.frame))
check(all(e >= -1e-9 for e in errs), 'and none is early — a frame is never anticipated')
drift = raw['epochs'][-1]['onsetMs'] - raw['epochs'][-1]['intendedOnsetMs']
check(abs(drift) <= d.frame + 1e-6,
      'the last stimulus has not drifted after 60 epochs (%.2f ms)' % drift)

print('\n=== the sleep-onset criterion stops the trial ===')
rts = [300.0] * 10 + [None] * 7 + [300.0] * 40
fired = []
sched2 = task.build_schedule_for(CFG)
d2 = FakeDisplay(rts=rts)
raw2 = task.run_trial(d2, sched2, CFG, on_sleep_onset=lambda i, t: fired.append((i, t)))
check(raw2['endReason'] == 'sleep_onset', 'it ended on the criterion')
check(len(raw2['epochs']) == 17, 'after 17 epochs — 10 hits then 7 misses (got %d)'
      % len(raw2['epochs']))
check(len(fired) == 1, 'the alarm callback fired exactly once')
check(raw2['sleptBeforeMax'] is True, 'and the trial is flagged as sleep onset')
check(raw2['sleepOnsetMs'] is not None and raw2['sleepOnsetMs'] > 0,
      'with a sleep-onset time of %.1f s' % (raw2['sleepOnsetMs'] / 1000.0))

print('\n=== six misses do not stop it ===')
rts = [300.0] * 5 + [None] * 6 + [300.0] * 60
sched3, d3, raw3 = run(CFG, rts)
check(raw3['endReason'] == 'max_duration', 'a run of six is under the criterion of seven')
check(len(raw3['epochs']) == 60, 'so the trial ran to the end')
ep = core.error_profiles([t['outcome'] for t in task.score_trial(raw3, CFG)['trials']], 7)
check(ep['ep3_6'] == 1 and ep['ep7plus'] == 0, 'and it is scored as one EP3-6 run')

print('\n=== dropped frames are detected, not silently absorbed ===')
sched4 = task.build_schedule_for(dict(CFG, maxMinutes=1))
d4 = FakeDisplay(rts=[300.0] * 20, drop_frames=(10, 50, 90, 130))
raw4 = task.run_trial(d4, sched4, dict(CFG, maxMinutes=1))
check(raw4['droppedFrames'] == 4,
      'all four injected drops were counted (got %d)' % raw4['droppedFrames'])
check(raw4['frames'] > 0, 'and the frame total is reported (%d frames)' % raw4['frames'])
check(len(raw4['epochs']) == 20, 'the trial still delivered all 20 stimuli')

print('\n=== PVT mode ===')
pvt = dict(CFG, mode='pvt', maxMinutes=3,
           isiSetMs=[2000, 4000, 6000, 8000, 10000], blockMs=30000)
sched5, d5, raw5 = run(pvt, [280.0] * 40)
check(sched5['nStimuli'] == 30, 'a 3-minute PVT presents 30 stimuli (got %d)'
      % sched5['nStimuli'])
check(sched5['method'] == 'block_permutation', 'built by block permutation')
check(len(raw5['epochs']) == 30, 'and 30 epochs ran')
check('counter' in d5.states, 'the PVT stimulus is the counter, not the dot')
check('dot' not in d5.states, 'and the dot never appears in PVT mode')
sched6, d6, raw6 = run(CFG, [300.0] * 60)
check('dot' in d6.states and 'counter' not in d6.states, 'BSRT mode keeps the dot')

print('\n=== extra presses are recorded but never overwrite the response ===')
cfg7 = dict(CFG, maxMinutes=1)
sched7 = task.build_schedule_for(cfg7)
# Two more presses in epoch 0, at 700 ms and 1400 ms after its onset. The
# second lands after the stimulus is gone but while the epoch is still open.
d7 = FakeDisplay(rts=[300.0] * 20, extras=[(0, 700.0), (0, 1400.0)])
raw7 = task.run_trial(d7, sched7, cfg7)
check(abs(raw7['epochs'][0]['rtMs'] - 300.0) < 1e-6,
      'the FIRST press is the epoch response, not a later one (%.3f ms)'
      % raw7['epochs'][0]['rtMs'])
check(len(raw7['presses']) == 22, 'all 22 presses are kept for the integrity check (got %d)'
      % len(raw7['presses']))
in_ep0 = [p for p in raw7['presses'] if p['epochIndex'] == 0]
check(len(in_ep0) == 3, 'three of them are attributed to epoch 0 (got %d)' % len(in_ep0))
g7 = task.score_trial(raw7, cfg7)
check(g7['integrity']['extraPresses'] == 2, 'and two count as extra presses (got %d)'
      % g7['integrity']['extraPresses'])
check(g7['totals']['hits'] == 20, 'the extras did not create or destroy any hits')

print('\n=== no responses at all ===')
sched8, d8, raw8 = run(dict(CFG, missCriterion=0, maxMinutes=1), [None] * 20)
check(len(raw8['epochs']) == 20, 'the trial runs to the end when the criterion is off')
check(all(e['rtMs'] is None for e in raw8['epochs']), 'and every epoch is a miss')
s8 = task.score_trial(raw8, dict(CFG, missCriterion=0, maxMinutes=1))
check(s8['totals']['hits'] == 0 and s8['totals']['misses'] == 20, 'scored as 20 misses')
check(s8['totals']['avgRt'] is None, 'with no RT average rather than a zero')

print('\n=== the experimenter can abort ===')
sched9 = task.build_schedule_for(CFG)
d9 = FakeDisplay(rts=[300.0] * 60, abort_at=200)
raw9 = task.run_trial(d9, sched9, CFG)
check(raw9['endReason'] == 'aborted', 'the abort is recorded as such')
check(len(raw9['epochs']) < 60, 'and the trial stopped early (%d epochs)' % len(raw9['epochs']))

print()
if failures:
    print('%d TASK TEST(S) FAILED' % len(failures))
    sys.exit(1)
print('ALL TRIAL-ENGINE TESTS PASSED')
