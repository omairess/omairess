"""A whole trial, from schedule to four written CSVs, with no screen involved.

Matching header names is not enough — a row that is one value short still
writes, and only breaks when someone reads it. This runs the simulated task end
to end and checks every row of every file against its header.
"""

import csv
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import bsrt_io as io_        # noqa: E402
import bsrt_norms as bn      # noqa: E402
import bsrt_task as task     # noqa: E402

failures = []


def check(cond, msg):
    print(('  ok    ' if cond else '  FAIL  ') + msg)
    if not cond:
        failures.append(msg)


class FakeDisplay:
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


CFG = {'mode': 'bsrt', 'isiMs': 3000, 'stimMs': 1000, 'hitWindowMs': 1000,
       'lapseMs': 500, 'missCriterion': 7, 'maxMinutes': 3, 'seed': 20260818,
       'removeFalseStarts': True, 'removeOutliers': True, 'sdMultiplier': 2,
       'falseStartMs': 100, 'alarm': True, 'fullscreen': True,
       'isiSetMs': [2000, 4000, 6000, 8000, 10000], 'blockMs': 30000}

# A plausible trial: mostly hits, a few misses, one false start, one late one.
rts = []
for i in range(60):
    if i in (17, 18, 33):
        rts.append(None)
    elif i == 25:
        rts.append(60.0)            # false start
    elif i == 41:
        rts.append(1400.0)          # responded, too late
    elif i % 7 == 0:
        rts.append(620.0)           # a lapse
    else:
        rts.append(280.0 + (i % 11) * 9)
sched = task.build_schedule_for(CFG)
disp = FakeDisplay(rts)
raw = task.run_trial(disp, sched, CFG)
raw['refreshHz'] = 1000.0 / disp.frame
raw['frameMadMs'] = 0.0
scored = task.score_trial(raw, CFG)

norms = bn.Norms.load()
report = bn.normative_report({
    'mode': 'bsrt', 'isiMs': 3000, 'hour': 3, 'elapsedMs': raw['elapsedMs'],
    'minuteBuckets': len(scored['perMinute']), 'trials': scored['trials'],
    'falseStartMs': 100, 'lapseMs': 500, 'missCriterion': 7,
}, norms)

rec = {
    'runId': 'test-run-1', 'language': 'en',
    'participant': {'participantId': 'P001', 'name': 'Alex', 'address': '',
                    'birthDate': '1985-03-07', 'education': 'MSc',
                    'sessionLabel': 'baseline', 'trialNumber': 1,
                    'date': '2026-08-18', 'time': '03:14:00'},
    'config': CFG, 'schedule': sched, 'raw': raw, 'scored': scored, 'norms': report,
    'kssWhen': 'both', 'kssBefore': 4, 'kssAfter': 7,
}

print('=== the trial ran ===')
check(len(scored['trials']) == 60, '60 epochs scored')
check(scored['totals']['misses'] == 4, '4 misses — 3 no-responses and 1 late (got %d)'
      % scored['totals']['misses'])
check(scored['totals']['falseStarts'] == 1, '1 false start')
check(scored['totals']['lateResponses'] == 1, '1 late response')
check(len(scored['perMinute']) == 3, '3 minutes')

print('\n=== every row matches its header ===')
out = tempfile.mkdtemp()
try:
    paths = io_.write_all(rec, out)
    check(len(paths) == 4, 'four files written')
    for p in paths:
        with open(p, newline='') as f:
            rows = list(csv.reader(f))
        header, body = rows[0], rows[1:]
        widths = {len(r) for r in body}
        name = os.path.basename(p)
        check(bool(body), '%-34s has %d data row(s)' % (name, len(body)))
        check(widths == {len(header)},
              '%-34s every row is %d wide%s'
              % (name, len(header),
                 '' if widths == {len(header)} else ' — found widths %s' % sorted(widths)))

    print('\n=== the values are the right ones ===')
    with open(os.path.join(out, 'bsrt_P001_baseline_t1_raw.csv'), newline='') as f:
        raws = list(csv.DictReader(f))
    check(raws[0]['epoch_index'] == '1', 'epoch index is 1-based')
    check(raws[0]['minute'] == '1', 'minute is 1-based')
    check(raws[17]['responded'] == '0' and raws[17]['rt_ms'] == '',
          'a missed epoch has no RT rather than a zero')
    check(raws[41]['outcome'] == 'miss' and raws[41]['late_response'] == '1',
          'the 1400 ms response is a late miss but keeps its RT: %s' % raws[41]['rt_ms'])
    check(raws[25]['false_start'] == '1', 'the 60 ms response is a false start')
    check(raws[0]['rt_source_used'] == 'psychopy_flip_to_kb',
          'the RT source names this build')

    with open(os.path.join(out, 'bsrt_P001_baseline_t1_summary.csv'), newline='') as f:
        summ = list(csv.DictReader(f))[0]
    check(summ['mode'] == 'bsrt', 'mode recorded')
    check(summ['schedule_seed'] == str(CFG['seed']), 'the seed is exported, so the schedule is reproducible')
    check(summ['scheduled_stimuli'] == '60', '60 scheduled stimuli')
    check(summ['correction'] == 'false_starts+outliers', 'the correction is named')
    check(summ['kss_before'] == '4' and summ['kss_after'] == '7', 'KSS answers exported')
    check(summ['frames_trial'] != '' and summ['dropped_frames_trial'] != '',
          'frame counts exported (%s frames, %s dropped)'
          % (summ['frames_trial'], summ['dropped_frames_trial']))
    check(summ['device_platform'] == 'psychopy', 'the platform says which build wrote this')
    check(summ['device_browser'] == '', 'browser-only columns are empty, not invented')
    check(summ['end_reason'] == 'max_duration', 'end reason recorded')

    with open(os.path.join(out, 'bsrt_P001_baseline_t1_norms.csv'), newline='') as f:
        nrows = list(csv.DictReader(f))
    check(len(nrows) == 21, 'one norms row per variable (got %d)' % len(nrows))
    check(nrows[0]['norm_available'] == '1', 'a comparison was available at 03:00')
    check(nrows[0]['norm_window_min'] == '3', 'compared over 3 minutes')
    check(nrows[0]['norm_ref_n_control'] != '' and nrows[0]['norm_ref_n_overnight'] != '',
          'the reference cell composition travels with the row (%s control, %s overnight)'
          % (nrows[0]['norm_ref_n_control'], nrows[0]['norm_ref_n_overnight']))
    bands = {r['band'] for r in nrows}
    check(bands <= {'green', 'orange', 'red', 'neutral', ''},
          'bands are the expected set: %s' % sorted(bands))
    ipr = [r for r in nrows if r['variable'] == 'RT interpercentile range']
    check(len(ipr) == 1 and ipr[0]['norm_ref_mean'] != '',
          'the IPR is compared against a reference mean')

    print('\n=== a PVT gets no normative comparison ===')
    rec2 = dict(rec, norms=bn.normative_report({
        'mode': 'pvt', 'isiMs': 3000, 'hour': 3, 'elapsedMs': raw['elapsedMs'],
        'minuteBuckets': 3, 'trials': scored['trials']}, norms))
    rows2 = io_.norms_rows(rec2)
    check(len(rows2) == 1, 'a single row saying so')
    check(len(rows2[0]) == len(io_.NORMS_HEADER),
          'still the full width (%d vs %d)' % (len(rows2[0]), len(io_.NORMS_HEADER)))
    check(rows2[0][12] == 'not_bsrt', 'with the reason: %s' % rows2[0][12])
finally:
    shutil.rmtree(out, ignore_errors=True)

print()
if failures:
    print('%d EXPORT TEST(S) FAILED' % len(failures))
    sys.exit(1)
print('ALL EXPORT TESTS PASSED')
