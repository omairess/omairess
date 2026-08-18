"""bsrt_core.py must score identically to scoring.js.

Three builds of one instrument are only one instrument if they agree on every
number. This drives the Python port and the JavaScript original with the same
inputs and compares the outputs field by field.

Run:  python3 BSRT-psychopy/tests/test_equivalence.py
Needs node on PATH and BSRT/scoring.js in place; nothing else.
"""

import json
import math
import os
import random
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.dirname(HERE))

import bsrt_core as C  # noqa: E402

SCORING_JS = os.path.join(ROOT, 'BSRT', 'scoring.js')
TOL = 1e-9

failures = []


def check(cond, msg):
    if cond:
        print('  ok    ' + msg)
    else:
        print('  FAIL  ' + msg)
        failures.append(msg)


def run_js(script):
    """Evaluate a snippet against scoring.js and return its JSON result."""
    with tempfile.NamedTemporaryFile('w', suffix='.js', delete=False) as f:
        f.write("const S = require(%s);\n" % json.dumps(SCORING_JS))
        f.write(script)
        path = f.name
    try:
        out = subprocess.run([node_bin(), path], capture_output=True, text=True)
        if out.returncode != 0:
            raise RuntimeError('node failed: ' + out.stderr[:800])
        return json.loads(out.stdout)
    finally:
        os.unlink(path)


def node_bin():
    for p in ('/opt/node22/bin/node', 'node'):
        if p == 'node' or os.path.exists(p):
            return p
    return 'node'


def close(a, b, path=''):
    """Deep comparison that treats None/null alike and floats within TOL.

    Returns the first differing path, or None when the two agree.
    """
    if isinstance(a, dict) and isinstance(b, dict):
        for k in set(a) | set(b):
            if k not in a:
                return path + '.' + k + ' (missing in python)'
            if k not in b:
                return path + '.' + k + ' (missing in js)'
            d = close(a[k], b[k], path + '.' + k)
            if d:
                return d
        return None
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return path + ' length %d vs %d' % (len(a), len(b))
        for i, (x, y) in enumerate(zip(a, b)):
            d = close(x, y, path + '[%d]' % i)
            if d:
                return d
        return None
    if a is None or b is None:
        return None if a is None and b is None else path + ' %r vs %r' % (a, b)
    if isinstance(a, bool) or isinstance(b, bool):
        return None if bool(a) == bool(b) else path + ' %r vs %r' % (a, b)
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if a == b:
            return None
        scale = max(1.0, abs(a), abs(b))
        return None if abs(a - b) <= TOL * scale else path + ' %r vs %r' % (a, b)
    return None if a == b else path + ' %r vs %r' % (a, b)


# --------------------------------------------------------------------------
print('=== the seeded schedule generator ===')
# The RNG has to be bit-identical or the same seed gives two different
# schedules, and a session stops being reproducible from its CSV.
seeds = [0, 1, 42, 12345, 2147483647, 4294967295]
js_rng = run_js("""
function makeRng(seed){let a=seed>>>0;return function(){a=(a+0x6D2B79F5)|0;
let t=Math.imul(a^(a>>>15),1|a);t=(t+Math.imul(t^(t>>>7),61|t))^t;
return ((t^(t>>>14))>>>0)/4294967296;};}
const out={};
for (const s of %s){const r=makeRng(s);out[s]=Array.from({length:20},()=>r());}
console.log(JSON.stringify(out));
""" % json.dumps(seeds))
py_rng = {}
for s in seeds:
    r = C.make_rng(s)
    py_rng[str(s)] = [r() for _ in range(20)]
d = close(py_rng, js_rng, 'rng')
check(d is None, 'mulberry32 matches over 120 draws from 6 seeds' + (' — ' + d if d else ''))

print('\n=== schedules ===')
cases = [
    {'mode': 'bsrt', 'isiMs': 3000, 'maxMs': 180000, 'blockMs': 30000,
     'isiSetMs': [3000], 'seed': 7},
    {'mode': 'bsrt', 'isiMs': 2500, 'maxMs': 60000, 'blockMs': 30000,
     'isiSetMs': [2500], 'seed': 99},
    {'mode': 'pvt', 'isiMs': 3000, 'maxMs': 180000, 'blockMs': 30000,
     'isiSetMs': [2000, 4000, 6000, 8000, 10000], 'seed': 12345},
    {'mode': 'pvt', 'isiMs': 3000, 'maxMs': 600000, 'blockMs': 30000,
     'isiSetMs': [2000, 4000, 6000, 8000, 10000], 'seed': 1},
    # A set that does NOT sum to the block length: must fall back to sampling.
    {'mode': 'pvt', 'isiMs': 3000, 'maxMs': 120000, 'blockMs': 30000,
     'isiSetMs': [3000, 5000, 7000], 'seed': 4},
]
js_sched = run_js("console.log(JSON.stringify(%s.map(o => S.buildSchedule(o))));"
                  % json.dumps(cases))
for i, (o, js) in enumerate(zip(cases, js_sched)):
    py = C.build_schedule(o)
    d = close(py, js, 'schedule[%d]' % i)
    label = '%s seed=%s max=%ds' % (o['mode'], o['seed'], o['maxMs'] // 1000)
    check(d is None, '%-26s %3d stimuli, method %s%s'
          % (label, py['nStimuli'], py['method'], ' — ' + d if d else ''))

print('\n=== scoring, on random trials ===')
rng = random.Random(20260818)


def make_epochs(n, miss_p, seed_shift=0):
    """A plausible trial: mostly hits, some misses, occasional late responses
    and false starts — the shapes that separate the two implementations."""
    r = random.Random(1000 + seed_shift)
    sched = C.build_schedule({'mode': 'bsrt', 'isiMs': 3000, 'maxMs': n * 3000,
                              'blockMs': 30000, 'isiSetMs': [3000], 'seed': 5})
    eps = []
    for i in range(min(n, sched['nStimuli'])):
        u = r.random()
        if u < miss_p:
            rt = None
        elif u < miss_p + 0.03:
            rt = r.uniform(1000.001, 2900)      # responded, but too late
        elif u < miss_p + 0.05:
            rt = r.uniform(1, 99)               # false start
        else:
            rt = r.uniform(180, 950)
        eps.append({'index': i, 'onsetMs': sched['onsets'][i],
                    'minute': sched['minutes'][i], 'block': sched['blocks'][i],
                    'isiBeforeMs': sched['isiBefore'][i],
                    'epochIsiMs': sched['epochIsi'][i], 'rtMs': rt})
    return eps


cfgs = [
    {'hitWindowMs': 1000, 'isiMs': 3000, 'lapseMs': 500, 'missCriterion': 7,
     'removeFalseStarts': False, 'removeOutliers': False},
    {'hitWindowMs': 1000, 'isiMs': 3000, 'lapseMs': 500, 'missCriterion': 7,
     'removeFalseStarts': True, 'removeOutliers': True, 'sdMultiplier': 2},
    {'hitWindowMs': 1000, 'isiMs': 3000, 'lapseMs': 355, 'missCriterion': 3,
     'removeFalseStarts': True, 'removeOutliers': False, 'falseStartMs': 150},
]
for ci, cfg in enumerate(cfgs):
    for n, miss_p in ((60, 0.02), (200, 0.15), (400, 0.35)):
        eps = make_epochs(n, miss_p, seed_shift=ci * 10 + n)
        presses = [{'tMs': e['onsetMs'] + (e['rtMs'] or 0), 'epochIndex': e['index']}
                   for e in eps if e['rtMs'] is not None]
        js = run_js("console.log(JSON.stringify(S.score(%s, %s, %s)));"
                    % (json.dumps(eps), json.dumps(cfg), json.dumps(presses)))
        py = C.score(eps, cfg, presses)
        d = close(py, js, 'score')
        check(d is None, 'cfg%d n=%-3d miss=%2d%%  hits=%-3d lapses=%-3d late=%-2d%s'
              % (ci, n, int(miss_p * 100), py['totals']['hits'], py['totals']['lapses'],
                 py['totals']['lateResponses'], ' — ' + d if d else ''))

print('\n=== edge cases ===')
edge = [
    ('no trials at all', [], cfgs[0]),
    ('every trial missed', [{'index': i, 'onsetMs': (i + 1) * 3000, 'minute': i // 20,
                             'isiBeforeMs': 3000, 'epochIsiMs': 3000, 'rtMs': None}
                            for i in range(40)], cfgs[0]),
    ('a single trial', [{'index': 0, 'onsetMs': 3000, 'minute': 0, 'isiBeforeMs': 3000,
                         'epochIsiMs': 3000, 'rtMs': 312.5}], cfgs[0]),
    ('two identical RTs (zero SD)', [{'index': i, 'onsetMs': (i + 1) * 3000, 'minute': 0,
                                      'isiBeforeMs': 3000, 'epochIsiMs': 3000, 'rtMs': 300.0}
                                     for i in range(2)], cfgs[1]),
    ('an anticipation at exactly 0', [{'index': i, 'onsetMs': (i + 1) * 3000, 'minute': 0,
                                       'isiBeforeMs': 3000, 'epochIsiMs': 3000,
                                       'rtMs': 0.0 if i == 0 else 300.0}
                                      for i in range(5)], cfgs[0]),
    ('exactly on the hit window', [{'index': i, 'onsetMs': (i + 1) * 3000, 'minute': 0,
                                    'isiBeforeMs': 3000, 'epochIsiMs': 3000,
                                    'rtMs': 1000.0 if i == 0 else 300.0}
                                   for i in range(5)], cfgs[0]),
    # 25 hits makes the decile exactly 2.5 observations: JS Math.round takes 3,
    # Python's banker's rounding would take 2.
    ('a decile of exactly 2.5', [{'index': i, 'onsetMs': (i + 1) * 3000, 'minute': 0,
                                  'isiBeforeMs': 3000, 'epochIsiMs': 3000,
                                  'rtMs': 200.0 + i}
                                 for i in range(25)], cfgs[0]),
]
for name, eps, cfg in edge:
    js = run_js("console.log(JSON.stringify(S.score(%s, %s, [])));"
                % (json.dumps(eps), json.dumps(cfg)))
    py = C.score(eps, cfg, [])
    d = close(py, js, 'score')
    check(d is None, name + (' — ' + d if d else ''))

print('\n=== the normative re-scoring ===')
for n, miss_p in ((60, 0.05), (160, 0.25)):
    eps = make_epochs(n, miss_p, seed_shift=77 + n)
    scored = C.score(eps, cfgs[0], [])
    for mins in (1, 3, 8):
        js = run_js("console.log(JSON.stringify(S.normativeSummary(%s, %d, %s)));"
                    % (json.dumps(scored['trials']), mins, json.dumps(cfgs[0])))
        py = C.normative_summary(scored['trials'], mins, cfgs[0])
        d = close(py, js, 'normativeSummary')
        check(d is None, 'n=%-3d over the first %d minute(s)%s'
              % (n, mins, ' — ' + d if d else ''))

print('\n=== the integrity check ===')
for name, presses in (
    ('genuine responding', [{'tMs': 3000 * (i + 1) + 300, 'epochIndex': i} for i in range(40)]),
    ('a held key', [{'tMs': i * 30, 'epochIndex': i // 100} for i in range(400)]),
    ('presses outside any epoch', [{'tMs': i * 700, 'epochIndex': None} for i in range(30)]),
    ('nothing at all', []),
):
    js = run_js("console.log(JSON.stringify(S.detectCheating(%s, 40, %s)));"
                % (json.dumps(presses), json.dumps(cfgs[0])))
    py = C.detect_cheating(presses, 40, cfgs[0])
    d = close(py, js, 'integrity')
    check(d is None, '%-26s suspected=%s%s'
          % (name, py['suspected'], ' — ' + d if d else ''))

print()
if failures:
    print('%d EQUIVALENCE CHECK(S) FAILED' % len(failures))
    sys.exit(1)
print('PYTHON AND JAVASCRIPT AGREE ON EVERY FIELD')
