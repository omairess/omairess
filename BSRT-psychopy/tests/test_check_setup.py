"""check_setup.py, run against the stub.

The setup checker is the first thing the user will run, so it should not be the
first thing that crashes. This drives it end to end with scripted keypresses.
Same caveat as test_psychopy_path.py: the stub is not PsychoPy.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'stub_psychopy'))
sys.path.insert(0, os.path.dirname(HERE))

import psychopy.core as pcore                       # noqa: E402
from psychopy.hardware import keyboard as pkb       # noqa: E402

import check_setup                                  # noqa: E402

failures = []


def check(cond, msg):
    print(('  ok    ' if cond else '  FAIL  ') + msg)
    if not cond:
        failures.append(msg)


FRAME = 1.0 / 60.0

# One key to start, then a press shortly after each of the three dots. The
# checker waits 1 s of blank frames before each dot, so the dots land roughly
# 1 s apart once it gets going; pressing every 100 ms from then on guarantees
# a press is waiting inside each 5 s response window.
keys = [(0.0, 'space')]
t = 0.5
while t < 20.0:
    keys.append((t, 'space'))
    t += 0.10
pkb.Keyboard.SCHEDULE = sorted(keys)
pcore.reset()

print('=== check_setup runs end to end ===')
rc = check_setup.main()
out = '\n'.join(check_setup.REPORT)

check(isinstance(rc, int), 'it returns an exit code (%r)' % rc)
check('PsychoPy imports' in out, 'it reports the PsychoPy version')
check('win.flip() returns a timestamp' in out, 'it checks the flip timestamp')
check('refresh rate measured' in out, 'it measures the refresh rate')
check('frames are not being dropped' in out, 'it checks for dropped frames')
check('share a clock' in out or 'no reaction times' in out,
      'it reaches the key/flip clock check')
check('measured RTs' in out, 'and reports the measured RTs')

print('\n=== the reaction-time verdict ===')
# The stub cannot press a key like a person, so the decision rule is tested
# directly rather than through a fake display that could only ever hit one
# branch of it.
for rts, want, label in [
    ([250.0, 310.0, 290.0], True, 'ordinary human reaction times pass'),
    ([81.0, 2999.0, 500.0], True, 'and so does everything just inside the band'),
    ([-950.0, -960.0], False, 'negative RTs fail — the clocks disagree'),
    ([66.7, 66.7, 66.7], False, 'impossibly fast RTs fail'),
    ([4000.0, 250.0], False, 'implausibly slow RTs fail'),
    ([], False, 'no RTs at all fails'),
]:
    got, detail, advice = check_setup.rt_verdict(rts)
    check(got is want, '%s (%s)' % (label, detail))
    check(want or advice, 'and a failure explains what it means')

report = os.path.join(os.path.dirname(HERE), 'bsrt_setup_report.txt')
check(os.path.exists(report), 'it wrote bsrt_setup_report.txt')
if os.path.exists(report):
    os.unlink(report)

print()
if failures:
    print('%d CHECK-SETUP TEST(S) FAILED' % len(failures))
    sys.exit(1)
print('CHECK_SETUP RUNS END TO END')
