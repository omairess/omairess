#!/usr/bin/env python
"""Check that PsychoPy is set up well enough to collect BSRT data.

    python check_setup.py

Run this FIRST, before any participant. It takes about a minute, asks you to
press a key a few times, and prints a verdict. It also writes
`bsrt_setup_report.txt` next to itself — send that file if anything fails.

It checks the four things that decide whether the reaction times are real:

  1. PsychoPy imports, and which keyboard backend it has.
  2. `win.flip()` returns a timestamp. The task treats that as the moment the
     stimulus reached the screen; if it returns nothing, every onset is wrong.
  3. Key timestamps and flip timestamps share a clock. A reaction time is one
     minus the other, so if they are on different clocks every RT is off by a
     constant — often by a lot, sometimes negative. This is the check that
     matters most and the easiest one to get silently wrong.
  4. The refresh rate is measured, and frames are not being dropped.

None of this is specific to BSRT; it is the same set of things any
reaction-time experiment depends on.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

REPORT = []


def say(line=''):
    print(line)
    REPORT.append(line)


def verdict(ok, label, detail=''):
    say('%s  %s%s' % ('PASS' if ok else 'FAIL', label, ('  — ' + detail) if detail else ''))
    return ok


def rt_verdict(rts):
    """Judge measured reaction times, and say what a bad one means.

    Separated from the screen code so the rule can be tested directly: the stub
    used in tests cannot press a key like a person, so driving this through a
    fake display would only ever exercise one branch.

    A human simple reaction time to a visual stimulus is roughly 200-500 ms.
    Below 80 ms nobody has reacted to anything; above 3 s the press did not
    belong to that stimulus. Either means the two clocks disagree.
    """
    if not rts:
        return False, 'no reaction times were measured', [
            'The dot appeared but no keypress was recorded within five',
            'seconds. Click the task window once so it has keyboard focus,',
            'then run this again.',
        ]
    negative = [r for r in rts if r < 0]
    if negative:
        return False, 'some reaction times are NEGATIVE', [
            'A negative reaction time means the key was stamped BEFORE the',
            'stimulus appeared, so the keyboard clock and the window clock',
            'have different origins. Every RT this build records would be',
            'wrong by that offset. Report this.',
        ]
    too_fast = [r for r in rts if r <= 80.0]
    if too_fast:
        return False, 'some reaction times are impossibly short', [
            'Under 80 ms is faster than a person can react to something they',
            'have seen, so these are not real responses to the dot. Either a',
            'stale keypress was counted, or the two clocks disagree.',
        ]
    too_slow = [r for r in rts if r >= 3000.0]
    if too_slow:
        return False, 'some reaction times are implausibly long', [
            'Over 3 s means the press was probably not a response to that',
            'dot at all. Try again, pressing as soon as you see it; if it',
            'persists, the two clocks disagree.',
        ]
    return True, 'all %d are plausible reaction times' % len(rts), []

def main():
    say('BSRT setup check')
    say('=' * 60)
    say()

    # ---- 1. PsychoPy is installed ----
    try:
        import psychopy
        from psychopy import visual, core, gui                     # noqa: F401
        from psychopy.hardware import keyboard
    except Exception as e:
        say('FAIL  PsychoPy could not be imported')
        say('      %s: %s' % (type(e).__name__, e))
        say()
        say('      You are probably running this with the wrong Python.')
        say('      Open the PsychoPy application, choose the Coder view, open')
        say('      this file there and press the green Run button instead.')
        _write()
        return 1

    ok = True
    ok &= verdict(True, 'PsychoPy imports', 'version %s' % psychopy.__version__)
    say('      python: %s' % sys.version.split()[0])

    try:
        import psychtoolbox                                        # noqa: F401
        verdict(True, 'psychtoolbox is available',
                'the keyboard can use hardware timestamps')
    except Exception:
        verdict(True, 'psychtoolbox is NOT installed',
                'PsychoPy will fall back to a less accurate keyboard')
        say('      This is not fatal, but it is the main reason to use this')
        say('      build. Install it if you can: pip install psychtoolbox')
    say()

    # ---- 2. a window, and a frame rate ----
    try:
        win = visual.Window(fullscr=False, size=(900, 650), color='black',
                            units='height', allowGUI=True, waitBlanking=True)
    except Exception as e:
        verdict(False, 'a window could not be opened', '%s: %s' % (type(e).__name__, e))
        _write()
        return 1
    ok &= verdict(True, 'a window opened')

    t = win.flip()
    flip_ok = isinstance(t, (int, float)) and t > 0
    ok &= verdict(flip_ok, 'win.flip() returns a timestamp',
                  ('got %r' % t) if not flip_ok else 'got %.3f s' % t)
    if not flip_ok:
        say('      Every stimulus onset in the task is this value. Stop here')
        say('      and report it — the task would mis-time everything.')

    hz = None
    try:
        hz = win.getActualFrameRate(nIdentical=20, nMaxFrames=240,
                                    nWarmUpFrames=20, threshold=1.0)
    except Exception as e:
        say('      getActualFrameRate raised %s' % type(e).__name__)
    if hz:
        ok &= verdict(hz > 20, 'refresh rate measured', '%.2f Hz (%.2f ms per frame)'
                      % (hz, 1000.0 / hz))
    else:
        verdict(True, 'refresh rate could not be measured',
                'the task will assume 60 Hz — check your display settings')
        hz = 60.0
    frame_ms = 1000.0 / hz
    say()

    # ---- 3. dropped frames over a quiet second ----
    times = []
    for _ in range(120):
        times.append(win.flip())
    gaps = [(b - a) * 1000.0 for a, b in zip(times, times[1:])]
    dropped = sum(1 for g in gaps if g > frame_ms * 1.5)
    ok &= verdict(dropped <= 2, 'frames are not being dropped',
                  '%d of %d intervals were long' % (dropped, len(gaps)))
    if dropped > 2:
        say('      Close other applications and try again. A few drops on a')
        say('      busy machine are normal; dozens are not.')
    say()

    # ---- 4. the clock check: keys and flips together ----
    kb = keyboard.Keyboard()
    dot = visual.Circle(win, radius=0.05, units='height', fillColor='red', lineColor=None)
    msg = visual.TextStim(win, height=0.04, color='white', units='height', wrapWidth=1.4)

    msg.text = ('Timing check\n\n'
                'A red dot will appear three times.\n'
                'Press the SPACE BAR as soon as you see it.\n\n'
                '[ press space to begin ]')
    msg.draw()
    win.flip()
    kb.clearEvents()
    waited = 0
    while not kb.getKeys(waitRelease=False):
        msg.draw()
        win.flip()
        waited += 1
        if waited > int(hz) * 60:
            verdict(False, 'no key was ever detected', 'the keyboard is not reporting presses')
            win.close()
            _write()
            return 1

    rts = []
    for i in range(3):
        # A blank pause, then the dot. Its onset is the flip that carried it.
        for _ in range(int(hz)):
            win.flip()
        kb.clearEvents()
        dot.draw()
        onset = win.flip()

        pressed = None
        for _ in range(int(hz) * 5):
            # Ignore anything stamped before the dot appeared: a key still in
            # the buffer would otherwise produce a negative reaction time and
            # look like a broken clock rather than a stale press.
            for k in kb.getKeys(waitRelease=False):
                if k.tDown >= onset:
                    pressed = k
                    break
            if pressed is not None:
                break
            dot.draw()
            win.flip()
        if pressed is None:
            say('      (no response to dot %d)' % (i + 1))
            continue
        rts.append((pressed.tDown - onset) * 1000.0)

    if not rts:
        ok &= verdict(False, 'no reaction times could be measured')
    else:
        say('      measured RTs: %s ms' % ', '.join('%.1f' % r for r in rts))
        good, detail, advice = rt_verdict(rts)
        ok &= verdict(good, 'key and flip timestamps share a clock', detail)
        for line in advice:
            say('      ' + line)
    say()

    # ---- 5. the sleep-onset alarm ----
    #
    # Reported from a real machine as "does not recognize sound while
    # everything is enabled". PsychoPy binds one audio backend at import and
    # gives up if it fails, so this tries each in turn and says which worked.
    import bsrt_psychopy as app
    alarm, notes = app._make_alarm({'alarm': True})
    for n in notes:
        say('      ' + n)
    if alarm is None:
        verdict(True, 'no audio backend is available',
                'the sleep-onset alarm will be visual only')
        say('      The task still works: the sleep-onset banner takes over the')
        say('      screen and waits to be dismissed. Only the sound is missing.')
    else:
        played = True
        try:
            alarm.play(loops=-1)
        except Exception as e:
            played = False
            say('      playing failed: %s: %s' % (type(e).__name__, e))
        if played:
            msg.text = ('Sound check\n\n'
                        'An alarm should be sounding now.\n\n'
                        'Press Y if you can hear it, N if you cannot.')
            kb.clearEvents()
            heard = None
            for _ in range(int(hz) * 30):
                msg.draw()
                win.flip()
                for k in kb.getKeys(waitRelease=False):
                    if k.name in ('y', 'n'):
                        heard = (k.name == 'y')
                if heard is not None:
                    break
            try:
                alarm.stop()
            except Exception:
                pass
            if heard is None:
                verdict(True, 'sound check skipped', 'no answer given')
            else:
                ok &= verdict(heard, 'the alarm is audible',
                              'you heard it' if heard else 'you did not hear it')
                if not heard:
                    say('      Check the output device and the volume, then run')
                    say('      this again. If it still fails, the task will fall')
                    say('      back to a visual-only alarm, which is safe but')
                    say('      easier to miss from another room.')
    try:
        win.close()
    except Exception:
        pass
    say()

    say('=' * 60)
    if ok:
        say('EVERYTHING PASSED — this machine can collect BSRT data.')
        say()
        say('Next: run  python bsrt_psychopy.py')
    else:
        say('SOMETHING FAILED — see above.')
        say()
        say('Send bsrt_setup_report.txt and it can be fixed against what your')
        say('machine actually does.')
    _write()
    return 0 if ok else 1


def _write():
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'bsrt_setup_report.txt')
    try:
        with open(path, 'w') as f:
            f.write('\n'.join(REPORT) + '\n')
        print('\nreport written to %s' % path)
    except Exception as e:
        print('\ncould not write the report: %s' % e)


if __name__ == '__main__':
    sys.exit(main())
