"""A stand-in for PsychoPy, so the real code path can be run without a screen.

WHAT THIS DOES AND DOES NOT PROVE. It exercises the branch of bsrt_psychopy.py
that talks to PsychoPy: the window/keyboard/sound calls, the instruction and KSS
screens, the countdown, the alarm-and-silence flow, and the export at the end.
It catches typos, wrong attribute names and flow bugs in OUR code.

It cannot prove that the PsychoPy API is used correctly, because it is not
PsychoPy — it is a stub written from the same reading of the docs as the code
under test, so a misreading would be reproduced in both. Only a run on a real
machine with PsychoPy installed settles that. See BSRT-psychopy/README.md, which
says so plainly.

The virtual clock advances one frame per flip, so timings are exact and the
tests are deterministic.
"""

__version__ = '0.0-stub'


class _Prefs(object):
    """Just enough of psychopy.prefs for the audio-backend selection."""

    def __init__(self):
        self.hardware = {'audioLib': ['ptb']}


prefs = _Prefs()

# Set by a test to imitate an audio backend that constructs fine and only fails
# when asked to play. It lives HERE rather than in sound.py because the code
# under test reloads sound.py to switch backends, and a reload would reset a
# flag defined there.
FAIL_SOUND = [False]
