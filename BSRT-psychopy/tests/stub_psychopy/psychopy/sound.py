"""A sound that records whether it was played and stopped."""

# The most recently created sound, so a test can check the alarm behaviour.
LAST = [None]


class Sound(object):
    def __init__(self, value=440, secs=0.5, stereo=True, **kw):
        LAST[0] = self
        self.value = value
        self.secs = secs
        self.playing = False
        self.plays = 0
        self.stops = 0

    # Set by a test to imitate an audio backend that constructs fine and only
    # fails when actually asked to play — the common macOS failure.
    FAIL_ON_PLAY = [False]

    def play(self, loops=0):
        if Sound.FAIL_ON_PLAY[0]:
            raise RuntimeError('no audio device available')
        self.playing = True
        self.plays += 1
        self.loops = loops

    def stop(self):
        self.playing = False
        self.stops += 1
