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

    def play(self, loops=0):
        import psychopy as _p
        if _p.FAIL_SOUND[0]:
            raise RuntimeError('no audio device available')
        self.playing = True
        self.plays += 1
        self.loops = loops

    def stop(self):
        self.playing = False
        self.stops += 1
