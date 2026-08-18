"""A keyboard that delivers presses scheduled on the virtual clock."""

from .. import core


class KeyPress(object):
    def __init__(self, name, tDown):
        self.name = name
        self.tDown = tDown
        self.rt = tDown
        self.duration = 0.01


class Keyboard(object):
    # (time in seconds, key name), consumed as the virtual clock passes them.
    SCHEDULE = []

    def __init__(self, clock=None, **kw):
        self.clock = clock or core.Clock()
        self._pending = list(Keyboard.SCHEDULE)

    def getKeys(self, keyList=None, waitRelease=False, clear=True):
        now = core.getTime()
        ready = [(t, n) for (t, n) in self._pending if t <= now]
        if keyList:
            ready = [(t, n) for (t, n) in ready if n in keyList]
        for item in ready:
            self._pending.remove(item)
        return [KeyPress(n, t) for (t, n) in ready]

    def clearEvents(self):
        # Real clearEvents empties the buffer, so anything already pressed is
        # discarded. Modelling that matters: code that relies on it would
        # otherwise look correct here and misbehave on a real keyboard.
        now = core.getTime()
        self._pending = [(t, n) for (t, n) in self._pending if t > now]
