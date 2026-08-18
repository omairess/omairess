"""Virtual clock. Time only moves when a frame is flipped."""

_now = [0.0]


def getTime():
    return _now[0]


def advance(seconds):
    _now[0] += seconds


def reset():
    _now[0] = 0.0


class Clock(object):
    def __init__(self):
        self._t0 = _now[0]

    def getTime(self):
        return _now[0] - self._t0

    def reset(self):
        self._t0 = _now[0]


class CountdownTimer(object):
    def __init__(self, start=0):
        self._end = _now[0] + start

    def getTime(self):
        return self._end - _now[0]


def wait(secs):
    advance(secs)


def quit():
    raise SystemExit(0)
