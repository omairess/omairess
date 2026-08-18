"""Windows and stimuli that record what they were asked to draw."""

from . import core

FRAME = 1.0 / 60.0


class _Snapshot(object):
    """What one stimulus looked like on one frame."""

    def __init__(self, stim):
        self.kind = type(stim).__name__
        self.text = getattr(stim, 'text', None)
        self.pos = getattr(stim, 'pos', None)
        self.stim = stim


class _Stim(object):
    def __init__(self, win, **kw):
        self.win = win
        self.__dict__.update(kw)
        self.drawn = 0

    def draw(self):
        self.drawn += 1
        self.win.pending.append(self)

    def setText(self, t):
        self.text = t


class TextStim(_Stim):
    def __init__(self, win, text='', **kw):
        _Stim.__init__(self, win, text=text, **kw)


class Circle(_Stim):
    pass


class Rect(_Stim):
    pass


class Line(_Stim):
    pass


class ImageStim(_Stim):
    pass


# The most recently created window, so a test can inspect what was drawn.
LAST = [None]


class Window(object):
    def __init__(self, fullscr=False, color='black', units='height',
                 allowGUI=False, waitBlanking=True, **kw):
        LAST[0] = self
        self.fullscr = fullscr
        self.color = color
        self.units = units
        self.monitorFramePeriod = FRAME
        self.pending = []
        self.frames = []            # what was on each flipped frame
        self.closed = False
        self.flips = 0

    def flip(self, clearBuffer=True):
        core.advance(FRAME)
        # A frame is a SNAPSHOT of what was on it, not a list of references.
        # The code under test caches and reuses its stimuli, so holding
        # references would make every recorded frame show the last text that
        # object ever carried.
        self.frames.append([_Snapshot(s) for s in self.pending])
        self.pending = []
        self.flips += 1
        return core.getTime()

    def getActualFrameRate(self, **kw):
        return 1.0 / FRAME

    def close(self):
        self.closed = True
