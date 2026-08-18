"""The look of the PsychoPy build: the browser app's palette, drawn on a window.

The other two builds are styled with CSS. PsychoPy has no CSS, so the same
design is rebuilt here out of rectangles and text — the same background, panel,
accent and band colours, taken from BSRT/styles.css, and the same layout ideas:
a titled card on a dark ground, muted secondary text, a key hint along the
bottom, and the norm bands carrying a colour AND a coloured edge AND a word, so
the results survive greyscale and red/green colour blindness exactly as they do
in the browser.

Colours are given in rgb255 throughout. PsychoPy's default colour space is a
-1..1 float space where '#5b8def' either works or silently does something else
depending on the version; naming the space removes the guess.

Everything here draws, and nothing here decides. The task logic lives in
bsrt_task.py and the scoring in bsrt_core.py, so a change to the look cannot
change a number.
"""

# From BSRT/styles.css — one palette across the three builds.
BG = (18, 20, 26)          # --bg      #12141a
PANEL = (27, 30, 38)       # --panel   #1b1e26
LINE = (46, 51, 64)        # --line    #2e3340
TEXT = (232, 234, 240)     # --text    #e8eaf0
MUTED = (154, 161, 177)    # --muted   #9aa1b1
ACCENT = (91, 141, 239)    # --accent  #5b8def
DANGER = (224, 102, 102)   # --danger  #e06666
LED = (255, 60, 40)        # --led     #ff3c28

BAND = {
    'green': (111, 207, 151),
    'orange': (224, 160, 102),
    'red': (235, 111, 111),
    'neutral': MUTED,
    '': MUTED,
}


class UI(object):
    """Drawing helpers bound to one window.

    Stimuli are created once and reused: building a TextStim costs milliseconds
    and rebuilding one per frame is the usual reason a PsychoPy screen stutters.
    """

    def __init__(self, win):
        from psychopy import visual

        self.win = win
        self.visual = visual
        self._text = {}
        self._rects = []
        self._n_rects = 0

    # -- primitives ---------------------------------------------------------

    def text(self, key, **kw):
        """A cached TextStim, so the same element is not rebuilt every frame."""
        if key not in self._text:
            opts = dict(text='', height=0.038, color=TEXT, colorSpace='rgb255',
                        units='height', wrapWidth=1.5, alignText='center',
                        anchorHoriz='center')
            opts.update(kw)
            self._text[key] = self.visual.TextStim(self.win, **opts)
        stim = self._text[key]
        for k, v in kw.items():
            try:
                setattr(stim, k, v)
            except Exception:
                pass
        return stim

    def rect(self, **kw):
        """A pooled Rect. Pooling keeps a table of 20 rows from building 20
        objects on every single frame."""
        if self._n_rects >= len(self._rects):
            self._rects.append(self.visual.Rect(
                self.win, width=1.0, height=0.1, units='height',
                fillColor=PANEL, colorSpace='rgb255', lineColor=None))
        r = self._rects[self._n_rects]
        self._n_rects += 1
        for k, v in kw.items():
            try:
                setattr(r, k, v)
            except Exception:
                pass
        return r

    def begin(self):
        """Start a frame: hand every pooled rectangle back."""
        self._n_rects = 0

    # -- composites ---------------------------------------------------------

    def panel(self, y, height, width=1.15, fill=PANEL):
        self.rect(width=width, height=height, pos=(0, y),
                  fillColor=fill, lineColor=LINE, lineWidth=1).draw()

    def title(self, s, y=0.36):
        self.text('title', text=s, height=0.055, color=TEXT, pos=(0, y)).draw()

    def subtitle(self, s, y=0.30):
        self.text('subtitle', text=s, height=0.032, color=MUTED, pos=(0, y)).draw()

    def hint(self, s, y=-0.44):
        """The key hint along the bottom, as in the browser build's footer."""
        self.text('hint', text=s, height=0.028, color=MUTED, pos=(0, y)).draw()

    def accent_rule(self, y, width=1.15):
        self.rect(width=width, height=0.004, pos=(0, y),
                  fillColor=ACCENT, lineColor=None).draw()

    def body(self, s, y=0.0, height=0.038, color=TEXT, width=1.5):
        self.text('body%s' % y, text=s, height=height, color=color,
                  pos=(0, y), wrapWidth=width).draw()

    def kv_row(self, label, value, y, color=TEXT, band=None, x=0.0, w=1.05):
        """One label/value line, with an optional coloured band edge.

        The edge is drawn as well as the colour, never instead of it, for the
        same reason the browser build draws a left border: colour alone is not
        readable to everyone and does not survive a greyscale printout.
        """
        if band:
            c = BAND.get(band, MUTED)
            self.rect(width=0.008, height=0.042, pos=(x - w / 2 - 0.012, y),
                      fillColor=c, lineColor=None).draw()
            color = c
        self.text('k%s' % y, text=label, height=0.030, color=MUTED,
                  pos=(x - w / 2, y), anchorHoriz='left', alignText='left',
                  wrapWidth=w).draw()
        self.text('v%s' % y, text=value, height=0.030, color=color,
                  pos=(x + w / 2, y), anchorHoriz='right', alignText='right',
                  wrapWidth=w).draw()


def clear(ui):
    """A blank themed frame."""
    ui.begin()


def draw_instructions(ui, title, lines, hint):
    ui.begin()
    ui.panel(0.0, 0.62)
    ui.title(title)
    ui.accent_rule(0.315)
    ui.body('\n'.join(lines), y=0.05, height=0.036, width=1.0)
    ui.hint(hint)


def draw_countdown(ui, n, hint):
    """A big number in an accent ring, rather than a bare digit."""
    ui.begin()
    ui.rect(width=0.30, height=0.30, pos=(0, 0.02), fillColor=PANEL,
            lineColor=ACCENT, lineWidth=3).draw()
    ui.text('cd', text='%d' % n, height=0.16, color=TEXT, pos=(0, 0.02)).draw()
    ui.hint(hint)


def draw_kss(ui, title, question, anchors, selected, hint):
    """The 9-point scale drawn as a row of numbered boxes with its anchors.

    A wall of nine numbered lines is what a plain text version looks like; the
    browser build shows a scale, so this does too. The anchor for whichever
    number is highlighted is shown underneath, which is also how the web
    version reads.
    """
    ui.begin()
    ui.panel(0.0, 0.72)
    ui.title(title, y=0.30)
    ui.accent_rule(0.255)
    ui.body(question, y=0.17, height=0.036, width=1.0)

    n = len(anchors)
    box = 0.085
    gap = 0.016
    total = n * box + (n - 1) * gap
    x0 = -total / 2.0 + box / 2.0
    for i in range(n):
        x = x0 + i * (box + gap)
        on = (selected == i + 1)
        ui.rect(width=box, height=box, pos=(x, 0.02),
                fillColor=ACCENT if on else PANEL,
                lineColor=ACCENT if on else LINE, lineWidth=2).draw()
        ui.text('kssn%d' % i, text='%d' % (i + 1), height=0.042,
                color=BG if on else TEXT, pos=(x, 0.02)).draw()

    # The two ends are always labelled, so the direction of the scale is never
    # ambiguous even before anything is highlighted.
    ui.text('kssLo', text=anchors[0], height=0.026, color=MUTED,
            pos=(x0, -0.075), wrapWidth=0.30).draw()
    ui.text('kssHi', text=anchors[-1], height=0.026, color=MUTED,
            pos=(x0 + (n - 1) * (box + gap), -0.075), wrapWidth=0.30).draw()

    if selected:
        ui.text('kssSel', text=anchors[selected - 1], height=0.038,
                color=ACCENT, pos=(0, -0.185), wrapWidth=1.0).draw()
    ui.hint(hint)


def draw_sleep_onset(ui, title, detail, hint):
    ui.begin()
    ui.panel(0.0, 0.42, fill=(40, 22, 22))
    ui.rect(width=1.15, height=0.006, pos=(0, 0.205),
            fillColor=DANGER, lineColor=None).draw()
    ui.text('soT', text=title, height=0.062, color=DANGER, pos=(0, 0.10)).draw()
    ui.text('soS', text=detail, height=0.032, color=TEXT, pos=(0, 0.01)).draw()
    ui.text('soH', text=hint, height=0.030, color=MUTED, pos=(0, -0.09)).draw()


def draw_results(ui, title, rows, norm_rows, footer, heading=''):
    """The end-of-trial screen: the headline numbers, then the norm bands.

    Deliberately short. The browser build can afford a full results screen
    because it is a page you can scroll; here the experimenter wants the few
    numbers that say whether the trial was any good, and the files hold the
    rest.
    """
    ui.begin()
    ui.panel(0.06, 0.78, width=1.25)
    ui.title(title, y=0.395)
    ui.accent_rule(0.355, width=1.25)

    y = 0.30
    for label, value in rows:
        ui.kv_row(label, value, y, w=1.10)
        y -= 0.048

    if norm_rows:
        y -= 0.015
        ui.rect(width=1.10, height=0.002, pos=(0, y + 0.020),
                fillColor=LINE, lineColor=None).draw()
        ui.text('nHdr', text=heading, height=0.026,
                color=MUTED, pos=(0, y - 0.008)).draw()
        y -= 0.055
        for label, value, band in norm_rows:
            ui.kv_row(label, value, y, band=band, w=1.10)
            y -= 0.048

    ui.hint(footer)
