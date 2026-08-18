"""The BSRT/PVT trial engine, with no PsychoPy import anywhere in it.

The engine drives a `Display` — anything exposing the handful of methods below
— one frame at a time. `bsrt_psychopy.py` supplies the real one; the tests
supply a fake that advances a virtual clock, which is how the whole task can be
exercised on a machine with no screen.

    frame_interval_ms()   nominal ms between flips, from the measured refresh
    flip(state)           draw `state` and block until it is on screen;
                          returns the time in ms at which it appeared
    poll()                [(t_ms, key)] for presses since the last call
    quit_requested()      True when the experimenter asked to abort

`state` is what should be visible on the coming frame:
    ('blank',)            nothing
    ('dot',)              the BSRT stimulus
    ('counter', ms)       the PVT counter, showing `ms`
    ('feedback', ms)      the PVT achieved-RT feedback

TIMING. A stimulus onset is the flip timestamp of the frame that actually
carried it, not the time the code asked for it, and a reaction time is the
keyboard timestamp minus that onset. Both come from the display, so on PsychoPy
they are the vsync time and the PTB keyboard clock. This is the whole reason for
a PsychoPy build: the browser can only ever approximate that.

EPOCH BOUNDARIES COME FROM THE INTENDED SCHEDULE, never from a measured onset.
Closing an epoch a fixed interval after the frame it actually appeared on lets
each epoch's presentation error feed into the start of the next, and the error
compounds: an early version of this loop drifted 25 ms per epoch and ran the
sixtieth stimulus of a 3-minute test 1.5 s late. Onsets are absolute positions
in the schedule, so a late frame costs that one stimulus and nothing after it.
The frame actually used is recorded per epoch, so presentation error is
measurable after the fact rather than assumed away.
"""

import bisect

import bsrt_core as core


def run_trial(display, schedule, cfg, on_sleep_onset=None):
    """Run one whole trial and return the raw record.

    `on_sleep_onset(epoch_index, t_ms)` is called once, if the miss criterion is
    reached, so the caller can start the alarm. The trial then stops.
    """
    stim_ms = cfg['stimMs']
    hit_window_ms = cfg['hitWindowMs']
    miss_criterion = cfg['missCriterion']
    is_pvt = cfg.get('mode') == 'pvt'
    feedback_ms = cfg.get('feedbackMs', 500) if is_pvt else 0

    onsets = schedule['onsets']
    epoch_isi = schedule['epochIsi']
    n = schedule['nStimuli']

    frame_int = display.frame_interval_ms()
    # A frame is 'dropped' when two flips are a whole multiple of the nominal
    # interval apart rather than one of it — the display missed a vsync.
    drop_threshold = frame_int * 1.5

    t0 = None
    last_flip = None
    frames = 0
    dropped = 0

    cur = None              # epoch being presented, or None between epochs
    nxt = 0                 # next epoch to start
    measured = {}           # epoch index -> measured onset, relative to t0
    first_rt = {}           # epoch index -> first response, relative to onset
    raw_presses = []        # (absolute time, key), assigned to epochs afterwards
    closed = []             # epoch indices in the order they finished

    consecutive_misses = 0
    ended = 'max_duration'
    sleep_onset_ms = None

    while True:
        if display.quit_requested():
            ended = 'aborted'
            break

        # Where the coming frame will land. Everything for this iteration is
        # decided against that time, so the stimulus is drawn on the first
        # frame at or after its intended onset rather than one frame later.
        t_next = 0.0 if last_flip is None else (last_flip - t0) + frame_int

        # ---- close any epoch whose intended window has ended ----
        if cur is not None and t_next >= onsets[cur] + epoch_isi[cur]:
            closed.append(cur)
            rt = first_rt.get(cur)
            is_hit = rt is not None and rt <= hit_window_ms
            consecutive_misses = 0 if is_hit else consecutive_misses + 1
            done = cur
            cur = None
            if miss_criterion > 0 and consecutive_misses >= miss_criterion:
                ended = 'sleep_onset'
                sleep_onset_ms = t_next
                if on_sleep_onset:
                    on_sleep_onset(done, t_next)
                break

        # ---- start the epoch whose intended onset has arrived ----
        if cur is None and nxt < n and t_next >= onsets[nxt]:
            cur = nxt
            nxt += 1

        if cur is None and nxt >= n:
            break

        # ---- what the coming frame should show ----
        state = ('blank',)
        if cur is not None:
            if cur not in measured:
                # First frame of this epoch: the stimulus goes up now, and its
                # onset is whatever time this flip actually lands at.
                state = ('counter', 0.0) if is_pvt else ('dot',)
            else:
                since = t_next - measured[cur]
                rt = first_rt.get(cur)
                if rt is not None and feedback_ms and since < rt + feedback_ms:
                    state = ('feedback', rt)
                elif rt is None and since < stim_ms:
                    state = ('counter', since) if is_pvt else ('dot',)

        # ---- present it ----
        t = display.flip(state)
        frames += 1
        if t0 is None:
            t0 = t
        elif t - last_flip > drop_threshold:
            dropped += int(round((t - last_flip) / frame_int)) - 1
        last_flip = t

        if cur is not None and cur not in measured and state[0] in ('dot', 'counter'):
            measured[cur] = t - t0

        # ---- collect responses ----
        for pt, key in display.poll():
            raw_presses.append((pt - t0, key))
            # The first press inside an epoch is that epoch's response. Later
            # ones are extras: recorded for the integrity check, but they never
            # overwrite it.
            j = epoch_of(onsets, epoch_isi, pt - t0)
            if j is not None and j in measured and j not in first_rt:
                first_rt[j] = (pt - t0) - measured[j]

    # Anything still open when the loop ended (the ceiling, or an abort) is
    # closed here so a response in the final epoch is not silently discarded.
    if cur is not None and ended != 'sleep_onset':
        closed.append(cur)

    epochs = []
    for k in closed:
        epochs.append({
            'index': k,
            'onsetMs': measured.get(k),
            'minute': schedule['minutes'][k],
            'block': schedule['blocks'][k],
            'isiBeforeMs': schedule['isiBefore'][k],
            'epochIsiMs': epoch_isi[k],
            'rtMs': first_rt.get(k),
            'intendedOnsetMs': onsets[k],
            'onsetErrorMs': None if measured.get(k) is None else measured[k] - onsets[k],
        })

    presses = []
    for rel, key in raw_presses:
        presses.append({'tMs': rel, 'epochIndex': epoch_of(onsets, epoch_isi, rel)})

    return {
        'epochs': epochs,
        'presses': presses,
        'elapsedMs': (last_flip - t0) if t0 is not None else 0,
        'endReason': ended,
        'sleepOnsetMs': sleep_onset_ms,
        'sleptBeforeMax': ended == 'sleep_onset',
        'frames': frames,
        'droppedFrames': dropped,
        'frameIntervalMs': frame_int,
    }


def epoch_of(onsets, epoch_isi, t_ms):
    """Which epoch's response window contains `t_ms`, or None between epochs.

    Windows are taken from the intended schedule, so a press is attributed the
    same way whatever the presentation error happened to be.
    """
    if not onsets or t_ms < onsets[0]:
        return None
    j = bisect.bisect_right(onsets, t_ms) - 1
    if j < 0:
        return None
    return j if t_ms < onsets[j] + epoch_isi[j] else None


def build_schedule_for(cfg):
    """The schedule for a configured trial."""
    return core.build_schedule({
        'mode': cfg.get('mode', 'bsrt'),
        'isiMs': cfg['isiMs'],
        'isiSetMs': cfg.get('isiSetMs', [2000, 4000, 6000, 8000, 10000]),
        'blockMs': cfg.get('blockMs', 30000),
        'maxMs': cfg['maxMinutes'] * 60000,
        'seed': cfg['seed'],
    })


def score_trial(raw, cfg):
    """Score a completed trial with the shared scoring code."""
    return core.score(raw['epochs'], cfg, raw['presses'])
