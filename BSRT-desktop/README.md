# BSRT (desktop) — Behavioral Sleep Resistance Task

An OSLER-paradigm sleep resistance task packaged as a native application for
**macOS and Windows**, with frame-locked stimulus presentation and a calibration
routine that measures and reports the timing characteristics of the machine it
is running on.

This is the accuracy-focused sibling of the browser version in `../BSRT`. The
paradigm, scoring, and CSV shape are the same. What differs is how time is
measured.

---

## Building the executables

You need [Node.js](https://nodejs.org) 18 or newer. Everything else is fetched
by `npm install`.

```bash
cd BSRT-desktop
npm install
npm start          # run it without packaging, to check it works
```

Then build for your platform:

```bash
npm run dist:mac   # -> dist/BSRT-1.0.0-arm64.dmg, -x64.dmg, and .zip
npm run dist:win   # -> dist/BSRT Setup 1.0.0.exe and a portable .exe
```

**You must build macOS apps on a Mac and Windows apps on Windows.** Code signing
and platform packaging are not portable. If you need both from one machine, run
the builds in CI (GitHub Actions with a `macos-latest` and a `windows-latest`
job) rather than trying to cross-compile.

### The unsigned-app warnings

These builds are not code-signed, because signing requires a paid Apple
Developer account and a Windows code-signing certificate. Participants will
therefore see a scary dialog on first launch:

- **macOS:** right-click the app → *Open* → *Open*. Or once, in Terminal:
  `xattr -cr /Applications/BSRT.app`
- **Windows:** SmartScreen shows "Windows protected your PC" → *More info* →
  *Run anyway*.

If you are distributing to participants rather than running the task yourself,
budget for signing certificates — talking people through a security warning over
the phone is a poor start to a testing session.

---

## Timing: what this actually fixes

This is the part that matters, so it is stated without marketing.

### What the browser version could not do

- Stimuli were scheduled with `setTimeout`, which is best-effort and drifts
  under load.
- The stimulus onset time was recorded as "when the code set the CSS class",
  which is not when the pixel changed. The gap is up to a full refresh interval
  plus compositing, and it varies.
- Responses were timestamped when JavaScript got round to handling the event,
  which includes any queuing delay in the renderer.
- Nothing knew the display's refresh rate, so none of the above could even be
  quantified.

### What this version does

1. **Frame-locked scheduling.** A single `requestAnimationFrame` loop drives the
   whole trial and counts real vsync intervals. An epoch is a fixed number of
   frames, so presentation is periodic by construction and cannot drift relative
   to the screen. Verified: across a trial, successive stimulus onsets stayed
   within a fraction of one frame of the nominal interval, with no accumulation.

2. **Onset times from the frame clock.** The recorded onset is the timestamp of
   the frame the stimulus was painted into, not the moment a class was added.

3. **Latency from measured timestamps.** Sleep latency is read off the observed
   frame times, never recomputed as `epoch index × nominal ISI`. The calibrated
   frame interval carries a small estimation error, and multiplying it by an
   epoch index would let that error accumulate across a 40-minute trial.

4. **Responses timestamped by the OS.** `event.timeStamp` carries the time the
   platform delivered the event, on the same clock as `performance.now()`.
   Reading it instead of calling `performance.now()` in the handler skips the
   renderer's queuing delay. The app *measures* that delay during calibration,
   reports it, and verifies the assumption — if a future Chromium changed the
   timebase, the app detects it, falls back to handler time, and records that it
   did so rather than silently reporting wrong numbers.

5. **No throttling.** `backgroundThrottling: false` plus the background-timer
   and renderer-backgrounding switches. Without these, Chromium quietly throttles
   `requestAnimationFrame` when it decides a window is idle or occluded, which
   would corrupt a long trial with no visible sign.

6. **Display sleep blocked at the OS level** via `powerSaveBlocker`, which is far
   more reliable than the web Wake Lock API — a necessity in a task whose whole
   point is that the participant stops moving.

### What remains, and cannot be fixed in software

**No software on a general-purpose OS can observe photon emission or physical
key travel.** What is left after the above is:

| Source | Magnitude | Character |
|---|---|---|
| Display pipeline depth (frames buffered between paint and panel) | 1–2 frames (17–33 ms at 60 Hz) | near-constant per machine |
| Panel response / pixel transition | ~1–20 ms, panel-dependent | near-constant |
| Input hardware polling — wired USB | 1–8 ms | jitter |
| Input hardware — Bluetooth keyboard | 10–30 ms, variable | **jitter, avoid** |
| OS → application event delivery | typically <1–3 ms | measured and removed |

So: **jitter is reduced to roughly the input polling interval — single-digit
milliseconds on a wired USB keyboard — and a constant offset of very roughly
20–60 ms remains, unmeasurable in software.**

That constant is the honest limit. It is also the least damaging kind of error:
it is identical for every trial on a given machine, so it cancels in
within-participant and within-device comparisons. It does *not* cancel when
comparing raw reaction times across different hardware, which is why every
number needed to characterise it is exported with the data.

**Sleep latency — the primary outcome — is unaffected by all of this.** It is
scored in 3-second epochs; a few tens of milliseconds is four orders of
magnitude below the measurement grain.

### Measuring the constant, if you need it

Set *Presentation pipeline* and *Measured hardware offset* under **Timing &
calibration**, and enable the **photodiode patch** — a white square in the
top-left corner that flips with the stimulus.

1. Tape a photodiode over the patch (Black Box Toolkit, Cedrus StimTracker, or
   any oscilloscope-plus-photodiode rig).
2. Have the rig also register a keypress, ideally via a response box that closes
   the same circuit.
3. Run a few dozen stimuli. The measured interval between light onset and key
   closure, minus the reaction time the app reports, is your constant offset.
4. Enter it as *Measured hardware offset*. Raw and corrected reaction times are
   both exported; nothing is overwritten.

Until you do this, leave the offset at 0 and report reaction times as raw. The
app never silently applies a correction you did not measure.

### The app tells you when the assumption is wrong

If more than 10% of responses come in under 100 ms, the results screen says so
and names the likely cause. Genuine reactions are essentially never that fast,
so a crop of them means the presentation-pipeline setting overestimates this
display's lag and every reaction time is shifted. Responses under 100 ms are
counted as **anticipations** and excluded from the reaction-time statistics
(standard PVT practice); the raw values are still exported per epoch so you can
apply a different rule in R.

---

## Before you run participants

Frame-locked timing assumes a **fixed** refresh rate. Several common settings
break that assumption, and the calibration check will show it as jitter:

- **Turn off variable refresh rate** — G-Sync, FreeSync, VRR — in the GPU
  control panel.
- **macOS ProMotion displays vary between 24 and 120 Hz.** Set a fixed refresh
  rate in *System Settings → Displays* rather than *ProMotion*.
- **Run on mains power.** Laptops throttle the GPU on battery, which shows up as
  dropped frames.
- **One display only.** Mirroring or extending can pin the compositor to the
  slowest attached screen. The app warns if it detects more than one.
- **Close other applications**, particularly anything doing video or screen
  capture.
- **Use a wired USB keyboard**, not Bluetooth. This is the single largest
  remaining source of *jitter* and it is entirely under your control.

Run the display check once on each machine before testing and keep the numbers.
If it grades `poor`, fix the machine rather than the data.

---

## What gets recorded

Every trial carries its own timing provenance, so a reviewer can see exactly
what the numbers mean. The summary CSV includes, alongside the usual protocol
and performance columns:

`refresh_hz_measured`, `refresh_hz_reported`, `frame_interval_ms`,
`frame_jitter_mad_ms`, `onset_quantisation_ms`, `dropped_frames_calibration`,
`dropped_frames_trial`, `calibration_grade`, `input_time_source`,
`input_dispatch_median_ms`, `presentation_offset_frames`, `hardware_offset_ms`,
`isi_achieved_ms`, `n_rt_valid`, `anticipations`, `platform`, `os_release`,
`electron_version`, `chrome_version`, `display_scale_factor`.

Note `refresh_hz_measured` versus `refresh_hz_reported`: the first is what the
compositor actually achieved, the second is what the OS claims. When they
disagree, trust the measured one and investigate the machine.

`isi_achieved_ms` is the real inter-stimulus interval, which is the requested
interval rounded to a whole number of frames — at 60 Hz, a requested 3000 ms
becomes 2999.7 ms or 3016.7 ms. Report the achieved value, not the requested one.

Data is stored locally per machine and exported through a normal save dialog.
Both CSV shapes are long-format with participant/session/trial keys on every
row, so files `rbind` cleanly in R.

---

## When you should use something else

If you need **sub-millisecond accuracy, hardware-triggered synchronisation, or
EEG/polysomnography event markers**, use [PsychoPy](https://psychopy.org) with a
Cedrus response box or a Black Box Toolkit instead. That is the field standard
and it can do things a Chromium-based app structurally cannot.

For this paradigm the extra precision buys little: sleep latency is scored in
3-second epochs, and sleepiness-related reaction-time slowing operates on the
scale of tens to hundreds of milliseconds, comfortably above this app's residual
jitter. But if a reviewer asks for photodiode-verified timing, the honest answer
is to verify it with a photodiode — which is why the patch is built in.

---

## Reference

Bennett, L. S., Stradling, J. R., & Davies, R. J. O. (1997). A behavioural test
to assess daytime sleepiness in obstructive sleep apnoea. *Journal of Sleep
Research*, 6(2), 142–145.
