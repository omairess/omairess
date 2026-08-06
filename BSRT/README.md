# BSRT — Behavioral Sleep Resistance Task

A software-based behavioral sleep resistance task, implementing the **OSLER
(Oxford Sleep Resistance Test)** paradigm as a standalone web app.

No build step, no server, no dependencies. Open `index.html` and it runs.

---

## The paradigm

OSLER is a behavioural alternative to the MWT (Maintenance of Wakefulness Test)
that requires no EEG. The participant sits in a dark, quiet room and responds to
a dim light that flashes at a fixed interval. Sleep onset is inferred
*behaviourally*, from a run of consecutive failures to respond.

Defaults reproduce the standard protocol:

| Parameter | Default | Meaning |
|---|---|---|
| Stimulus interval | 3000 ms | onset-to-onset spacing of flashes |
| Stimulus duration | 1000 ms | how long the light stays on |
| Miss criterion | 7 | consecutive missed flashes (= 21 s) scored as sleep onset |
| Maximum duration | 40 min | ceiling; trials reaching it are censored |

All four are editable on the setup screen, so the task can also be run under
modified parameters.

---

## Outcome measures

**Sleep latency** is reported under two explicitly labelled definitions, because
the literature is not consistent on this point:

- `latency_first_miss_ms` — time to the onset of the **first** flash in the
  terminating run of misses. This is the app's headline figure.
- `latency_criterion_ms` — time at which the criterion was **confirmed**
  (i.e. `latency_first_miss_ms + 7 × 3000 ms`).

Trials that reach the ceiling without a scored sleep onset are **censored**, not
recorded as a latency of 40 minutes. `slept_before_max = 0` marks these — treat
them as right-censored in analysis (survival methods) rather than as observed
latencies.

Alongside latency, each trial records: response rate, missed-stimulus count,
longest run of misses, mean / median / SD of reaction time, count of slow
responses (>1000 ms), extra taps, and a count of times the page lost focus
mid-trial.

---

## Data

Trials are saved to **`localStorage` in that browser only**. Nothing is
transmitted anywhere. Clearing site data deletes them — **export after every
session**.

Two CSV shapes, both long-format and ready for `read.csv()` in R:

- **Epoch-level** — one row per flash: `epoch_index`, `responded`, `rt_ms`,
  `extra_responses`, plus participant/session/trial keys on every row so
  multiple files `rbind` cleanly.
- **Summary-level** — one row per trial, with all scored measures and the
  parameters the trial was run under.

"Export all" produces a single combined file across every saved trial.

### Research use

This app keeps data on the participant's device, which is privacy-friendly by
default. But once you are collecting data *from participants for research*,
consent and data-governance obligations (GDPR and your ethics approval) apply to
the exported CSVs regardless of where they were generated. This is a research
instrument, not a diagnostic device.

---

## A more accurate desktop version exists

If reaction-time precision matters to your study, use `../BSRT-desktop` instead.
It is the same paradigm packaged as a native macOS/Windows application, with
frame-locked stimulus presentation, a display calibration routine, and per-trial
timing provenance exported with the data. This browser version remains the right
choice when you need participants to run the task on their own devices without
installing anything.

## Timing accuracy

Stimuli are scheduled against **absolute targets** (`startTime + n × ISI`) rather
than by chaining relative delays, so presentation does not drift over a
40-minute trial. Reaction times are measured with `performance.now()` from the
moment the stimulus is rendered.

Browser timing is nonetheless less precise than dedicated hardware. Reaction
times carry display and input latency of a few tens of milliseconds. This is
immaterial for the primary outcome — sleep latency is scored in 3-second
epochs — but treat the RT distributions as relative measures, and do not compare
them across different devices or browsers.

The app requests a **screen wake lock** and **fullscreen** at trial start (both
are best-effort; browsers may decline). Navigating away mid-trial is counted and
reported.

---

## Running it

**Locally:** open `index.html` in a browser. That is the whole procedure.

**Hosting it** (so participants can load it on their own device): drag this
folder onto <https://app.netlify.com/drop>. It goes live at a `*.netlify.app`
URL. Sign in to claim the site — unclaimed anonymous deploys are temporary. The
same folder also works on GitHub Pages or Cloudflare Pages; all asset paths are
relative.

To update a deployed site, change the files and drag the folder onto the drop
zone on the site's Deploys page again.

---

## Reference

Bennett, L. S., Stradling, J. R., & Davies, R. J. O. (1997). A behavioural test
to assess daytime sleepiness in obstructive sleep apnoea. *Journal of Sleep
Research*, 6(2), 142–145.
