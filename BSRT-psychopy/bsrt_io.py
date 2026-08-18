"""The four CSV exports, with the same columns as the browser and desktop builds.

An analysis script should not care which build produced a file, so the headers
here are the headers there — `tests/test_headers.py` extracts them from
BSRT/app.js and asserts every column matches, in order. Add a column in one
build and that test fails until it is added in the others.

Columns that describe a browser (screen_w, device_browser, page_host,
cross_origin_isolated, and the rest) are written empty rather than dropped: the
shape stays identical so the files still rbind, and empty means "this build
cannot know", which is different from zero. The PsychoPy-specific timing figures
go in the columns that already exist for them, documented in README.md.
"""

import csv
import io
import os

PARTICIPANT_COLS = [
    'run_id', 'date', 'time', 'language', 'participant_id', 'name', 'address',
    'birth_date', 'educational_level', 'session_label', 'trial_number',
]

RAW_HEADER = PARTICIPANT_COLS + [
    'epoch_index', 'block', 'minute', 'onset_ms', 'epoch_isi_ms', 'isi_before_ms',
    'responded', 'rt_event_ms', 'rt_handler_ms', 'input_delay_ms', 'rt_source_used',
    'rt_ms', 'rs_per_sec', 'outcome', 'lapse', 'late_response', 'false_start',
    'extra_responses',
]

PM_HEADER = PARTICIPANT_COLS + [
    'minute', 'trials', 'hits', 'misses', 'lapses', 'late_responses', 'hit_ratio',
    'avg_rt', 'median_rt', 'stdev_rt', 'fastest10_rt', 'slowest10_rt',
    'corr_avg_rt', 'corr_median_rt', 'corr_stdev_rt', 'corr_fastest10_rt',
    'corr_slowest10_rt', 'avg_rs', 'median_rs', 'stdev_rs', 'fastest10_rs',
    'slowest10_rs', 'corr_avg_rs', 'corr_median_rs', 'corr_stdev_rs',
    'corr_fastest10_rs', 'corr_slowest10_rs', 'n_rt', 'n_rt_corrected',
    'false_starts_removed', 'outliers_removed', 'velocity_rs', 'acceleration_rs',
    'velocity_rt', 'acceleration_rt',
]

SUMMARY_HEADER = PARTICIPANT_COLS + [
    'sleep_onset_ms', 'sleep_onset_criterion_ms', 'slept_before_max', 'end_reason',
    'elapsed_ms', 'total_trialrun', 'hit_ratio', 'total_hits', 'total_miss',
    'total_lapse', 'late_responses', 'total_false_starts', 'ep_1_2', 'ep_3_6',
    'ep_7plus', 'longest_miss_run', 'total_avg_rt', 'total_median_rt',
    'total_stdev_rt', 'fastest10_rt', 'slowest10_rt', 'corr_total_avg_rt',
    'corr_total_median_rt', 'corr_total_stdev_rt', 'corr_fastest10_rt',
    'corr_slowest10_rt', 'total_avg_rs', 'total_median_rs', 'total_stdev_rs',
    'fastest10_rs', 'slowest10_rs', 'corr_total_avg_rs', 'corr_total_median_rs',
    'corr_total_stdev_rs', 'corr_fastest10_rs', 'corr_slowest10_rs', 'n_rt',
    'n_rt_corrected', 'false_starts_removed', 'outliers_removed', 'rs_undefined',
    'rs_slope_per_min', 'rt_slope_per_min', 'mode', 'isi_ms', 'isi_set_s', 'block_s',
    'schedule_method', 'schedule_seed', 'scheduled_stimuli', 'stim_ms',
    'hit_window_ms', 'lapse_threshold_ms', 'miss_criterion', 'max_minutes',
    'correction', 'false_start_threshold_ms', 'sd_multiplier', 'kss_when',
    'kss_before', 'kss_after', 'total_presses', 'extra_presses', 'burst_max',
    'rapid_pairs', 'cheating_suspected', 'cheating_reasons', 'alarm_enabled',
    'extra_responses', 'page_blur_count', 'device_browser', 'device_platform',
    'screen_w', 'screen_h', 'device_pixel_ratio', 'refresh_hz_measured',
    'frame_interval_ms', 'onset_quantisation_ms', 'frame_mad_ms',
    'timer_resolution_ms', 'cross_origin_isolated', 'input_dispatch_median_ms',
    'input_dispatch_mad_ms', 'input_dispatch_n', 'input_stamp_usable',
    'rt_source_setting', 'mean_rt_event_ms', 'mean_rt_handler_ms',
    'mean_input_delay_ms', 'median_input_delay_ms', 'n_rt_compared',
    'n_stamp_rejected', 'frames_trial', 'dropped_frames_trial', 'dropped_rate_trial',
    'ran_fullscreen', 'touch_used', 'cpu_cores', 'device_memory_gb', 'page_protocol',
    'page_host', 'page_local', 'page_network_ms', 'connection_type', 'device_grade',
]

NORMS_HEADER = PARTICIPANT_COLS + [
    'norm_available', 'norm_reason', 'norm_source', 'norm_hour', 'norm_window_min',
    'norm_test_min', 'norm_protocol_min', 'norm_below_protocol', 'norm_truncated',
    'norm_ref_sessions', 'norm_ref_n_hour_bin', 'norm_hour_max_min',
    'norm_ref_n_control', 'norm_ref_n_overnight', 'norm_ref_n_ended_early',
    'norm_mixed_studies', 'section', 'variable', 'unit', 'higher_is', 'norm_value',
    'norm_ref_mean', 'norm_ref_sd', 'z_worse', 'band', 'zero_variance_reference',
]


def rnd(v, places):
    """Round for export, leaving None as an empty cell rather than a zero."""
    return None if v is None else round(v, places)


def _flag(b):
    return 1 if b else 0


def participant_vals(rec):
    p = rec['participant']
    return [rec['runId'], p.get('date', ''), p.get('time', ''), rec.get('language', 'en'),
            p.get('participantId', ''), p.get('name', ''), p.get('address', ''),
            p.get('birthDate', ''), p.get('education', ''), p.get('sessionLabel', ''),
            p.get('trialNumber', '')]


def raw_rows(rec):
    pv = participant_vals(rec)
    scored = rec['scored']
    extra = {}
    for p in rec['raw']['presses']:
        k = p.get('epochIndex')
        if k is not None:
            extra[k] = extra.get(k, 0) + 1
    out = []
    for i, t in enumerate(scored['trials']):
        e = rec['raw']['epochs'][i]
        out.append(pv + [
            t['index'] + 1,                       # 1-based, as on screen
            None if t['block'] is None else t['block'] + 1,
            t['minute'] + 1,
            rnd(t['onsetMs'], 3), t['epochIsiMs'], t['isiBeforeMs'],
            _flag(t['rtMs'] is not None),
            # This build measures one reaction time, from the vsync timestamp of
            # the frame that carried the stimulus to the keyboard's own stamp.
            # There is no second clock to correct against, so the browser's
            # event/handler pair is left empty rather than filled with copies.
            rnd(t['rtMs'], 3), None, None, 'psychopy_flip_to_kb',
            rnd(t['rtMs'], 3), rnd(t['rsPerSec'], 6), t['outcome'],
            t['lapse'], t['lateResponse'], t['falseStart'],
            max(0, extra.get(t['index'], 0) - 1),
        ])
    return out


def pm_rows(rec):
    pv = participant_vals(rec)
    scored = rec['scored']
    drs = scored['dynamicsRs']
    drt = scored['dynamicsRt']
    out = []
    for i, m in enumerate(scored['perMinute']):
        out.append(pv + [
            m['minute'], m['trials'], m['hits'], m['misses'], m['lapses'],
            m['lateResponses'], rnd(m['hitRatio'], 6),
            rnd(m['avgRt'], 3), rnd(m['medianRt'], 3), rnd(m['sdRt'], 3),
            rnd(m['fastest10Rt'], 3), rnd(m['slowest10Rt'], 3),
            rnd(m['corrAvgRt'], 3), rnd(m['corrMedianRt'], 3), rnd(m['corrSdRt'], 3),
            rnd(m['corrFastest10Rt'], 3), rnd(m['corrSlowest10Rt'], 3),
            rnd(m['avgRs'], 6), rnd(m['medianRs'], 6), rnd(m['sdRs'], 6),
            rnd(m['fastest10Rs'], 6), rnd(m['slowest10Rs'], 6),
            rnd(m['corrAvgRs'], 6), rnd(m['corrMedianRs'], 6), rnd(m['corrSdRs'], 6),
            rnd(m['corrFastest10Rs'], 6), rnd(m['corrSlowest10Rs'], 6),
            m['n'], m['nCorrected'], m['nFalseStartsRemoved'], m['nOutliersRemoved'],
            rnd(drs['velocity'][i], 6), rnd(drs['acceleration'][i], 6),
            rnd(drt['velocity'][i], 3), rnd(drt['acceleration'][i], 3),
        ])
    return out


def summary_rows(rec):
    """The one-row-per-trial summary.

    Built as a name -> value mapping and then emitted in SUMMARY_HEADER order,
    rather than as a 116-long positional list. A positional list silently
    shifts every column after a missing value — which is exactly what happened
    when this was first written, putting the RT source into the
    input_stamp_usable column and making the frame counts read as zero. Here a
    missing column is a KeyError and an unknown one is caught immediately.
    """
    pv = participant_vals(rec)
    s = rec['scored']
    t = s['totals']
    ep = s['errorProfiles']
    g = s['integrity']
    cfg = rec['config']
    raw = rec['raw']

    correction = 'none'
    if cfg['removeFalseStarts'] and cfg['removeOutliers']:
        correction = 'false_starts+outliers'
    elif cfg['removeFalseStarts']:
        correction = 'false_starts'
    elif cfg['removeOutliers']:
        correction = 'outliers'

    extra_responses = sum(r[-1] for r in raw_rows(rec)) if s['trials'] else 0
    frames = raw['frames']

    v = {
        'sleep_onset_ms': rnd(raw.get('sleepOnsetMs'), 3),
        'sleep_onset_criterion_ms': cfg['missCriterion'] * cfg['isiMs'],
        'slept_before_max': _flag(raw['sleptBeforeMax']),
        'end_reason': raw['endReason'],
        'elapsed_ms': rnd(raw['elapsedMs'], 3),
        'total_trialrun': t['trials'],
        'hit_ratio': rnd(t['hitRatio'], 6),
        'total_hits': t['hits'],
        'total_miss': t['misses'],
        'total_lapse': t['lapses'],
        'late_responses': t['lateResponses'],
        'total_false_starts': t['falseStarts'],
        'ep_1_2': ep['ep1_2'], 'ep_3_6': ep['ep3_6'], 'ep_7plus': ep['ep7plus'],
        'longest_miss_run': ep['longestRun'],
        'total_avg_rt': rnd(t['avgRt'], 3),
        'total_median_rt': rnd(t['medianRt'], 3),
        'total_stdev_rt': rnd(t['sdRt'], 3),
        'fastest10_rt': rnd(t['fastest10Rt'], 3),
        'slowest10_rt': rnd(t['slowest10Rt'], 3),
        'corr_total_avg_rt': rnd(t['corrAvgRt'], 3),
        'corr_total_median_rt': rnd(t['corrMedianRt'], 3),
        'corr_total_stdev_rt': rnd(t['corrSdRt'], 3),
        'corr_fastest10_rt': rnd(t['corrFastest10Rt'], 3),
        'corr_slowest10_rt': rnd(t['corrSlowest10Rt'], 3),
        'total_avg_rs': rnd(t['avgRs'], 6),
        'total_median_rs': rnd(t['medianRs'], 6),
        'total_stdev_rs': rnd(t['sdRs'], 6),
        'fastest10_rs': rnd(t['fastest10Rs'], 6),
        'slowest10_rs': rnd(t['slowest10Rs'], 6),
        'corr_total_avg_rs': rnd(t['corrAvgRs'], 6),
        'corr_total_median_rs': rnd(t['corrMedianRs'], 6),
        'corr_total_stdev_rs': rnd(t['corrSdRs'], 6),
        'corr_fastest10_rs': rnd(t['corrFastest10Rs'], 6),
        'corr_slowest10_rs': rnd(t['corrSlowest10Rs'], 6),
        'n_rt': t['n'], 'n_rt_corrected': t['nCorrected'],
        'false_starts_removed': t['nFalseStartsRemoved'],
        'outliers_removed': t['nOutliersRemoved'],
        'rs_undefined': t['nRsUndefined'],
        'rs_slope_per_min': rnd(s['dynamicsRs']['slope'], 6),
        'rt_slope_per_min': rnd(s['dynamicsRt']['slope'], 3),

        'mode': cfg.get('mode', 'bsrt'),
        'isi_ms': cfg['isiMs'],
        'isi_set_s': ('/'.join(str(x / 1000.0) for x in cfg.get('isiSetMs', []))
                      if cfg.get('mode') == 'pvt' else ''),
        'block_s': cfg.get('blockMs', 30000) / 1000.0,
        'schedule_method': rec['schedule']['method'],
        'schedule_seed': rec['schedule']['seed'],
        'scheduled_stimuli': rec['schedule']['nStimuli'],
        'stim_ms': cfg['stimMs'],
        'hit_window_ms': cfg['hitWindowMs'],
        'lapse_threshold_ms': cfg['lapseMs'],
        'miss_criterion': cfg['missCriterion'],
        'max_minutes': cfg['maxMinutes'],
        'correction': correction,
        'false_start_threshold_ms': cfg.get('falseStartMs', 100),
        'sd_multiplier': cfg.get('sdMultiplier', 2),

        'kss_when': rec.get('kssWhen', 'none'),
        'kss_before': rec.get('kssBefore', ''),
        'kss_after': rec.get('kssAfter', ''),

        'total_presses': g['totalPresses'],
        'extra_presses': g['extraPresses'],
        'burst_max': g['burstMax'],
        'rapid_pairs': g['rapidPairs'],
        'cheating_suspected': _flag(g['suspected']),
        'cheating_reasons': '; '.join(g['reasons']),
        'alarm_enabled': _flag(cfg.get('alarm', True)),
        'extra_responses': extra_responses,

        # What this build can actually measure about its own timing.
        'page_blur_count': 0,
        'device_platform': 'psychopy',
        'refresh_hz_measured': rnd(raw.get('refreshHz'), 3),
        'frame_interval_ms': rnd(raw['frameIntervalMs'], 4),
        # On a frame-locked build the quantisation of an onset IS the frame.
        'onset_quantisation_ms': rnd(raw['frameIntervalMs'], 4),
        'frame_mad_ms': rnd(raw.get('frameMadMs'), 4),
        'rt_source_setting': 'psychopy_flip_to_kb',
        'frames_trial': frames,
        'dropped_frames_trial': raw['droppedFrames'],
        'dropped_rate_trial': rnd(raw['droppedFrames'] / frames, 6) if frames else '',
        'ran_fullscreen': _flag(cfg.get('fullscreen', True)),
        'page_local': 1,
    }

    # Everything the browser build measures about a browser is left empty here.
    # The column stays so the files still rbind; empty means "this build cannot
    # know", which is not the same as zero.
    for col in SUMMARY_HEADER:
        v.setdefault(col, '')

    unknown = set(v) - set(SUMMARY_HEADER)
    if unknown:
        raise KeyError('summary values with no column: %s' % sorted(unknown))

    return [pv + [v[c] for c in SUMMARY_HEADER[len(PARTICIPANT_COLS):]]]


def norms_rows(rec):
    pv = participant_vals(rec)
    nm = rec.get('norms')
    if not nm or not nm['available']:
        return [pv + [0, nm['reason'] if nm else 'no_norms'] + [''] * 24]
    prov = nm['provenance'] or [None, None, None]
    ctx = [1, '', nm['source'], nm['hour'], nm['windowMinutes'], nm['testMinutes'],
           '/'.join(str(p) for p in nm['protocols']) or nm['protocolMinutes'],
           _flag(nm['belowProtocol']), _flag(nm['truncated']),
           nm['sessions'], nm['n'], nm['hourMaxLength'],
           prov[0], prov[1], prov[2], _flag(nm['mixedStudies'])]
    out = []
    for row in nm['rows']:
        c = row['comparison']
        out.append(pv + ctx + [
            row['section'], row['label'], row['unit'],
            'not judged' if row['dir'] == 0 else ('better' if row['dir'] > 0 else 'worse'),
            rnd(c['value'], 5) if c else None,
            rnd(c['mean'], 5) if c else None,
            rnd(c['sd'], 5) if c else None,
            rnd(c['z'], 4) if c and c['z'] is not None else None,
            c['band'] if c else '',
            _flag(c['degenerate']) if c else 0,
        ])
    return out


def to_csv(header, rows):
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator='\n')
    w.writerow(header)
    for r in rows:
        w.writerow(['' if v is None else v for v in r])
    return buf.getvalue()


EXPORTS = (
    ('raw', RAW_HEADER, raw_rows),
    ('perminute', PM_HEADER, pm_rows),
    ('summary', SUMMARY_HEADER, summary_rows),
    ('norms', NORMS_HEADER, norms_rows),
)


def write_all(rec, out_dir):
    """Write all four CSVs for one trial. Returns the paths written.

    Four separate files rather than one, matching the other builds, because the
    raw file is the only irreplaceable one and it should be obvious which that
    is.
    """
    stem = file_stem(rec)
    os.makedirs(out_dir, exist_ok=True)
    paths = []
    for name, header, fn in EXPORTS:
        path = os.path.join(out_dir, '%s_%s.csv' % (stem, name))
        with open(path, 'w', newline='') as f:
            f.write(to_csv(header, fn(rec)))
        paths.append(path)
    return paths


def file_stem(rec):
    p = rec['participant']
    return 'bsrt_%s_%s_t%s' % (safe(p.get('participantId') or 'NA'),
                               safe(p.get('sessionLabel') or 'NA'),
                               p.get('trialNumber', 1))


def safe(s):
    return ''.join(c if c.isalnum() or c in '-_' else '_' for c in str(s))
