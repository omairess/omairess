"""Build BSRT/norms.js from the raw trial data of both normative studies.

Run:  python3 gen_norms2.py [--write]

Without --write it prints the tables and quality checks and touches nothing.
"""

import os, sys, json, math, collections, statistics as st
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from normlib import read_sessions, summarise, VAR_KEYS

# Where the two raw-trial workbooks live. They carry dates of birth, so they are
# deliberately NOT committed to this repository — point this at your own copy.
U = os.environ.get('BSRT_NORM_DIR', os.path.expanduser('~/bsrt-norm-source')) + '/'
CTRL_FILE = 'BSRT_CTRL_GROUP.xlsx'
IH_FILE = 'BSRT_outcome_calculations_IH_2024.xlsx'

# The overnight study recorded a session number, not a clock time. The four
# sessions ran at fixed hours, given by the study.
IH_SESSION_HOUR = {1: 23, 2: 1, 3: 3, 4: 5}

MAX_LEN = 40
TRIALS_PER_MIN = 20
ROUND = 4

VAR_META = [
    ('trials', 'Total trials run', 'Test overview', 'count', 0),
    ('hitRatio', 'Hit ratio', 'Test overview', '%', 1),
    ('hits', 'Hits', 'Test overview', 'count', 1),
    ('misses', 'Misses', 'Test overview', 'count', -1),
    ('lapses', 'Lapses (>500 ms)', 'Test overview', 'count', -1),
    ('falseStarts', 'False starts (<100 ms)', 'Test overview', 'count', -1),
    ('ep12', 'Error profile EP1–2', 'Test overview', 'runs', -1),
    ('ep36', 'Error profile EP3–6', 'Test overview', 'runs', -1),
    ('ep7', 'Error profile EP7+', 'Test overview', 'runs', -1),
    ('rtMean', 'RT average', 'Reaction time', 'ms', -1),
    ('rtMedian', 'RT median', 'Reaction time', 'ms', -1),
    ('rtSd', 'RT within-test SD', 'Reaction time', 'ms', -1),
    ('rtFast10', 'RT 10% fastest (mean)', 'Reaction time', 'ms', -1),
    ('rtSlow10', 'RT 10% slowest (mean)', 'Reaction time', 'ms', -1),
    ('rtIpr', 'RT interpercentile range', 'Reaction time', 'ms', -1),
    ('rsMean', 'RS average', 'Reaction speed', '1000/RT', 1),
    ('rsMedian', 'RS median', 'Reaction speed', '1000/RT', 1),
    ('rsSd', 'RS within-test SD', 'Reaction speed', '1000/RT', -1),
    ('rsFast10', 'RS 10% fastest (mean)', 'Reaction speed', '1000/RT', 1),
    ('rsSlow10', 'RS 10% slowest (mean)', 'Reaction speed', '1000/RT', 1),
    ('rsIpr', 'RS interpercentile range', 'Reaction speed', '1000/RT', -1),
]
assert [v[0] for v in VAR_META] == VAR_KEYS


def load():
    ctrl = read_sessions(U + CTRL_FILE, 'Feuil1', 7, 807,
                         lambda r: r[2].hour)
    for s in ctrl:
        s['study'] = 'ctrl'
        s['protocol'] = 8
    ih = read_sessions(U + IH_FILE, 'Sheet1', 2, 802,
                       lambda r: IH_SESSION_HOUR[r[1]])
    for s in ih:
        s['study'] = 'ih'
        s['protocol'] = 40
        s['id'] = 'IH%s' % s['id']
    return ctrl, ih


def build(sessions):
    """stats[hour][length-1] -> flat [mean, sd] * 21, or None if the cell is empty."""
    by_hour = collections.defaultdict(list)
    for s in sessions:
        by_hour[s['hour']].append(s)

    stats, ns, prov = [], [], []
    for h in range(24):
        row_s, row_n, row_p = [], [], []
        for m in range(1, MAX_LEN + 1):
            got = [(s, summarise(s['raw'], m, TRIALS_PER_MIN)) for s in by_hour.get(h, [])]
            used = [(s, x) for s, x in got if x]
            # A session that stopped at sleep onset before this length is
            # counted separately: it is missing because the participant fell
            # asleep, which is not the same as never having been tested.
            ended = sum(1 for s, x in got
                        if x is None and len(s['raw']) < s['protocol'] * TRIALS_PER_MIN)
            if not used:
                row_s.append(None); row_n.append(None); row_p.append(None)
                continue
            flat = []
            for k in VAR_KEYS:
                vals = [x[k] for _, x in used if x.get(k) is not None]
                mean = round(st.mean(vals), ROUND) if vals else None
                sd = round(st.stdev(vals), ROUND) if len(vals) > 1 else (0.0 if vals else None)
                flat.append(mean); flat.append(sd)
            row_s.append(flat)
            row_n.append(len(used))
            row_p.append([sum(1 for s, _ in used if s['study'] == 'ctrl'),
                          sum(1 for s, _ in used if s['study'] == 'ih'),
                          ended])
        stats.append(row_s); ns.append(row_n); prov.append(row_p)
    return stats, ns, prov


def main():
    ctrl, ih = load()
    allses = ctrl + ih
    print('CTRL %d sessions, %d participants (8-min protocol)'
          % (len(ctrl), len({s['id'] for s in ctrl})))
    print('IH   %d sessions, %d participants (40-min protocol)'
          % (len(ih), len({s['id'] for s in ih})))

    stats, ns, prov = build(allses)

    print('\nsessions per hour bin (ctrl + ih), and how deep the bin goes:')
    for h in range(24):
        lens = [m + 1 for m in range(MAX_LEN) if ns[h][m]]
        if not lens:
            continue
        p = prov[h][0]
        print('  %02d:00  n=%-3d (ctrl %2d, ih %2d)  lengths 1-%d'
              % (h, ns[h][0], p[0], p[1], max(lens)))

    print('\ncensoring at the long lengths (overnight study only):')
    for m in (8, 10, 20, 30, 40):
        tot = sum(ns[h][m - 1] or 0 for h in range(24))
        end = sum((prov[h][m - 1] or [0, 0, 0])[2] for h in range(24))
        print('  %2d min: %3d sessions contribute, %3d already ended at sleep onset' % (m, tot, end))

    print('\nthe two studies where they overlap (hours 23, 01, 03, 05; RT average, ms):')
    print('  %-6s %-22s %-22s' % ('length', 'control group', 'overnight study'))
    ci = VAR_KEYS.index('rtMean')
    for m in (1, 3, 5, 8):
        a = [x['rtMean'] for s in ctrl if s['hour'] in (23, 1, 3, 5)
             for x in [summarise(s['raw'], m)] if x and x.get('rtMean')]
        b = [x['rtMean'] for s in ih
             for x in [summarise(s['raw'], m)] if x and x.get('rtMean')]
        print('  %-6d %7.1f (SD %5.1f) n=%-3d %7.1f (SD %5.1f) n=%-3d'
              % (m, st.mean(a), st.stdev(a), len(a), st.mean(b), st.stdev(b), len(b)))

    print('\nmiss rate, same hours (%% of trials):')
    for m in (3, 8):
        a = [x['misses'] / x['trials'] * 100 for s in ctrl if s['hour'] in (23, 1, 3, 5)
             for x in [summarise(s['raw'], m)] if x]
        b = [x['misses'] / x['trials'] * 100 for s in ih
             for x in [summarise(s['raw'], m)] if x]
        print('  %2d min: control %.2f%%   overnight %.2f%%' % (m, st.mean(a), st.mean(b)))

    data = {
        'source': '%s (control group, 8-min protocol) + %s (overnight study, 40-min protocol)'
                  % (CTRL_FILE, IH_FILE),
        'sessions': len(allses),
        'participants': len({s['id'] for s in allses}),
        'studies': [
            {'key': 'ctrl', 'label': 'Control group', 'file': CTRL_FILE,
             'sessions': len(ctrl), 'participants': len({s['id'] for s in ctrl}),
             'protocolMinutes': 8, 'hours': 'all 24'},
            {'key': 'ih', 'label': 'Overnight study', 'file': IH_FILE,
             'sessions': len(ih), 'participants': len({s['id'] for s in ih}),
             'protocolMinutes': 40, 'hours': '23, 01, 03, 05'},
        ],
        'protocolMinutes': 8,
        'trialsPerMinute': TRIALS_PER_MIN,
        'hours': 24,
        'maxLength': MAX_LEN,
        'vars': [{'key': k, 'label': l, 'section': s, 'unit': u, 'dir': d}
                 for k, l, s, u, d in VAR_META],
        'n': ns,
        'prov': prov,
        'stats': stats,
    }

    if '--write' in sys.argv:
        emit(data)
        print('\nwritten')
    else:
        print('\n(dry run — pass --write to emit norms.js)')


HEADER = '''\
'use strict';

/* BSRT preliminary normative reference — two studies, computed from raw trials.
 *
 * SHARED FILE — a byte-identical copy lives at BSRT-desktop/renderer/norms.js.
 * Generated by gen_norms2.py from the two raw-trial workbooks; do not hand-edit.
 *
 * WHAT THIS IS. %(sessions)d test sessions from %(participants)d people, binned by the HOUR THE
 * SESSION STARTED (24 bins) and by CUMULATIVE TEST LENGTH, where length m is
 * the first 20*m trial opportunities of the same session. A 3-minute test is
 * compared against the first 3 minutes of the reference sessions, never
 * against their whole run.
 *
 * TWO STUDIES, pooled:
 *   control group     %(nctrl)3d sessions, %(pctrl)2d participants, 8-minute protocol, all 24 hours
 *   overnight study   %(nih)3d sessions, %(pih)2d participants, 40-minute protocol, hours 23, 01, 03, 05
 *
 * So hours 23, 01, 03 and 05 pool both studies for lengths 1-8 and carry the
 * overnight study alone from 9 to 40 minutes; every other hour stops at 8.
 * `prov[hour][length-1]` is [control sessions, overnight sessions, sessions
 * that had already ended at sleep onset], so any cell's composition is
 * recoverable and the pooling is reversible.
 *
 * WHAT IT IS NOT. These are preliminary, and the two samples are not
 * interchangeable: the overnight study misses about %(missih).0f%% of trials against the
 * control group's %(missctrl).1f%% over the same hours, and %(ended)d of its %(nih)d sessions reached the
 * sleep-onset criterion and stopped early. Where a cell pools both, it pools
 * two different populations tested under two different expected durations.
 * Treat a z-score here as a rough position, not a clinical cut-off.
 *
 * CENSORING. A session contributes to every length it completed and to none
 * after it. Sessions that stopped at sleep onset therefore drop out of the
 * longer lengths, so the 20- and 40-minute cells describe the people who
 * stayed awake. That is a survivor bias in the direction of alertness: it
 * makes a sleepy participant look worse than the true population would. The
 * count of dropped-out sessions is the third entry of each `prov` cell.
 *
 * CONVENTIONS, identical to normativeSummary() in scoring.js, so a
 * participant and the reference are scored by the same rule:
 *   miss       no response before the next stimulus (recorded as 3000)
 *   false start  RT < 100 ms; counted, then excluded from RT statistics
 *   valid RT   any response over 100 ms, whatever its latency
 *   lapse      valid RT > 500 ms
 *   decile     mean of the fastest/slowest ceil(10%% x valid trials)
 *   IPR        |slowest decile mean - fastest decile mean|
 *   SD         sample SD across the sessions in the bin
 *
 * VERIFIED. Recomputing the previous 24x8 published tables from the control
 * group's raw trials reproduced all 8064 cells to within 5e-5, the rounding
 * of the published file. The overnight study's own derived columns were
 * reproduced from its raw trials for all six of its test lengths.
 *
 * stats[hour][length-1] is a flat array of mean, sd pairs in `vars` order, or
 * null where that hour was never tested at that length.
 */

(function (root, factory) {
  var api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.BSRTNorms = api;
})(typeof self !== 'undefined' ? self : this, function () {

  var DATA = %(json)s;

  function varIndex(key) {
    for (var v = 0; v < DATA.vars.length; v++) {
      if (DATA.vars[v].key === key) return v;
    }
    return -1;
  }

  function cell(hour, lengthMin) {
    if (!(hour >= 0 && hour < 24)) return null;
    if (!(lengthMin >= 1 && lengthMin <= DATA.maxLength)) return null;
    return DATA.stats[hour][lengthMin - 1];
  }

  /* Mean and SD for one variable, or null when the bin or variable is unknown. */
  function lookup(hour, lengthMin, key) {
    var flat = cell(hour, lengthMin);
    if (!flat) return null;
    var i = varIndex(key);
    if (i < 0) return null;
    return {
      mean: flat[i * 2],
      sd: flat[i * 2 + 1],
      n: DATA.n[hour][lengthMin - 1],
      dir: DATA.vars[i].dir
    };
  }

  /* The longest length this hour was actually tested at, or 0 if never. */
  function maxLengthFor(hour) {
    if (!(hour >= 0 && hour < 24)) return 0;
    var row = DATA.stats[hour];
    for (var m = row.length; m > 0; m--) {
      if (row[m - 1]) return m;
    }
    return 0;
  }

  return {
    data: DATA,
    vars: DATA.vars,
    maxLength: DATA.maxLength,
    protocolMinutes: DATA.protocolMinutes,
    sessions: DATA.sessions,
    studies: DATA.studies,
    lookup: lookup,
    maxLengthFor: maxLengthFor,
    /* [control sessions, overnight sessions, sessions already ended] or null. */
    provenance: function (hour, lengthMin) {
      if (!(hour >= 0 && hour < 24)) return null;
      if (!(lengthMin >= 1 && lengthMin <= DATA.maxLength)) return null;
      return DATA.prov[hour][lengthMin - 1];
    },
    nFor: function (hour, lengthMin) {
      if (!(hour >= 0 && hour < 24)) return null;
      if (!(lengthMin >= 1 && lengthMin <= DATA.maxLength)) return null;
      return DATA.n[hour][lengthMin - 1];
    }
  };
});
'''


def emit(data):
    ctrl_m = data['studies'][0]
    ih_m = data['studies'][1]
    ended = sum((p[2] if p else 0) for h in range(24) for p in [data['prov'][h][MAX_LEN - 1]])
    # headline miss rates for the header, over the shared hours
    body = HEADER % {
        'sessions': data['sessions'], 'participants': data['participants'],
        'nctrl': ctrl_m['sessions'], 'pctrl': ctrl_m['participants'],
        'nih': ih_m['sessions'], 'pih': ih_m['participants'],
        'missih': MISS_IH, 'missctrl': MISS_CTRL, 'ended': ended,
        'json': json.dumps(data, separators=(',', ':'), ensure_ascii=False),
    }
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    for p in (os.path.join(here, 'BSRT', 'norms.js'),
              os.path.join(here, 'BSRT-desktop', 'renderer', 'norms.js')):
        with open(p, 'w') as f:
            f.write(body)
        print('  wrote', p)


MISS_CTRL = 0.0
MISS_IH = 0.0

if __name__ == '__main__':
    ctrl, ih = load()
    MISS_CTRL = st.mean([x['misses'] / x['trials'] * 100 for s in ctrl if s['hour'] in (23, 1, 3, 5)
                         for x in [summarise(s['raw'], 8)] if x])
    MISS_IH = st.mean([x['misses'] / x['trials'] * 100 for s in ih
                       for x in [summarise(s['raw'], 8)] if x])
    main()
