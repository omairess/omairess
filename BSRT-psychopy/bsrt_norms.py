"""The normative tables, and the comparison against them.

Port of norms.js and normativeReport() from scoring.js. The tables themselves
are `norms.json`, written by tools/gen_norms2.py in the same run that writes
norms.js — one generated artefact, three builds, so a z-score does not depend on
which front end collected the data.

See BSRT/README.md for what the reference is and, more importantly, what it is
not: two pooled studies of unequal sleepiness, a survivor bias at the longer
lengths, and hour bins holding a dozen sessions each.
"""

import json
import os

import bsrt_core as core

_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_PATH = os.path.join(_HERE, 'norms.json')


class Norms(object):
    def __init__(self, data):
        self.data = data
        self.vars = data['vars']
        self.max_length = data['maxLength']
        self.sessions = data['sessions']
        self.protocol_minutes = data['protocolMinutes']
        self._index = {v['key']: i for i, v in enumerate(self.vars)}

    @classmethod
    def load(cls, path=None):
        with open(path or DEFAULT_PATH) as f:
            return cls(json.load(f))

    def cell(self, hour, length_min):
        if not (0 <= hour < 24):
            return None
        if not (1 <= length_min <= self.max_length):
            return None
        return self.data['stats'][hour][length_min - 1]

    def lookup(self, hour, length_min, key):
        """Mean and SD for one variable, or None when the bin is empty."""
        flat = self.cell(hour, length_min)
        if flat is None:
            return None
        i = self._index.get(key)
        if i is None:
            return None
        return {'mean': flat[i * 2], 'sd': flat[i * 2 + 1],
                'n': self.data['n'][hour][length_min - 1], 'dir': self.vars[i]['dir']}

    def max_length_for(self, hour):
        """The longest length this hour was actually tested at, or 0 if never.

        Only the four overnight hours go past 8 minutes, so this is what decides
        the comparison window rather than the table's global maximum.
        """
        if not (0 <= hour < 24):
            return 0
        row = self.data['stats'][hour]
        for m in range(len(row), 0, -1):
            if row[m - 1]:
                return m
        return 0

    def n_for(self, hour, length_min):
        if not (0 <= hour < 24) or not (1 <= length_min <= self.max_length):
            return None
        return self.data['n'][hour][length_min - 1]

    def provenance(self, hour, length_min):
        """[control sessions, overnight sessions, sessions already ended]."""
        if not (0 <= hour < 24) or not (1 <= length_min <= self.max_length):
            return None
        return self.data['prov'][hour][length_min - 1]


def normative_report(o, norms):
    """The whole normative comparison for one trial, or a refusal with a reason.

    Eligibility is deliberately strict. The reference sessions are BSRT runs at
    a 3000 ms interval, which is 20 stimuli a minute; a PVT averages ten, and a
    BSRT at another interval delivers a different number again. Rather than emit
    a plausible-looking z-score against the wrong reference, the report says
    which condition failed.
    """
    out = {
        'available': False, 'reason': None, 'hour': o.get('hour'), 'mode': o.get('mode'),
        'isiMs': o.get('isiMs'), 'windowMinutes': None, 'testMinutes': None,
        'truncated': False, 'belowProtocol': False,
        'protocolMinutes': norms.protocol_minutes if norms else None,
        'hourMaxLength': None, 'provenance': None, 'protocols': [],
        'mixedStudies': False, 'censored': 0,
        'sessions': norms.sessions if norms else None,
        'participants': norms.data['participants'] if norms else None,
        'source': norms.data['source'] if norms else None,
        'n': None, 'summary': None, 'rows': [],
    }
    if not norms:
        out['reason'] = 'no_norms'
        return out
    if o.get('mode') != 'bsrt':
        out['reason'] = 'not_bsrt'
        return out
    if o.get('isiMs') != 3000:
        out['reason'] = 'isi_mismatch'
        return out
    hour = o.get('hour')
    if hour is None or not (0 <= hour < 24):
        out['reason'] = 'no_hour'
        return out

    hour_max = norms.max_length_for(hour)
    if not hour_max:
        out['reason'] = 'no_hour'
        return out
    out['hourMaxLength'] = hour_max

    completed = int((o.get('elapsedMs') or 0) // 60000)
    buckets = completed if o.get('minuteBuckets') is None else o['minuteBuckets']
    window = min(completed, buckets, hour_max)
    out['testMinutes'] = completed
    if window < 1:
        out['reason'] = 'too_short'
        return out
    out['windowMinutes'] = window
    out['truncated'] = completed > window

    prov = norms.provenance(hour, window)
    if prov:
        out['provenance'] = prov
        out['censored'] = prov[2] or 0
        if prov[0]:
            out['protocols'].append(8)
        if prov[1]:
            out['protocols'].append(40)
        out['mixedStudies'] = prov[0] > 0 and prov[1] > 0
    out['belowProtocol'] = any(window < p for p in out['protocols'])

    summary = core.normative_summary(o['trials'], window, o)
    if not summary:
        out['reason'] = 'too_short'
        return out

    out['summary'] = summary
    out['n'] = norms.n_for(hour, window)
    out['available'] = True

    for v in norms.vars:
        norm = norms.lookup(hour, window, v['key'])
        cmp_ = core.compare_to_norm(summary.get(v['key']), norm)
        out['rows'].append({'key': v['key'], 'label': v['label'], 'section': v['section'],
                            'unit': v['unit'], 'dir': v['dir'], 'comparison': cmp_})

    out['nOrange'] = sum(1 for r in out['rows']
                         if r['comparison'] and r['comparison']['band'] == 'orange')
    out['nRed'] = sum(1 for r in out['rows']
                      if r['comparison'] and r['comparison']['band'] == 'red')
    return out
