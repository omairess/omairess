"""Extract the on-screen strings the PsychoPy build needs from BSRT/i18n.js.

Run:  python3 tools/gen_strings.py

The PsychoPy build shows a small subset of the interface: instructions, the
countdown hint, the KSS, the alarm banner and the end-of-trial message. Those
strings already exist, in four languages, in i18n.js — retyping them here would
let the two drift, and the KSS anchors in particular belong to a validated
instrument and must be identical wherever they are shown.

So they are pulled straight out of i18n.js by running it under node, and written
to BSRT-psychopy/strings.json. Add a key to KEYS below and re-run.
"""

import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
I18N = os.path.join(ROOT, 'BSRT', 'i18n.js')
OUT = os.path.join(ROOT, 'BSRT-psychopy', 'strings.json')

LANGS = ['en', 'fr', 'nl', 'de']

KEYS = [
    'kss.heading', 'kss.instruction', 'kss.question', 'kss.continue',
    'kss.beforeTitle', 'kss.afterTitle',
    'task.hintBsrt', 'task.hintPvt', 'countdown.hint',
]


def node_bin():
    return '/opt/node22/bin/node' if os.path.exists('/opt/node22/bin/node') else 'node'


SCRIPT = """
const L = require(%s);
const langs = %s, keys = %s;
const out = {};
for (const lang of langs) {
  L.setLanguage(lang);
  const d = { _anchors: L.kssAnchors(lang) };
  for (const k of keys) d[k] = L.t(k);
  out[lang] = d;
}
console.log(JSON.stringify(out));
"""


def main():
    script = SCRIPT % (json.dumps(I18N), json.dumps(LANGS), json.dumps(KEYS))
    proc = subprocess.run([node_bin(), '-e', script], capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        return 1
    data = json.loads(proc.stdout)

    missing = []
    for lang in LANGS:
        for k in KEYS:
            v = data[lang].get(k)
            # i18n.t() falls back to English rather than returning the key, so a
            # missing translation shows up as an identical string, not a crash.
            if not v or v == k:
                missing.append('%s/%s' % (lang, k))
        if len(data[lang]['_anchors']) != 9:
            missing.append('%s: %d KSS anchors, expected 9'
                           % (lang, len(data[lang]['_anchors'])))
    if missing:
        sys.stderr.write('missing strings: %s\n' % ', '.join(missing))
        return 1

    with open(OUT, 'w') as f:
        json.dump(data, f, ensure_ascii=False, indent=1, sort_keys=True)
    print('wrote %s (%d languages, %d keys plus the 9 KSS anchors)'
          % (OUT, len(LANGS), len(KEYS)))
    for lang in LANGS:
        print('  %s  %s' % (lang, data[lang]['_anchors'][0]))
    return 0


if __name__ == '__main__':
    sys.exit(main())
