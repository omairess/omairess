"""The PsychoPy CSVs must have the same columns as the browser build's.

An analysis script should not care which build wrote a file. This extracts the
header lists straight out of BSRT/app.js and compares them, in order, with the
ones in bsrt_io.py — so adding a column to one build fails here until it is
added to the other, rather than being discovered months later in R.
"""

import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.dirname(HERE))

import bsrt_io as io_  # noqa: E402

APP_JS = os.path.join(ROOT, 'BSRT', 'app.js')

EXTRACT = r"""
const fs = require('fs');
const src = fs.readFileSync(process.argv[2], 'utf8');
function grab(name) {
  const start = src.indexOf('var ' + name + ' = ');
  if (start < 0) throw new Error('missing ' + name);
  let i = src.indexOf('=', start) + 1, depth = 0, out = '';
  for (; i < src.length; i++) {
    const c = src[i];
    if ('([{'.includes(c)) depth++;
    if (')]}'.includes(c)) depth--;
    if (c === ';' && depth === 0) break;
    out += c;
  }
  return out;
}
const names = ['PARTICIPANT_COLS', 'RAW_HEADER', 'PM_HEADER', 'SUMMARY_HEADER', 'NORMS_HEADER'];
let code = '';
for (const n of names) code += 'var ' + n + ' = ' + grab(n) + ';\n';
code += 'JSON.stringify({' + names.map(n => n + ':' + n).join(',') + '})';
console.log(eval(code));
"""

failures = []


def check(cond, msg):
    print(('  ok    ' if cond else '  FAIL  ') + msg)
    if not cond:
        failures.append(msg)


def node_bin():
    return '/opt/node22/bin/node' if os.path.exists('/opt/node22/bin/node') else 'node'


with tempfile.NamedTemporaryFile('w', suffix='.js', delete=False) as f:
    f.write(EXTRACT)
    extractor = f.name
try:
    proc = subprocess.run([node_bin(), extractor, APP_JS], capture_output=True, text=True)
    if proc.returncode != 0:
        print('could not read headers from app.js:\n' + proc.stderr[:800])
        sys.exit(1)
    web = json.loads(proc.stdout)
finally:
    os.unlink(extractor)

print('=== column names, against BSRT/app.js ===')
mine = {
    'PARTICIPANT_COLS': io_.PARTICIPANT_COLS,
    'RAW_HEADER': io_.RAW_HEADER,
    'PM_HEADER': io_.PM_HEADER,
    'SUMMARY_HEADER': io_.SUMMARY_HEADER,
    'NORMS_HEADER': io_.NORMS_HEADER,
}
for name in ('PARTICIPANT_COLS', 'RAW_HEADER', 'PM_HEADER', 'SUMMARY_HEADER', 'NORMS_HEADER'):
    a, b = mine[name], web[name]
    if a == b:
        check(True, '%-17s %3d columns, identical' % (name, len(a)))
        continue
    only_py = [c for c in a if c not in b]
    only_js = [c for c in b if c not in a]
    detail = []
    if len(a) != len(b):
        detail.append('%d vs %d columns' % (len(a), len(b)))
    if only_py:
        detail.append('only in psychopy: ' + ', '.join(only_py[:6]))
    if only_js:
        detail.append('only in web: ' + ', '.join(only_js[:6]))
    if not detail:
        first = next(i for i, (x, y) in enumerate(zip(a, b)) if x != y)
        detail.append('order differs at %d: %s vs %s' % (first, a[first], b[first]))
    check(False, '%-17s %s' % (name, '; '.join(detail)))

print()
if failures:
    print('%d HEADER CHECK(S) FAILED' % len(failures))
    sys.exit(1)
print('ALL FOUR EXPORTS MATCH THE BROWSER BUILD COLUMN FOR COLUMN')
