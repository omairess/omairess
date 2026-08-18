"""Run every BSRT-psychopy test suite and report a single verdict.

    python3 BSRT-psychopy/tests/run_all.py
"""

import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

SUITES = [
    ('test_equivalence.py', 'the Python scoring against scoring.js'),
    ('test_task.py', 'the trial loop on a virtual clock'),
    ('test_headers.py', 'CSV columns against the browser build'),
    ('test_export.py', 'a whole trial, written to four files'),
    ('test_psychopy_path.py', 'the PsychoPy branch, against a stub'),
    ('test_check_setup.py', 'the setup checker, and its verdict rule'),
]


def main():
    failed = []
    for name, what in SUITES:
        print('=' * 70)
        print('%s  —  %s' % (name, what))
        print('=' * 70)
        rc = subprocess.call([sys.executable, os.path.join(HERE, name)])
        if rc != 0:
            failed.append(name)
        print()

    print('=' * 70)
    if failed:
        print('FAILED: %s' % ', '.join(failed))
        return 1
    print('ALL %d SUITES PASSED' % len(SUITES))
    print()
    print('Note: none of this proves the build works against a REAL PsychoPy')
    print('installation — see "What is verified, and what is not" in README.md.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
