#!/usr/bin/env python3
"""List initialised locals of the parallel routines that may be read before they
are written in a call.

usage: uninit_scan.py SRC_DIR MOD_DIR [--allowed allowed.txt] [--fc gfortran]

omp_generate.py gives the OpenMP build a copy of these locals that is set at the
start of every call. That changes results when the serial code reads such a local
before writing it and so uses the value from the previous call (the previous plant,
hru or object). This compiles each routine with the initialisers removed and no
reset, and lets gfortran's -Wmaybe-uninitialized point at the reads. MOD_DIR contains
the .mod files of a build of SRC_DIR (any build directory).

Candidates already reviewed are listed in allowed.txt as `local` (harmless) or
`saved` (kept as they are); the rest are printed. gfortran does not see every path,
so a clean scan does not replace the comparison runs in check.sh.
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ompaudit
import omp_generate


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('src')
    ap.add_argument('mods')
    ap.add_argument('--allowed', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'allowed.txt'))
    ap.add_argument('--fc', default='gfortran')
    a = ap.parse_args()

    allowed = ompaudit.read_allowed(a.allowed)
    known = set(allowed['saved']) | set(allowed['local'])
    res, tree = ompaudit.analyse(a.src, omp_generate.ROOTS)
    units, byname, modules = ompaudit.build(a.src)

    # files and the initialised, written locals of each routine in them
    cand = {}
    for n, r in res.items():
        u = byname[n]
        written = omp_generate.written_names(u)
        for s in r['saved']:
            if s['init'] is not None and 'save' not in s['attrs'] and s['name'] in written:
                cand.setdefault(u.path, {}).setdefault(s['line'], set()).add(s['name'])

    found = []
    with tempfile.TemporaryDirectory() as tmp:
        for path, decls in sorted(cand.items()):
            lines = open(path).read().split('\n')
            names = set()
            for ln, vs in decls.items():
                new, _ = omp_generate.rewrite_decl(lines[ln - 1], vs)
                if new:
                    lines[ln - 1] = new
                    names |= vs
            f = os.path.join(tmp, os.path.basename(path))
            open(f, 'w').write('\n'.join(lines))
            p = subprocess.run([a.fc, '-cpp', '-ffree-line-length-none', '-O2', '-Wmaybe-uninitialized',
                                '-Wuninitialized', '-I' + a.mods, '-J' + tmp, '-c', f, '-o', os.devnull],
                               capture_output=True, text=True)
            loc = None
            for line in p.stderr.split('\n'):
                m = re.match(r'^.*?([\w.]+\.f90):(\d+):\d+:', line)
                if m:
                    loc = (m.group(1), int(m.group(2)))
                    continue
                m = re.search(r"[‘'](\w+)[’'] (may be|is) used uninitialized", line)
                if m and loc and m.group(1).lower() in names:
                    # the routine that owns the line
                    owner = [u for u in units if u.path == path and u.kind != 'module' and u.line <= loc[1]]
                    owner = max(owner, key=lambda u: u.line).name if owner else '?'
                    key = owner + ':' + m.group(1).lower()
                    if key not in known:
                        found.append((loc[0], loc[1], key))
    for f, ln, key in sorted(set(found)):
        print(f'local   {f}:{ln} {key} may be read before it is written in a call (review, then add it to allowed.txt)')
    print(f'{len(set(found))} new candidates')
    return 1 if found else 0


if __name__ == '__main__':
    sys.exit(main())
