#!/usr/bin/env python3
"""Write the source the OpenMP build compiles.

usage: omp_generate.py SRC_DIR OUT_DIR [--allowed allowed.txt] [--report]

Copies every .f90 file of SRC_DIR to OUT_DIR. In the routines that hru batches
run in parallel (the call trees of hru_control and obj_inflow), a local declared
with an initialiser (`real :: x = 0.`) has the implicit SAVE attribute, so all
threads would share one copy. The copy written to OUT_DIR declares it without
the initialiser and assigns the same value at the start of the routine, so each
call (and each thread) has its own copy.

Locals that must keep their value between calls are listed as `saved` in
allowed.txt and are left as they are. Locals that are never changed are
constants and are left as they are too.

A file is only rewritten when its content changes, so builds stay incremental.
"""
import argparse
import glob
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fscan
import ompaudit

ROOTS = ['hru_control', 'obj_inflow']


def split_code_comment(line):
    q = None
    for i, ch in enumerate(line):
        if q:
            if ch == q:
                q = None
        elif ch in ('"', "'"):
            q = ch
        elif ch == '!':
            return line[:i].rstrip(), line[i:]
    return line.rstrip(), ''


def rewrite_decl(line_text, names):
    """drop the initialisers of `names` from a one-line declaration;
    returns (new line, {name: (op, value, name as written)}) or (None, reason)"""
    code, comment = split_code_comment(line_text)
    if '::' not in code:
        return None, 'no ::'
    head, body = code.split('::', 1)
    out, inits = [], {}
    for e in fscan.split_top(body):
        m = re.match(r'^(\s*)([a-z_]\w*)(\s*\([^=]*?\))?(\s*\*\s*\d+)?\s*(=>|=)\s*(.*?)(\s*)$', e, re.I)
        if m and m.group(2).lower() in names:
            inits[m.group(2).lower()] = (m.group(5), m.group(6).strip(), m.group(2))
            out.append(m.group(1) + m.group(2) + (m.group(3) or '') + (m.group(4) or ''))
        else:
            out.append(e)
    if set(inits) != set(names):
        return None, 'entity not found'
    new = (head + '::' + ','.join(out)).rstrip()
    if comment:
        new = new + ' ' * max(1, line_text.index(comment) - len(new)) + comment
    return new, inits


def written_names(u):
    """names a routine assigns, reads into, allocates or passes as an actual argument"""
    out = {b for _, b, _, _, _ in u.assigns}
    for _, _, acts in u.callargs + u.funcargs:
        for act in acts:
            m = re.match(r'\s*([a-z_]\w*)', act, re.I)
            if m:
                out.add(m.group(1).lower())
    return out


def convert(src, allowed, res, byname):
    """returns {path: new text} for the files that change, and a list of skipped locals"""
    edits = {}      # path -> {decl line: set(names)}
    inserts = {}    # path -> [(executable line, unit name, [names])]
    skipped = []
    for n, r in res.items():
        u = byname[n]
        written = written_names(u)
        names = [s['name'] for s in r['saved']
                 if s['init'] is not None and 'save' not in s['attrs'] and s['name'] in written
                 and n + ':' + s['name'] not in allowed]
        if not names:
            continue
        if u.executable_start is None:
            skipped += [(u.path, n, v, 'no executable statements') for v in names]
            continue
        for v in names:
            edits.setdefault(u.path, {}).setdefault(u.decls[v]['line'], set()).add(v)
        inserts.setdefault(u.path, []).append((u.executable_start, n, names))

    out = {}
    for path in edits:
        lines = open(path).read().split('\n')
        initval = {}
        for ln, names in edits[path].items():
            txt = lines[ln - 1]
            if split_code_comment(txt)[0].endswith('&') or (ln < len(lines) and lines[ln].lstrip().startswith('&')):
                skipped += [(path, '?', v, 'declaration continues on the next line') for v in names]
                continue
            new, inits = rewrite_decl(txt, names)
            if new is None:
                skipped += [(path, '?', v, inits) for v in names]
                continue
            lines[ln - 1] = new
            for v, x in inits.items():
                initval[(ln, v)] = x
        # insert the assignments, last routine first so line numbers stay valid
        for exec_line, unit, names in sorted(inserts[path], key=lambda x: -x[0]):
            u = byname[unit]
            ind = re.match(r'^(\s*)', lines[exec_line - 1]).group(1)
            block = []
            for v in names:
                x = initval.get((u.decls[v]['line'], v))
                if x:
                    block.append(f'{ind}{x[2]} {x[0]} {x[1]}')
            lines[exec_line - 1:exec_line - 1] = block
        out[path] = '\n'.join(lines)
    return out, skipped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('src')
    ap.add_argument('out')
    ap.add_argument('--allowed', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'allowed.txt'))
    ap.add_argument('--report', action='store_true', help='list the converted and skipped locals')
    a = ap.parse_args()

    allowed = ompaudit.read_allowed(a.allowed)['saved']
    res, tree = ompaudit.analyse(a.src, ROOTS)
    units, byname, modules = ompaudit.build(a.src)
    changed, skipped = convert(a.src, allowed, res, byname)
    os.makedirs(a.out, exist_ok=True)
    nwritten = 0
    for f in sorted(glob.glob(os.path.join(a.src, '*.f90'))):
        text = changed.get(f)
        if text is None:
            text = open(f).read()
        dst = os.path.join(a.out, os.path.basename(f))
        if os.path.exists(dst) and open(dst).read() == text:
            continue
        with open(dst, 'w') as fh:
            fh.write(text)
        nwritten += 1
    # drop files that are no longer in the source
    keep = {os.path.basename(f) for f in glob.glob(os.path.join(a.src, '*.f90'))}
    for f in glob.glob(os.path.join(a.out, '*.f90')):
        if os.path.basename(f) not in keep:
            os.remove(f)
    print(f'omp_generate: {len(changed)} files converted, {nwritten} written to {a.out}')
    for s in skipped:
        print('omp_generate: left as is (review): %s %s %s (%s)' % s, file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
