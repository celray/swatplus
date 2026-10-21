#!/usr/bin/env python3
"""Audit the call tree under given root routines for OpenMP hazards.

usage: ompaudit.py SRC_DIR ROOT [ROOT ...] [--json out.json] [--report out.md]

Reports, for every routine reachable from the roots:
  saved locals   locals with an initialiser or SAVE (shared between threads)
  module writes  module variables assigned (base name), with the module
  io             write/print/open/... statements
"""
import sys, re, json, argparse
from collections import defaultdict
sys.path.insert(0, __file__.rsplit('/', 1)[0])
import fscan


def build(src):
    units = fscan.scan_dir(src)
    byname = {}
    modules = {}
    for u in units:
        if u.kind == 'module':
            modules[u.name] = u
        else:
            byname.setdefault(u.name, u)
    return units, byname, modules


def module_exports(modules):
    """module -> {var: declaring module}; transitive through use."""
    cache = {}

    def exp(m, seen=()):
        if m in cache:
            return cache[m]
        if m not in modules or m in seen:
            return {}
        u = modules[m]
        out = {}
        for mod, only in u.uses.items():
            e = exp(mod, seen + (m,))
            for k, v in e.items():
                if only is None or k in only:
                    out[k] = v
        for k, d in u.decls.items():
            if k == '__save__':
                continue
            if 'parameter' in d['attrs']:
                continue
            out[k] = m
        cache[m] = out
        return out

    for m in modules:
        exp(m)
    return cache


def visible(u, exports):
    vis = {}
    for mod, only in u.uses.items():
        e = exports.get(mod, {})
        for k, v in e.items():
            if only is None or k in only:
                vis[k] = v
    return vis


_dw = {}


def dummy_written(cu):
    if cu.name in _dw:
        return _dw[cu.name]
    out = set()
    for d in cu.args:
        dd = cu.decls.get(d)
        if dd and any(a.replace(' ', '') in ('intent(out)', 'intent(inout)') for a in dd['attrs']):
            out.add(d)
    for line, base, lhs, rhs, stmt in cu.assigns:
        if base in cu.args:
            out.add(base)
    _dw[cu.name] = out
    return out


def callees(u, byname, funcnames):
    out = set()
    for _, n in u.calls:
        out.add(n)
    for ident in u.refs:
        if ident in funcnames and not (ident in u.decls and u.decls[ident]['dims']) and ident != u.name:
            out.add(ident)
    return {c for c in out if c in byname}


def reach(roots, byname, funcnames, stop=()):
    seen = {}
    todo = [(r, None) for r in roots]
    while todo:
        n, par = todo.pop()
        if n in seen or n not in byname:
            continue
        seen[n] = par
        if n in stop:
            continue
        for c in callees(byname[n], byname, funcnames):
            if c not in seen:
                todo.append((c, n))
    return seen


def analyse(src, roots, stop=()):
    units, byname, modules = build(src)
    exports = module_exports(modules)
    funcnames = {u.name for u in units if u.kind == 'function'}
    tree = reach(roots, byname, funcnames, stop)
    res = {}
    for n in sorted(tree):
        u = byname[n]
        locs = set(u.decls) | set(u.args)
        if u.result:
            locs.add(u.result)
        if u.parent is not None and u.parent.kind != 'module':
            locs |= set(u.parent.decls)  # host association (internal procedures)
        vis = visible(u, exports)
        if u.parent is not None and u.parent.kind == 'module':
            # module procedure: sees its module's variables and uses
            pv = visible(u.parent, exports)
            pv.update({k: u.parent.name for k, d in u.parent.decls.items()
                       if k != '__save__' and 'parameter' not in d['attrs']})
            for k, v in pv.items():
                vis.setdefault(k, v)
        saved = []
        has_save_all = '__save__' in u.decls and u.decls['__save__']['text'].strip().lower() == 'save'
        for k, d in u.decls.items():
            if k == '__save__':
                continue
            if 'parameter' in d['attrs']:
                continue
            if k in u.args:
                continue
            if d['init'] is not None or 'save' in d['attrs'] or has_save_all:
                saved.append(dict(name=k, line=d['line'], init=d['init'], attrs=d['attrs'],
                                  type=d['type'], dims=d['dims']))
        mw = defaultdict(list)
        for line, base, lhs, rhs, stmt in u.assigns:
            if base in locs:
                continue
            if base in vis:
                mw[base].append(dict(line=line, lhs=lhs, rhs=rhs[:160], module=vis[base]))
        # writes through actual arguments of calls and function references
        argrefs = [(l, n, a) for l, n, a in u.callargs] + \
                  [(l, n, a) for l, n, a in u.funcargs if n in funcnames and not (n in u.decls and u.decls[n]['dims'])]
        for line, callee, actuals in argrefs:
            cu = byname.get(callee)
            if cu is None:
                continue
            wd = dummy_written(cu)
            for pos, act in enumerate(actuals):
                if '=' in act and not act.startswith('='):
                    kw, act = act.split('=', 1)
                    kw = kw.strip().lower()
                    if kw not in cu.args:
                        continue
                    dummy = kw
                else:
                    if pos >= len(cu.args):
                        continue
                    dummy = cu.args[pos]
                if dummy not in wd:
                    continue
                m = re.match(r'\s*([a-z_]\w*)\s*(\(|%|$)', act, re.I)
                if not m:
                    continue
                base = m.group(1).lower()
                if base in locs or base not in vis:
                    continue
                mw[base].append(dict(line=line, lhs=act.strip(), rhs='<arg of %s:%s>' % (callee, dummy),
                                     module=vis[base]))
        res[n] = dict(path=u.path.rsplit('/', 1)[1], line=u.line, kind=u.kind,
                      parent=tree[n], saved=saved, modwrites=dict(mw),
                      io=[dict(line=l, stmt=s[:200]) for l, s in u.io],
                      omp=u.omp, exec_start=u.executable_start)
    return res, tree


HRU_INDEX = {'j', 'ihru', 'iihru', 'hru_id', 'icmd'}


def threadprivate_vars(src, built=None):
    """module variables named in !$omp threadprivate directives"""
    units, byname, modules = built or build(src)
    out = set()
    for u in list(modules.values()) + units:
        text = ' '.join(t for _, t in u.omp)
        for m in re.finditer(r'threadprivate\s*\(([^)]*)\)', text.replace('!$omp&', ' ').replace('&', ' '), re.I):
            out |= {x.strip().lower() for x in m.group(1).split(',') if x.strip()}
    return out


def read_allowed(path):
    """allowed.txt: '<kind> <name>  reason' lines; kinds: shared, saved, local, io"""
    allowed = {'shared': {}, 'saved': {}, 'local': {}, 'io': {}}
    for line in open(path):
        line = line.split('#', 1)[0].rstrip()
        if not line.strip():
            continue
        parts = line.split(None, 2)
        if len(parts) < 3:
            sys.exit('allowed.txt: kind, name and reason needed: ' + line)
        allowed[parts[0]][parts[1].lower()] = parts[2]
    return allowed


def per_object(lhs):
    """true if a subscript of the designator is the current hru or object"""
    for m in re.finditer(r'\(([^()]*)\)', lhs):
        for sub in m.group(1).split(','):
            if sub.strip().lower() in HRU_INDEX:
                return True
    return False


def check(src, roots, allowed_path):
    res, tree = analyse(src, roots)
    units, byname, modules = build(src)
    tp = threadprivate_vars(src)
    allowed = read_allowed(allowed_path)
    problems = []
    for n, r in sorted(res.items()):
        u = byname[n]
        written = {b for _, b, _, _, _ in u.assigns}
        for _, _, acts in u.callargs:
            for act in acts:
                m = re.match(r'\s*([a-z_]\w*)', act, re.I)
                if m:
                    written.add(m.group(1).lower())
        for sv in r['saved']:
            if sv['name'] not in written:
                continue        # never changed: a constant, safe to share
            key = n + ':' + sv['name']
            if key not in allowed['saved']:
                problems.append(f"saved   {r['path']}:{sv['line']} {key} keeps its value between calls "
                                f"(make it automatic with fix_saved_locals.py, or allow it)")
        for v, ws in r['modwrites'].items():
            if v in tp or v in allowed['shared']:
                continue
            for w in ws:
                if per_object(w['lhs']):
                    continue
                problems.append(f"shared  {r['path']}:{w['line']} {w['module']}.{v} written as {w['lhs'][:60]}")
                break
        if r['io'] and n not in allowed['io']:
            problems.append(f"io      {r['path']}:{r['io'][0]['line']} {n} writes files ({len(r['io'])} statements)")
    return problems


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('src')
    ap.add_argument('roots', nargs='+')
    ap.add_argument('--stop', nargs='*', default=[])
    ap.add_argument('--json')
    ap.add_argument('--check', metavar='ALLOWED', help='report hazards not in this allowed list; exit 1 if any')
    a = ap.parse_args()
    if a.check:
        problems = check(a.src, a.roots, a.check)
        for pr in problems:
            print(pr)
        print(f'{len(problems)} problems')
        sys.exit(1 if problems else 0)
    res, tree = analyse(a.src, a.roots, a.stop)
    if a.json:
        with open(a.json, 'w') as fh:
            json.dump(res, fh, indent=1)
    nsaved = sum(len(r['saved']) for r in res.values())
    mods = defaultdict(set)
    for n, r in res.items():
        for v, ws in r['modwrites'].items():
            mods[ws[0]['module']].add(v)
    nio = sum(len(r['io']) for r in res.values())
    print(f'{len(res)} routines reachable, {nsaved} saved locals, '
          f'{sum(len(v) for v in mods.values())} module variables written, {nio} io statements')
    for m in sorted(mods):
        print(f'  {m}: {len(mods[m])}')


if __name__ == '__main__':
    main()
