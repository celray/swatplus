#!/usr/bin/env python3
"""Lightweight scanner for the SWAT+ Fortran sources.

Builds, per program unit: uses, local declarations (with initialisers),
assignments (LHS base name), calls, function references and I/O statements.
Module units also record their variable declarations.

The parse is line based (comments stripped, continuations joined) and is
tuned to the SWAT+ coding style; it is not a general Fortran parser.
"""
import re, os, glob, json, sys
from collections import defaultdict

TYPE_RE = re.compile(r'^(integer|real|logical|character|complex|double\s+precision|type\s*\(|class\s*\()', re.I)
UNIT_START = re.compile(r'^(?:(?:recursive|pure|elemental|impure|module)\s+)*'
                        r'(?:(?:integer|real|logical|character|double\s+precision)(?:\s*\([^)]*\))?\s+)?'
                        r'(subroutine|function)\s+(\w+)', re.I)
MODULE_START = re.compile(r'^module\s+(?!procedure\b)(\w+)\s*$', re.I)
END_RE = re.compile(r'^end\s*(subroutine|function|module|program)?\b', re.I)
PROGRAM_START = re.compile(r'^program\s+(\w+)', re.I)
CONTAINS_RE = re.compile(r'^contains\s*$', re.I)
TYPEDEF_START = re.compile(r'^type\s*(?:,\s*[\w\s,()]*)?(?:::)?\s*(\w+)\s*$', re.I)
END_TYPE = re.compile(r'^end\s*type\b', re.I)
INTERFACE_START = re.compile(r'^(abstract\s+)?interface\b', re.I)
END_INTERFACE = re.compile(r'^end\s*interface\b', re.I)
USE_RE = re.compile(r'^use\s*(?:,\s*intrinsic\s*::)?\s*(\w+)\s*(?:,\s*only\s*:\s*(.*))?$', re.I)
CALL_RE = re.compile(r'\bcall\s+(\w+)', re.I)
IO_RE = re.compile(r'^(?:\d+\s+)?(write|print|open|close|read|inquire|rewind|backspace|flush)\b\s*(\(|\*|,)?', re.I)
KEYWORDS = set('''if then else elseif endif end do enddo while select case default call return exit cycle
stop go goto continue where elsewhere forall associate block allocate deallocate nullify
format write print read open close inquire rewind backspace flush and or not eq ne lt le gt ge true false
eqv neqv'''.split())
INTRINSICS = set('''abs acos aimag aint alog alog10 amax1 amin1 anint asin atan atan2 ceiling char cmplx cos cosh
dble dim dmax1 dmin1 dot_product dprod exp float floor iabs ichar index int len len_trim log log10 max max0
max1 maxloc maxval min min0 min1 minloc minval mod modulo nint real sign sin sinh size sqrt sum tan tanh trim
adjustl adjustr allocated associated present huge tiny epsilon count any all lbound ubound shape reshape
transpose matmul merge pack unpack spread achar iachar ishft iand ior ieor btest ibset ibclr not
selected_real_kind selected_int_kind kind range precision digits product cshift eoshift isnan
date_and_time cpu_time system_clock random_number random_seed get_command_argument command_argument_count
getenv get_environment_variable execute_command_line nearest exponent fraction scale spacing rrspacing
null c_loc c_f_pointer repeat scan verify lge lgt lle llt idnint idint ifix sngl new_line move_alloc
storage_size is_iostat_end is_iostat_eor norm2 erf erfc gamma log_gamma bessel_j0 hypot amod dsqrt dexp dlog
dabs dsin dcos dtan datan dlog10 floor ichar'''.split())


def strip_comment(line):
    out = []
    q = None
    for ch in line:
        if q:
            out.append(ch)
            if ch == q:
                q = None
        else:
            if ch in ('"', "'"):
                q = ch
                out.append(ch)
            elif ch == '!':
                break
            else:
                out.append(ch)
    return ''.join(out)


def logical_lines(path):
    """Yield (lineno, text, raw) with comments stripped and continuations joined.
    Keeps OpenMP sentinels (!$omp) as separate 'omp' lines."""
    with open(path, errors='replace') as fh:
        raw = fh.readlines()
    buf = ''
    start = None
    for i, r in enumerate(raw, 1):
        s = r.rstrip('\n')
        st = s.strip()
        if st.lower().startswith('!$omp'):
            yield i, st, 'omp'
            continue
        if st.startswith('#'):
            yield i, st, 'cpp'
            continue
        t = strip_comment(s).rstrip()
        if not t.strip():
            if buf:
                continue
            continue
        tt = t.strip()
        if buf:
            if tt.startswith('&'):
                tt = tt[1:]
            buf += ' ' + tt
        else:
            buf = tt
            start = i
        if buf.endswith('&'):
            buf = buf[:-1]
            continue
        # split on ; outside strings
        for part in split_semicolons(buf):
            if part.strip():
                yield start, part.strip(), 'code'
        buf = ''


def split_semicolons(s):
    parts, cur, q = [], [], None
    for ch in s:
        if q:
            cur.append(ch)
            if ch == q:
                q = None
        elif ch in ('"', "'"):
            q = ch
            cur.append(ch)
        elif ch == ';':
            parts.append(''.join(cur))
            cur = []
        else:
            cur.append(ch)
    parts.append(''.join(cur))
    return parts


def blank_strings(s):
    out, q = [], None
    for ch in s:
        if q:
            if ch == q:
                q = None
                out.append(ch)
            else:
                out.append('_')
        elif ch in ('"', "'"):
            q = ch
            out.append(ch)
        else:
            out.append(ch)
    return ''.join(out)


def split_top(s, sep=','):
    parts, depth, cur = [], 0, []
    for ch in s:
        if ch in '([':
            depth += 1
        elif ch in ')]':
            depth -= 1
        if ch == sep and depth == 0:
            parts.append(''.join(cur))
            cur = []
        else:
            cur.append(ch)
    parts.append(''.join(cur))
    return parts


def match_paren(s, i):
    """s[i] == '(' ; return index of matching ')'"""
    depth = 0
    for j in range(i, len(s)):
        if s[j] == '(':
            depth += 1
        elif s[j] == ')':
            depth -= 1
            if depth == 0:
                return j
    return len(s) - 1


ASSIGN_LHS = re.compile(r'^([a-z_]\w*)\s*(\(|%|=(?!=))', re.I)


def find_assignment(stmt):
    """Return (lhs_text, rhs_text) if stmt is an assignment, else None."""
    s = stmt
    m = ASSIGN_LHS.match(s)
    if not m:
        return None
    name = m.group(1).lower()
    if name in ('if', 'do', 'select', 'case', 'where', 'forall', 'call', 'print', 'write', 'read',
                'allocate', 'deallocate', 'return', 'else', 'elseif', 'go', 'goto', 'open', 'close',
                'format', 'inquire', 'associate', 'block', 'exit', 'cycle', 'stop', 'nullify', 'end',
                'while', 'rewind', 'backspace', 'flush', 'use', 'implicit', 'contains', 'data'):
        return None
    # walk the designator: name ( ... ) % name ( ... ) ...
    i = 0
    n = len(s)
    while i < n:
        mm = re.match(r'\s*[a-z_]\w*\s*', s[i:], re.I)
        if not mm:
            return None
        i += mm.end()
        while i < n and s[i] == '(':
            i = match_paren(s, i) + 1
            while i < n and s[i] == ' ':
                i += 1
        if i < n and s[i] == '%':
            i += 1
            continue
        break
    rest = s[i:].lstrip()
    if rest.startswith('=') and not rest.startswith('=='):
        return s[:i].strip(), rest[1:].strip()
    if rest.startswith('=>'):
        return s[:i].strip(), rest[2:].strip()
    return None


IF_PREFIX = re.compile(r'^(?:\d+\s+)?(?:else\s*)?if\s*\(', re.I)


def strip_if_prefix(stmt):
    """'if (cond) stmt' -> (cond, stmt) ; returns (None, stmt) otherwise"""
    m = IF_PREFIX.match(stmt)
    if not m:
        return None, stmt
    j = match_paren(stmt, m.end() - 1)
    rest = stmt[j + 1:].strip()
    if rest.lower().startswith('then') or rest == '':
        return stmt[m.end():j], ''
    return stmt[m.end():j], rest


class Unit:
    def __init__(self, kind, name, path, line, parent=None):
        self.kind, self.name, self.path, self.line = kind, name.lower(), path, line
        self.parent = parent
        self.uses = {}          # module -> None (all) or set(only names)
        self.decls = {}         # var -> dict(line, init, attrs, text)
        self.assigns = []       # (line, base, lhs, rhs, stmt)
        self.calls = []         # (line, name)
        self.callargs = []      # (line, name, [actual args])
        self.funcargs = []      # (line, name, [actual args]) every name( ... ) reference
        self.refs = defaultdict(list)   # identifier -> lines (rhs/expression references)
        self.io = []            # (line, stmt)
        self.omp = []           # (line, text)
        self.args = []
        self.result = None
        self.contains = []
        self.typedefs = set()
        self.executable_start = None

    def to_json(self):
        return dict(kind=self.kind, name=self.name, path=self.path, line=self.line,
                    parent=self.parent.name if self.parent else None,
                    uses={k: (sorted(v) if v is not None else None) for k, v in self.uses.items()},
                    decls=self.decls, assigns=self.assigns, calls=self.calls,
                    io=self.io, omp=self.omp, args=self.args, result=self.result,
                    executable_start=self.executable_start)


DECL_SPLIT = re.compile(r'::')


def parse_decl(stmt):
    """Return (typespec, attrs, [(name, dims, init)]) or None"""
    if not TYPE_RE.match(stmt):
        return None
    low = stmt.lower()
    if low.startswith('type') and not re.match(r'^type\s*\(', low):
        return None
    if re.match(r'^(integer|real|logical|character|double\s+precision|complex)(\s*\*\s*\d+)?(\s*\([^)]*\))?\s*function\b', low):
        return None
    if '::' in stmt:
        head, body = stmt.split('::', 1)
    else:
        # old style: real x, y
        m = re.match(r'^((?:integer|real|logical|character|double\s+precision)(?:\s*\*\s*\d+)?(?:\s*\([^)]*\))?)\s+(.*)$', stmt, re.I)
        if not m:
            return None
        head, body = m.group(1), m.group(2)
    parts = split_top(head)
    typespec = parts[0].strip()
    attrs = [p.strip().lower() for p in parts[1:]]
    out = []
    for item in split_top(body):
        item = item.strip()
        if not item:
            continue
        init = None
        # name(dims) = init  or name = init or name => null()
        m = re.match(r'^([a-z_]\w*)\s*(\([^=]*?\))?\s*(\*\s*\d+)?\s*(=>|=)?\s*(.*)$', item, re.I)
        if not m:
            continue
        nm = m.group(1).lower()
        dims = m.group(2)
        if m.group(4):
            init = m.group(5)
        out.append((nm, dims, init))
    return typespec, attrs, out


def scan_file(path):
    units = []
    stack = []
    in_type = 0
    in_iface = 0
    for line, stmt, kind in logical_lines(path):
        cur = stack[-1] if stack else None
        if kind == 'omp':
            if cur:
                cur.omp.append((line, stmt))
            continue
        if kind == 'cpp':
            continue
        s = blank_strings(stmt)
        low = s.lower()
        if in_iface:
            if END_INTERFACE.match(low):
                in_iface -= 1
            continue
        if INTERFACE_START.match(low) and not low.startswith('interface operator') and True:
            in_iface += 1
            continue
        if in_type:
            if END_TYPE.match(low):
                in_type -= 1
            continue
        m = MODULE_START.match(low)
        if m and not stack:
            u = Unit('module', m.group(1), path, line)
            units.append(u)
            stack.append(u)
            continue
        m = PROGRAM_START.match(low)
        if m and not stack:
            u = Unit('program', m.group(1), path, line)
            units.append(u)
            stack.append(u)
            continue
        m = UNIT_START.match(low)
        if m and not low.startswith('end'):
            u = Unit(m.group(1).lower(), m.group(2), path, line, parent=cur)
            am = re.search(r'\b' + m.group(2) + r'\s*\(([^)]*)\)', s, re.I)
            if am:
                u.args = [a.strip().lower() for a in am.group(1).split(',') if a.strip()]
            rm = re.search(r'result\s*\(\s*(\w+)\s*\)', low)
            u.result = rm.group(1) if rm else (u.name if u.kind == 'function' else None)
            if cur:
                cur.contains.append(u)
            units.append(u)
            stack.append(u)
            continue
        if re.match(r'^end\s*((subroutine|function|module|program)(\s+\w+)?)?\s*$', low):
            if stack:
                stack.pop()
            continue
        if cur is None:
            continue
        if TYPEDEF_START.match(low) and not re.match(r'^type\s*\(', low):
            in_type += 1
            tm = TYPEDEF_START.match(low)
            cur.typedefs.add(tm.group(1))
            continue
        if CONTAINS_RE.match(low):
            continue
        m = USE_RE.match(low)
        if m:
            mod = m.group(1)
            if m.group(2) is not None:
                names = set()
                for x in m.group(2).split(','):
                    x = x.strip()
                    if '=>' in x:
                        x = x.split('=>')[0].strip()   # local name
                    if x:
                        names.add(x)
                prev = cur.uses.get(mod, set())
                cur.uses[mod] = None if prev is None else (prev | names)
            else:
                cur.uses[mod] = None
            continue
        if low.startswith('implicit'):
            continue
        d = parse_decl(s)
        if d and cur.executable_start is None:
            typespec, attrs, items = d
            for nm, dims, init in items:
                cur.decls[nm] = dict(line=line, init=init, attrs=attrs, type=typespec.lower(),
                                     dims=dims, text=stmt)
            continue
        if re.match(r'^(external|intrinsic|save|data|parameter|dimension|equivalence|common|namelist|private|public|allocatable|target|pointer|optional|intent|sequence|procedure)\b', low):
            if low.startswith('save') and cur.executable_start is None:
                cur.decls.setdefault('__save__', dict(line=line, init=None, attrs=['save'], type='', dims=None, text=stmt))
            continue
        if re.match(r'^\d+\s+format\b', low) or low.startswith('format'):
            continue
        if cur.kind != 'module' and cur.executable_start is None:
            cur.executable_start = line
        # executable statement
        for cm in CALL_RE.finditer(low):
            cur.calls.append((line, cm.group(1)))
            k = cm.end()
            while k < len(s) and s[k] == ' ':
                k += 1
            if k < len(s) and s[k] == '(':
                e = match_paren(s, k)
                cur.callargs.append((line, cm.group(1), [a.strip() for a in split_top(s[k + 1:e])]))
            else:
                cur.callargs.append((line, cm.group(1), []))
        for fm in re.finditer(r'\b([a-z_]\w*)\s*\(', low):
            nm = fm.group(1)
            if nm in KEYWORDS or nm in INTRINSICS:
                continue
            k = fm.end() - 1
            e = match_paren(s, k)
            cur.funcargs.append((line, nm, [a.strip() for a in split_top(s[k + 1:e])]))
        cond, body = strip_if_prefix(s)
        bl = body.lower()
        iom = IO_RE.match(bl)
        if iom and (iom.group(2) or bl.split()[0] in ('flush', 'close', 'rewind')):
            cur.io.append((line, stmt))
        a = find_assignment(body)
        if a:
            lhs, rhs = a
            base = re.match(r'\s*([a-z_]\w*)', lhs, re.I).group(1).lower()
            cur.assigns.append((line, base, lhs, rhs, stmt))
        # do-loop index assignment
        dm = re.match(r'^(?:\w+\s*:\s*)?do\s+([a-z_]\w*)\s*=', bl)
        if dm:
            cur.assigns.append((line, dm.group(1).lower(), dm.group(1), '', stmt))
        # read statement targets / allocate targets
        rm = re.match(r'^read\s*\((.*?)\)\s*(.*)$', bl)
        if rm:
            for it in split_top(rm.group(2)):
                im = re.match(r'\s*\(*\s*([a-z_]\w*)', it)
                if im:
                    cur.assigns.append((line, im.group(1), it.strip(), '<read>', stmt))
        am = re.match(r'^(allocate|deallocate)\s*\((.*)\)\s*$', bl)
        if am:
            for it in split_top(am.group(2)):
                im = re.match(r'\s*([a-z_]\w*)', it)
                if im and '=' not in it:
                    cur.assigns.append((line, im.group(1), it.strip(), '<' + am.group(1) + '>', stmt))
        for idm in re.finditer(r'\b([a-z_]\w*)\s*(\()?', low):
            cur.refs[idm.group(1)].append(line)
    return units


def scan_dir(src):
    allunits = []
    for f in sorted(glob.glob(os.path.join(src, '*.f90'))):
        allunits.extend(scan_file(f))
    return allunits


if __name__ == '__main__':
    units = scan_dir(sys.argv[1])
    print(len(units), 'units')
    kinds = defaultdict(int)
    for u in units:
        kinds[u.kind] += 1
    print(dict(kinds))
