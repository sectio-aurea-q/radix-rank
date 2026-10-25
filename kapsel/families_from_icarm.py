"""Faserungen aus der ICARM-Datenbank (database.json) lesen und als Familien-Dateien für kapsel-scan schreiben.
Jede Familie wird geprüft: die angegebene Spezialisierung T muss die eingetragene Kurve erzeugen (gleiches j).
Aufruf: python3 families_from_icarm.py <database.json> <zielordner>"""
import sys, os, re, json
from fractions import Fraction as Fr
import cypari2
pari = cypari2.Pari()

def py(expr):
    e = expr.strip().rstrip('.,;').replace('−', '-').replace('^', '**')
    e = re.sub(r'(\d)\s*T', r'\1*T', e); e = re.sub(r'(\d)\s*\(', r'\1*(', e); e = re.sub(r'\)\s*\(', ')*(', e)
    e = re.sub(r'\)\s*T', ')*T', e)
    return e

def parse(c):
    c = c.replace('\r', '')
    if not re.search(r'y\^2\s*[+-]\s*\(\s*Lx\s*\+\s*B\s*\)\s*y\s*=\s*x\(x\s*\+\s*D\)\(x\s*\+\s*E\)', c): return None
    if not re.search(r'B\s*=\s*pE\s*\+\s*qD\s*\+\s*pq\(L\s*-\s*p\s*-\s*q\)', c): return None
    f = {}
    for k in 'LpqDE':
        m = re.search(r'(?m)^\s*' + k + r'\s*=\s*([-+0-9T^*/() \s]+?)\s*[,.]?\s*$', c)
        if not m: return None
        f[k] = py(m.group(1))
    return f

def coeffs(expr):
    P = pari(expr.replace('**', '^'))
    v = [int(x) for x in pari.Vec(P)] if str(pari.type(P)) == 't_POL' else [int(P)]
    return v

def ainvs(f, u, v):
    T = Fr(u, v); env = {'T': T}
    L, p, q, D, E = (eval(f[k], {'__builtins__': {}}, env) for k in 'LpqDE')
    B = p * E + q * D + p * q * (L - p - q)
    a = [L, D + E, B, D * E]
    s = Fr(v) ** 2
    while True:
        r = [a[i] * s ** (i + 1) for i in range(4)]
        if all(x.denominator == 1 for x in r): return [int(x) for x in r] + [0]
        s *= v

if __name__ == '__main__':
    db = json.load(open(sys.argv[1])); out = sys.argv[2]; os.makedirs(out, exist_ok=True)
    seen = {}
    for x in db['curves']:
        c = x['commentary'] or ''
        f = parse(c); m = re.search(r'T\s*=\s*(-?\d+)\s*/\s*(\d+)', c)
        if not f or not m: continue
        u, v = int(m.group(1)), int(m.group(2))
        try:
            fam = {k: coeffs(f[k]) for k in 'LpqDE'}
            E1 = pari.ellinit(ainvs(f, u, v)); E2 = pari.ellinit([int(a) for a in x['ainvs']])
            if E1[12] != E2[12]: continue
        except Exception:
            continue
        key = json.dumps(fam, sort_keys=True)
        if key in seen: continue
        name = f"icarm{x['id']}"
        seen[key] = name
        fam.update(name=name, source=f"ICARM #{x['id']} (rank >= {x['rank_lower_bound']}, T = {u}/{v}, {x['submitter']})",
                   record_T=[u, v], record_rank=int(x['rank_lower_bound']))
        json.dump(fam, open(os.path.join(out, name + '.json'), 'w'), indent=1)
        print(name, fam['source'])
