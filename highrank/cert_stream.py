"""Memory-lean Arb certification: same formula as radix-rank/scripts/certify.py, prime sum streamed from binary files."""
import json, math, sys, time, numpy as np
sys.path.insert(0, 'radix-rank/scripts')
import certify as C
from flint import arb, fmpq, fmpq_poly
cid, outfile, key = sys.argv[1], sys.argv[2], sys.argv[3]
d = json.loads(open(outfile).readline()); h = d[key]
c = json.load(open(f'icarm{cid}_3e7.json')); N = int(c['N_str']); bad = set(int(p) for p in json.load(open('database.json'))['curves'] and [])
bad = set(int(q) for q in c['bad'])
L = repr(h['l']); g = C.g_poly(h['coeffs'], L); gc = [fmpq(x) for x in g.coeffs()]
g0 = arb(g(fmpq(0))); W = g0 * (arb(N).log() - 2 * (2 * arb.pi()).log()) + C.arch_term(g, L)
twoL = 2 * float(C.fq(L))
p_all = np.fromfile(f'icarm{cid}.p.u64', dtype='<u8'); a_all = np.fromfile(f'icarm{cid}.ap.i64', dtype='<i8')
n = int(np.searchsorted(p_all, math.exp(twoL) * 1.000001)) + 1
S = arb(0); t = time.time()
for p, ap in zip(p_all[:n].tolist(), a_all[:n].tolist()):
    lp = arb(p).log()
    if float(lp.lower()) >= twoL: break
    s1 = arb(ap) / arb(p).sqrt(); s_prev, s_cur, k = arb(2), s1, 1
    while True:
        ln = k * lp
        if float(ln.lower()) >= twoL: break
        if float(ln.upper()) >= twoL: raise ValueError("prime power too close to window edge")
        ck = s1 ** k if p in bad else s_cur
        S += lp * ck / arb(p) ** (arb(k) / 2) * C.horner_arb(gc, ln)
        if p not in bad: s_prev, s_cur = s_cur, s1 * s_cur - s_prev
        k += 1
W -= 2 * S
f = sum((C.fq(x) * fmpq_poly.legendre_p(2 * k) for k, x in enumerate(h['coeffs'])), fmpq_poly(0)); F = f.integral()
b = W / arb(C.fq(L) * (F(fmpq(1)) - F(fmpq(-1)))) ** 2
target = d["rank"] + (2 if key.endswith("r2") else 1)
print(f"#{cid} {key} L={h['l']} primes={h['primes']} Arb {b.str(10, radius=True)} < {target}: "
      f"{'PROVEN' if b.upper() < target else 'NOT proven'} ({time.time()-t:.0f}s)", flush=True)
