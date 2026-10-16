"""
Rigorous (Arb) certification of a rank bound from compact .rap files, in parallel over all cores.

Same formula and arithmetic as scripts/certify.py (256-bit balls):
    r_an(E) <= W_E(f) / f^(0)^2,   W_E(f) = g(0)(log N - 2 log 2pi) + Arch(g) - 2 sum_n Lambda(n) c_n n^{-1/2} g(log n)
Each worker sums the prime part of whole .rap blocks as an Arb ball; the balls are passed back EXACTLY
(midpoint and radius as integer mantissa/exponent) and added, so the final enclosure is as rigorous as a
single sequential sum.

Usage: python3 cert_rap.py <curve.json> <rap_prefix> <stream_out.json> <key> [workers]
  curve.json : {"label", "N": "<exact conductor>", "bad": [...], "rank"}
  key        : radix_hit | std_hit | radix_end | std_end   (fields of the radix-stream output)
"""
import glob, json, math, os, sys, time
from multiprocessing import Pool
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import certify as C          # g_poly, arch_term, horner_arb, fq  (sets flint.ctx.prec = 256)
import flint
from flint import arb, fmpq, fmpq_poly

G = {}


def read_rap(path):
    with open(path, 'rb') as f:
        b = f.read()
    assert b[:4] == b'RAP1', path
    lo, hi, n = np.frombuffer(b, dtype='<u8', count=3, offset=4)
    gaps = np.frombuffer(b, dtype='<u2', count=int(n), offset=28).astype(np.uint64)
    ps = int(lo) + np.cumsum(gaps)
    aps = np.frombuffer(b, dtype='<i4', count=int(n), offset=28 + 2 * int(n))
    return ps, aps


def init(L, coeffs, bad):
    flint.ctx.prec = 256
    g = C.g_poly(coeffs, L)
    G['gc'] = [fmpq(x) for x in g.coeffs()]
    G['twoL'] = 2 * float(C.fq(L))
    G['bad'] = set(bad)


def to_wire(x):
    return x.mid().man_exp(), x.rad().man_exp()


def from_wire(w):
    (mm, me), (rm, re) = w
    return arb(arb(mm) * arb(2) ** me, arb(rm) * arb(2) ** re)


def block_sum(path):
    gc, twoL, bad = G['gc'], G['twoL'], G['bad']
    ps, aps = read_rap(path)
    n = int(np.searchsorted(ps, math.exp(twoL) * 1.000001)) + 1
    S, cnt = arb(0), 0
    for p, ap in zip(ps[:n].tolist(), aps[:n].tolist()):
        lp = arb(p).log()
        if float(lp.lower()) >= twoL:
            break
        cnt += 1
        s1 = arb(ap) / arb(p).sqrt()
        s_prev, s_cur, k = arb(2), s1, 1
        while True:
            ln = k * lp
            if float(ln.lower()) >= twoL:
                break
            if float(ln.upper()) >= twoL:
                raise ValueError(f"prime power {p}^{k} too close to the window edge")
            ck = s1 ** k if p in bad else s_cur
            S += lp * ck / arb(p) ** (arb(k) / 2) * C.horner_arb(gc, ln)
            if p not in bad:
                s_prev, s_cur = s_cur, s1 * s_cur - s_prev
            k += 1
    return to_wire(S), cnt


def main():
    cv = json.load(open(sys.argv[1]))
    prefix, out, key = sys.argv[2], json.load(open(sys.argv[3])), sys.argv[4]
    workers = int(sys.argv[5]) if len(sys.argv) > 5 else os.cpu_count()
    h = out[key]
    assert h is not None, f"{key}: no hit in {sys.argv[3]}"
    L = repr(h['l'])
    twoL = 2 * float(C.fq(L))
    blocks = []
    for f in glob.glob(prefix + '.*.rap'):
        lo, hi = map(int, f[len(prefix) + 1:-4].split('-'))
        if lo < math.exp(twoL) * 1.000001:
            assert os.path.exists(f[:-4] + '.ok'), f
            blocks.append((lo, hi, f))
    blocks.sort()
    end = 0
    for lo, hi, _ in blocks:
        assert lo == end, f"gap in the a_p data at {end}"
        end = hi
    assert end > math.exp(twoL), f"a_p data end at {end:.3e} < e^(2L) = {math.exp(twoL):.3e}"
    t = time.time()
    with Pool(workers, initializer=init, initargs=(L, h['coeffs'], [int(b) for b in cv['bad']])) as pool:
        res = pool.map(block_sum, [f for _, _, f in blocks], chunksize=1)
    S = arb(0)
    cnt = 0
    for w, c in res:
        S += from_wire(w)
        cnt += c
    init(L, h['coeffs'], [])
    g = C.g_poly(h['coeffs'], L)
    N = int(cv['N'])
    W = arb(g(fmpq(0))) * (arb(N).log() - 2 * (2 * arb.pi()).log()) + C.arch_term(g, L) - 2 * S
    f = sum((C.fq(x) * fmpq_poly.legendre_p(2 * k) for k, x in enumerate(h['coeffs'])), fmpq_poly(0))
    F = f.integral()
    b = W / arb(C.fq(L) * (F(fmpq(1)) - F(fmpq(-1)))) ** 2
    target = out['target']
    res = dict(label=cv['label'], key=key, L=h['l'], primes=cnt, arb=b.str(10, radius=True), target=target,
               proven=bool(b.upper() < target), workers=workers, secs=round(time.time() - t))
    print(json.dumps(res), flush=True)


if __name__ == '__main__':
    main()
