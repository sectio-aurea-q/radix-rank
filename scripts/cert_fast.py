"""
Fast rigorous certification: same bound as scripts/cert_rap.py, with the prime sum enclosed by the Rust program
cert_sum (f64 interval arithmetic, all cores) instead of a Python loop over Arb balls.

Exact parts stay in Arb / exact rationals, as in certify.py:
  * the test function f and g = f*f (exact rational polynomials),
  * the Chebyshev coefficients of g on [0, 2L] (exact rationals), passed to Rust as outward-rounded f64 intervals,
  * 1/L and 2L (exact rationals -> outward f64 intervals), log 2 (f64 interval, checked against Arb here),
  * conductor term, archimedean term, normalisation.
The returned interval [lo, hi] for the prime sum is turned into an Arb ball containing it, so the final enclosure is
as rigorous as cert_rap.py's (both are reported side by side by --compare).

Usage: python3 cert_fast.py <curve.json> <rap_prefix> <stream_out.json> <key> [--compare]
"""
import json, math, os, subprocess, sys, tempfile, time
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import certify as C
import flint
from flint import arb, fmpq, fmpq_poly, fmpz

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CERT_SUM = os.path.join(ROOT, 'target', 'release', 'cert_sum')


def fq_to_frac(q):
    return Fraction(int(q.p), int(q.q))


def interval(q):
    """outward-rounded f64 interval [lo, hi] containing the exact rational q"""
    fr = fq_to_frac(q) if isinstance(q, fmpq) else Fraction(q)
    f = float(fr)
    lo = f if Fraction(f) <= fr else math.nextafter(f, -math.inf)
    hi = f if Fraction(f) >= fr else math.nextafter(f, math.inf)
    assert Fraction(lo) <= fr <= Fraction(hi)
    return lo, hi


def hx(x):
    return format(int.from_bytes(__import__('struct').pack('>d', x), 'big'), '016x')


def unhx(s):
    return __import__('struct').unpack('>d', int(s, 16).to_bytes(8, 'big'))[0]


def chebyshev(poly_x):
    """exact Chebyshev coefficients of a rational polynomial in x on [-1, 1]"""
    cs = [fmpq(0)] * (poly_x.degree() + 1)
    for n, c in enumerate(poly_x.coeffs()):
        c = fmpq(c)
        if c == 0:
            continue
        # x^n = 2^(1-n) sum_{k <= n/2} binom(n,k) T_{n-2k}, the T_0 term halved
        for k in range(n // 2 + 1):
            w = fmpq(math.comb(n, k), 2 ** (n - 1)) if n > 0 else fmpq(1)
            if n > 0 and 2 * k == n:
                w = w / 2
            cs[n - 2 * k] += c * w
    return cs


def ln2_interval():
    f = math.log(2)
    lo, hi = math.nextafter(f, -math.inf), math.nextafter(f, math.inf)
    l2 = arb(2).log()
    assert l2 > arb(fmpq(*Fraction(lo).as_integer_ratio())) and l2 < arb(fmpq(*Fraction(hi).as_integer_ratio()))
    return lo, hi


def main():
    cv = json.load(open(sys.argv[1]))
    prefix, out, key = sys.argv[2], json.load(open(sys.argv[3])), sys.argv[4]
    compare = '--compare' in sys.argv
    h = out[key]
    L = repr(h['l'])
    Lq = C.fq(L)
    g = C.g_poly(h['coeffs'], L)                               # exact, in u on [0, 2L]
    gx = g(fmpq_poly([Lq, Lq]))                                # g(L (x + 1)), x in [-1, 1]
    cheb = chebyshev(gx)
    # check: the Chebyshev series reproduces g exactly at a few rational points
    for xr in (fmpq(0), fmpq(1, 3), fmpq(-5, 7), fmpq(1)):
        t = [fmpq(1), xr]
        while len(t) < len(cheb):
            t.append(2 * xr * t[-1] - t[-2])
        assert sum(c * tj for c, tj in zip(cheb, t)) == gx(xr)
    job = dict(cheb=[[hx(a), hx(b)] for a, b in (interval(c) for c in cheb)],
               inv_l=[hx(x) for x in interval(1 / Lq)], two_l=[hx(x) for x in interval(2 * Lq)],
               ln2=[hx(x) for x in ln2_interval()],
               bad=sorted(int(b) for b in cv['bad'] if int(b) < 2 ** 63))
    t = time.time()
    with tempfile.NamedTemporaryFile('w', suffix='.json', delete=False) as tf:
        json.dump(job, tf)
    r = subprocess.run([CERT_SUM, tf.name, prefix], check=True, capture_output=True, text=True)
    os.unlink(tf.name)
    s = json.loads(r.stdout)
    lo, hi = unhx(s['lo']), unhx(s['hi'])
    S = arb(fmpq(*Fraction(lo).as_integer_ratio())).union(arb(fmpq(*Fraction(hi).as_integer_ratio())))
    secs = time.time() - t
    N = int(cv['N'])
    W = arb(g(fmpq(0))) * (arb(N).log() - 2 * (2 * arb.pi()).log()) + C.arch_term(g, L) - 2 * S
    f = sum((C.fq(x) * fmpq_poly.legendre_p(2 * k) for k, x in enumerate(h['coeffs'])), fmpq_poly(0))
    F = f.integral()
    b = W / arb(Lq * (F(fmpq(1)) - F(fmpq(-1)))) ** 2
    target = out['target']
    res = dict(label=cv['label'], key=key, L=h['l'], primes=s['primes'], terms=s['terms'],
               arb=b.str(10, radius=True), target=target, proven=bool(b.upper() < target),
               prime_sum_width=hi - lo, secs=round(secs, 1), method='cert_sum (f64 intervals)')
    print(json.dumps(res), flush=True)


if __name__ == '__main__':
    main()
