"""
Beweisfeste Prüfung (Arb-Ballarithmetik via python-flint) einer Rangschranke:

   r_an(E) <= W_E(f) / f^(0)^2        (unter GRH, für jede reelle gerade f mit Träger in [-L, L])

f = sum_k c_k P_{2k}(x/L) mit exakt rationalen c_k (die Gleitkommazahlen aus dem Rust-Kern),
g = f * f als exaktes Polynom in Q[u] auf [0, 2L];
alle transzendenten Größen (log, sqrt, Euler-gamma, archimedisches Integral) als Arb-Bälle.
Das archimedische Integral: [0, d] über die Bernoulli-Reihe von u/(e^u-1) mit Restgliedschranke,
[d, 2L] mit acb.integral (rigoros).
"""
import json, math, sys
from fractions import Fraction
import flint
from flint import arb, acb, fmpq, fmpq_poly

flint.ctx.prec = 256

def fq(x):
    fr = Fraction(x) if not isinstance(x, str) else Fraction(x)
    return fmpq(fr.numerator, fr.denominator)

def g_poly(coeffs, L):
    """g(u) = int f(x) f(x-u) dx auf [0, 2L] als Polynom in u (exakt rational)."""
    f = sum((fq(c) * fmpq_poly.legendre_p(2 * k) for k, c in enumerate(coeffs)), fmpq_poly(0))
    s = fmpq_poly([0, 1])
    sigma_minus_1 = fmpq_poly([-1, 1])
    G = fmpq_poly(0)
    deriv, fact = f, fmpq(1)
    m = 0
    while not deriv.is_zero():
        F = (f * deriv).integral()                       # Stammfunktion in s
        term = (F(fmpq(1)) - F(sigma_minus_1))           # int_{sigma-1}^{1} f(s) f^(m)(s) ds
        G += fmpq_poly([0, -1]) ** m * term * (fmpq(1) / fact)
        m += 1; fact *= m; deriv = deriv.derivative()
    # g(u) = L * G(u / L)
    return fq(L) * G(fmpq_poly([0, fmpq(1) / fq(L)]))

def horner_acb(coefs, z):
    acc = acb(0)
    for c in reversed(coefs):
        acc = acc * z + acb(arb(c))
    return acc

def horner_arb(coefs, x):
    acc = arb(0)
    for c in reversed(coefs):
        acc = acc * x + arb(c)
    return acc

def arch_term(g, L):
    g0 = g(fmpq(0))
    H = (fmpq_poly([g0]) - g)                             # g(0) - g(u), Nullstelle bei u = 0
    H = fmpq_poly(H.coeffs()[1:]) if H.length() > 1 else fmpq_poly(0)   # (g0 - g(u)) / u
    Hc = [fmpq(c) for c in H.coeffs()]
    d, M = fmpq(1, 4), 80
    # [0, d]: u/(e^u-1) = sum_n B_n u^n / n!,  |B_n/n!| <= 4 (2pi)^-n
    ser = fmpq_poly([fmpq.bernoulli(n) * fmpq(1, math.factorial(n)) for n in range(M)])
    P = (H * ser).truncate(H.length() + M).integral()
    I0 = arb(P(d))
    absH = sum(abs(arb(c)) * arb(d) ** i for i, c in enumerate(Hc))   # max |H| auf [0,d]
    rad = arb(4) * (arb(d) / (2 * arb.pi())) ** M / (1 - arb(d) / (2 * arb.pi())) * absH * arb(d)
    I0 = I0 + arb(0, rad.upper())
    # [d, 2L]: rigoroses Integral des analytischen Integranden H(z) * z / (e^z - 1)
    f = lambda z, analytic: horner_acb(Hc, z) * z / z.expm1()
    I1 = acb.integral(f, acb(arb(d)), acb(2 * arb(fq(L)))).real
    gam = arb.const_euler()
    tail = (1 - (-2 * arb(fq(L))).exp()).log()
    return 2 * (-gam * arb(g0) + I0 + I1 - arb(g0) * tail)

def weights(curve, L):
    """(log n, Lambda(n) c_n / sqrt n) als Arb-Bälle für n = p^k < e^{2L}."""
    out = []
    X = 2 * float(fq(L)) + 1e-9
    N = curve["N"]
    for p, ap in zip(curve["primes"], curve["ap"]):
        lp = arb(p).log()
        if float(lp.lower()) >= X: break
        s1 = arb(ap) / arb(p).sqrt()
        s_prev, s_cur, k = arb(2), s1, 1
        while True:
            ln = k * lp
            if float(ln.lower()) >= 2 * float(fq(L)): break
            if float(ln.upper()) >= 2 * float(fq(L)):
                raise ValueError("Primzahlpotenz zu nah am Fensterrand – anderes L wählen")
            ck = s1 ** k if N % p == 0 else s_cur
            out.append((ln, lp * ck / arb(p) ** (arb(k) / 2)))
            if N % p != 0: s_prev, s_cur = s_cur, s1 * s_cur - s_prev
            k += 1
    return out

def certify(curve, L, coeffs):
    g = g_poly(coeffs, L)
    gc = [fmpq(c) for c in g.coeffs()]
    g0 = arb(g(fmpq(0)))
    cond = arb(curve["N"]).log() - 2 * (2 * arb.pi()).log()
    W = g0 * cond + arch_term(g, L)
    for ln, w in weights(curve, L):
        W -= 2 * w * horner_arb(gc, ln)
    f = sum((fq(c) * fmpq_poly.legendre_p(2 * k) for k, c in enumerate(coeffs)), fmpq_poly(0))
    F = f.integral()
    fhat0 = fq(L) * (F(fmpq(1)) - F(fmpq(-1)))
    return W / arb(fhat0) ** 2

if __name__ == "__main__":
    curves = {json.loads(l)["label"]: json.loads(l) for l in open(sys.argv[1])}
    for line in open(sys.argv[2]):
        d = json.loads(line); cv = curves[d["label"]]
        for key in ["std_r1", "radix_r1"]:
            h = d[key]
            if h is None: continue
            b = certify(cv, repr(h["l"]), h["coeffs"])
            target = d["rank"] + 1
            ok = b.upper() < target
            print(f"{d['label']:13s} {key:9s} L={h['l']:<5} Gleitkomma {h['bound']:.10f}  Arb {b.str(12, radius=True)}  "
                  f"< {target}: {'BEWIESEN' if ok else 'nicht bewiesen'}")
