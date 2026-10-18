"""Die Radix-Fläche: der optimale Kern f-hat(t)^2 über Höhe t und Fensterbreite L.

Für jede Fensterbreite L wird die Gram-Matrix W der Weil-Form in der Basis P_{2k}(x/L)
aus der zertifizierten Routine scripts/certify.py gewonnen -- per Polarisierung, also
ohne eine zweite Implementierung der Form:

    W_ij = ( B(e_i + e_j) - B(e_i) - B(e_j) ) / 2 ,

wobei B der unnormierte Zähler W_E(f) ist. Der optimale Kern ist c = W^{-1} v mit
v = (2L, 0, ..., 0), normiert auf f-hat(0) = 1, und die Fläche ist f-hat(t)^2 mit

    f-hat(t) = sum_k c_k * 2L * (-1)^k * j_{2k}(L t)     (sphärische Besselfunktion)

Die Nullstellen von L(E,s) kommen aus PARI. Ausgabe: paper/radix_surface.json
"""
import json, os, sys
import numpy as np
from scipy.special import spherical_jn
from scipy.linalg import solve
from flint import arb, fmpq, fmpq_poly

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import certify as C

K = 10


def numerator(curve, L, coeffs):
    """W_E(f) ohne Normierung, als float -- dieselbe Rechnung wie certify.certify()."""
    g = C.g_poly(coeffs, L)
    gc = [fmpq(c) for c in g.coeffs()]
    W = arb(g(fmpq(0))) * (arb(curve["N"]).log() - 2 * (2 * arb.pi()).log()) + C.arch_term(g, L)
    for ln, w in C.weights(curve, L):
        W -= 2 * w * C.horner_arb(gc, ln)
    return float(W)


def gram(curve, L):
    e = lambda i: [1 if j == i else 0 for j in range(K)]
    B = [numerator(curve, L, e(i)) for i in range(K)]
    W = np.empty((K, K))
    for i in range(K):
        W[i, i] = B[i]
        for j in range(i + 1, K):
            s = numerator(curve, L, [a + b for a, b in zip(e(i), e(j))])
            W[i, j] = W[j, i] = (s - B[i] - B[j]) / 2
    return W


def fhat(t, L, c):
    t = np.atleast_1d(np.asarray(t, dtype=float))
    out = np.zeros_like(t)
    for k, ck in enumerate(c):
        out += ck * 2 * L * (-1) ** k * spherical_jn(2 * k, L * t)
    return out


if __name__ == '__main__':
    label = sys.argv[1] if len(sys.argv) > 1 else '5077a1'
    curve = next(json.loads(l) for l in open(os.path.join(ROOT, 'data', 'reference_curves.jsonl'))
                 if json.loads(l)['label'] == label)
    Ls = np.round(np.linspace(0.6, 3.0, 41), 3)
    ts = np.round(np.linspace(0, 30, 151), 3)
    H, bounds = [], []
    for L in Ls:
        W = gram(curve, repr(float(L)))
        v = np.zeros(K); v[0] = 2 * L
        c = solve(W, v); c = c / (c @ v)          # f-hat(0) = 1
        H.append(np.round(fhat(ts, L, c) ** 2, 6).tolist())
        bounds.append(round(float(1 / (v @ solve(W, v))), 4))
        print(f"  L={L:5.3f}  bound={bounds[-1]:9.4f}", flush=True)
    try:
        import cypari2
        pari = cypari2.Pari(); pari.allocatemem(2 * 10 ** 9, silent=True)
        E = pari.ellinit(curve['a'])
        zs = sorted({round(float(z), 4) for z in pari.lfunzeros(pari.lfuncreate(E), 30) if float(z) > 1e-6})
    except Exception as exc:
        print("Nullstellen nicht berechnet:", exc); zs = []
    out = os.path.join(ROOT, 'paper', 'radix_surface.json')
    json.dump(dict(curve=label, rank=curve['rank'], N=curve['N'],
                   L=Ls.tolist(), t=ts.tolist(), h=H, zeros=zs, bound=bounds), open(out, 'w'))
    print("geschrieben:", out, " Nullstellen:", zs[:6])
