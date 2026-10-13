"""Independent implementation of Bober's bound (Bober 2013, eq. 3; archimedean term in Spicer's closed form),
written directly from the papers, sharing no code with radix-rank.
   sum_gamma f(gamma;D) = -eta/(pi D) + (pi^2/6 - Li2(exp(-2 pi D)))/(2 pi^2 D^2) + log(sqrt(N)/(2 pi))/(pi D)
                          - 1/(pi D) * sum_{p^k <= exp(2 pi D)} log p * c(p^k) / p^{k/2} * (1 - k log p / (2 pi D))
   c(p^k) = alpha^k + beta^k with alpha + beta = a_p, alpha beta = p (good p);  a_p^k (bad p)."""
import math, mpmath
mpmath.mp.dps = 30
def bober_bound(logN, primes, ap, bad, D):
    X = 2 * math.pi * D
    pr = mpmath.mpf(0)
    terms = []
    for p, a in zip(primes, ap):
        lp = math.log(p)
        if lp > X: break
        k, c_prev, c_cur = 1, 2, a          # integer power sums, un-normalised
        pk = p
        while k * lp <= X:
            c = (a ** k) if p in bad else c_cur
            terms.append(lp * c / pk * (1 - k * lp / X))   # (normalised c / p^{k/2}) / p^{k/2} = c / p^k
            if p not in bad:
                c_prev, c_cur = c_cur, a * c_cur - p * c_prev
            k += 1; pk *= p
    S = math.fsum(terms)
    piD = mpmath.pi * D
    arch = -mpmath.euler / piD + (mpmath.pi ** 2 / 6 - mpmath.polylog(2, mpmath.e ** (-2 * piD))) / (2 * mpmath.pi ** 2 * D ** 2)
    cond = (mpmath.mpf(logN) / 2 - mpmath.log(2 * mpmath.pi)) / piD
    return float(arch + cond - S / piD), len([p for p in primes if math.log(p) <= X])
