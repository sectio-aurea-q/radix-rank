# ICARM #4239 (rank ≥ 32): Radix vs Bober

Curve: ICARM leaderboard #4239, rank ≥ 32, submitted 15 Oct 2026 by Z. Klagsbrun (Elkies–Klagsbrun–Li; specialization of a
rank-17 fibration of the K3 surface of discriminant 1152). Smallest conductor of any curve of rank ≥ 32 on the board.

Conductor recomputed from local data (PARI `elllocalred` at all 15 bad primes): matches the leaderboard,
ln N = 414.0857. Root number +1 = (−1)^32, so with parity a bound < 34 gives analytic rank 32 (under GRH).

a_p via PARI for all primes p < 10⁹ (50 847 534 primes). K = 10, same explicit formula, archimedean term and Arb certification as
`records_273_302.md`. "Bober" = fixed Fejér kernel, "Radix" = optimal f in span{P₀,…,P₁₈}(x/L).

| | L | primes used | certified bound (Arb) |
|---|---|---|---|
| Bober, all p < 10⁹ | 10.36 | 50 690 199 | [35.67841612 ± 2.8e-9] |
| Radix, all p < 10⁹ | 10.36 | 50 690 199 | **[35.30281011 ± 1.2e-9]** |
| Radix, matching Bober's p < 10⁹ value | 9.94 | **22 863 324** | [35.67523812 ± 3.1e-9] |

Target (< 34) is **not reached** with p < 10⁹ by either kernel, as for #302.

* At equal primes the Radix bound is 0.376 lower.
* Radix reaches Bober's best value (all primes < 10⁹) with 22.9 M primes (p < e^19.88 ≈ 4.3·10⁸): **55 % fewer primes (factor 2.2)**, window shorter by ΔL = 0.42.
* This is a smaller gain than at ranks 28–31 (factor ≈ 6.5 at the target r + 2). There the comparison is made where Bober's bound flattens just
  above the target; here both bounds are still ≈ 1.3–1.7 above it and decreasing at the same rate (≈ 0.85 per unit of L at L = 10.36).
* Rough linear extrapolation of the certified trajectory (not certified): bound < 34 needs L ≳ 11.9 for Radix (p ≳ 2·10¹⁰) and L ≳ 12.3
  for Bober (p ≳ 5·10¹⁰); the true requirement is likely larger since the slope decreases.

Reproduce: `rootno_4239.py` (conductor + root number), `ap_blocks_4239.py` (a_p in 10⁷ blocks, parallel workers), the optimiser with
`RADIX_STREAM=1 RADIX_LAST=1` (new switch: also report the last window of the grid, for endpoint certification) and `cert_4239.py`
(streamed Arb certification).
