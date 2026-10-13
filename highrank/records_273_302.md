# The 2026 record curves #273 (rank ≥ 30) and #302 (rank ≥ 31)

Curves found by Claude with L. Alpöge and A. Howell (ICARM leaderboard, Aug 2026). For #273 the leaderboard commentary
states that Bober's method with Δ = 4.25 gives an analytic-rank bound below 32, i.e. primes up to e^{8.5π} ≈ 3.9·10¹¹
(≈ 1.5·10¹⁰ primes); with root number +1 this gives rank 30 under GRH and BSD.

Conductors recomputed from local data (PARI); root numbers +1 (#273) and −1 (#302) = (−1)^r. a_p for all p < 10⁹.

| curve | r | ln N | target | Radix | Bober (fixed kernel) |
|---|---|---|---|---|---|
| #273 | 30 | 339.35 | < 32 | **13 461 681 primes (p < 2.4·10⁸), L = 9.66, Arb [31.99478129 ± 4.5e-9]** | not reached with all 50 847 534 primes < 10⁹: Arb [32.08095018 ± 1.5e-9] at L = 10.36 |
| #302 | 31 | 375.22 | < 33 | not reached with p < 10⁹ (33.566 at L = 10.36) | not reached: Arb [34.02864458 ± 1.3e-9] at L = 10.36 |

For #273 the optimal kernel therefore needs fewer than a quarter of the primes the fixed kernel needs at the very least
(certified), and about 1/1100 of the primes used with Δ = 4.25.
