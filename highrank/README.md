# High-rank test: radix-rank vs. Bober on the ICARM leaderboard

Question (Z. Klagsbrun, Oct 2026): on curves of large conductor, does the per-curve optimal test function
reduce the number of primes needed to certify the rank, compared with the fixed Fejér kernel of Bober (2013)?

**Curves:** for each rank 28–31, the curve with the smallest recorded conductor on the
[ICARM Elliptic Curve Rank Leaderboard](https://elliptic-rank.icarm.cloud) (snapshot 12 Oct 2026, 4201 curves).
Conductors and root numbers recomputed from local data with PARI; root number = (−1)^r in all four cases.

**Method:** same explicit formula, archimedean term and certification for both kernels. "Bober" = f constant on [−L, L]
(g = f∗f the Fejér/triangle kernel); "Radix" = optimal f in span{P₀,…,P₁₈}(x/L), which contains Bober's kernel.
a_p via PARI for all p < 10⁹ (50 847 534 primes). L-grid step 0.01. Every number below is certified with Arb.

| ICARM | rank | log N | bound < r+2 (with parity): Bober | Radix | bound < r+1: Bober | Radix |
|---|---|---|---|---|---|---|
| #3762 | 28 | 242.9 | 122 776 | **19 235** | 21 600 954 | **234 855** |
| #3835 | 29 | 265.4 | 316 068 | **48 890** | > 50 847 534 (30.0029 at p < 10⁹) | **717 477** |
| #3897 | 30 | 279.9 | 387 811 | **59 820** | > 50 847 534 (31.0898) | **691 123** |
| #3900 | 31 | 325.3 | 4 340 969 | **665 760** | > 50 847 534 (32.4306) | **13 209 786** |

With parity the Radix kernel needs 84–85 % fewer primes (factor ≈ 6.5) on all four curves; it reaches each
target with a window about ΔL ≈ 1 shorter, i.e. a prime bound smaller by roughly e². For the bound below r+1 the gap is
larger (factor 92 at rank 28), because the Fejér bound flattens just above the target.

Scripts: `export_icarm.py`, `ap_chunks.py` (a_p), `rootno.py` (conductor + root number), `cert_stream.py` (Arb certification,
streamed). The Rust core gained `logN`/`bad` fields (conductors beyond 64 bit), binary a_p input, and an O(N·m) Chebyshev-moment
prime sum; on the 7 reference curves it reproduces the published results exactly.
