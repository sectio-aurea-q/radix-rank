# Validation against published values

Is the fixed-kernel baseline in radix-rank really Bober's method? We re-implemented Bober's bound
directly from the papers (`bober_formula.py`: Bober 2013 eq. 3, archimedean term in Spicer's closed form,
no code shared with radix-rank) and reproduced every published value we could find:

| curve | Δ | published | independent | radix-rank fixed kernel (L = πΔ) |
|---|---|---|---|---|
| E20 | 2.0 | 21.70 (Bober 2013) | 21.6907 | 21.6907 |
| E21 | 2.5 | 22.68 | 22.6727 | 22.6727 |
| E22 | 2.0 | 23.71 | 23.7048 | 23.7048 |
| E23 | 2.5 | 24.49 | 24.4834 | 24.4834 |
| E24 | 2.5 | 25.57 | 25.5682 | 25.5682 |
| E27 | 3.1 | < 29 (Klagsbrun–Sherman–Weigandt) | 28.9881 | 28.9881 |
| E28 | 3.2 | 31.30 (Bober 2013) | 31.2984 | 31.2984 |

(Bober's table rounds up.) The Rust float values agree with the independent formula to 2·10⁻⁷, the Arb certification to 10⁻¹³.
Conductors are taken from the ICARM leaderboard; log N agrees with Bober's table to all printed digits.

## Radix vs Bober on the historic record curves (all certified with Arb)

| curve | primes for bound < r+2: Bober | Radix | bound < r+1: Radix |
|---|---|---|---|
| E20 | 11 988 | **5 515** | 47 114 |
| E21 | 181 084 | **83 298** | — |
| E22 | 11 988 | **5 133** | 37 128 |
| E23 | 77 395 | **23 068** | 353 398 |
| E24 | 103 973 | **32 069** | 433 671 |
| E27 | 15 078 540 | **4 340 969** | — |

Bober's minimal window for E27 (L = 9.72, Δ = 3.09) matches the Δ = 3.1 used by Klagsbrun–Sherman–Weigandt.
For E28 neither kernel reaches < 30 with p < 10⁹ (L = 10.36: Bober 31.09, Radix 30.82); the gain is curve-dependent,
from a factor ≈ 2 (E21, E20) to ≈ 6.5 (ICARM curves, see `../highrank/`).
