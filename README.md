# radix-rank

**Optimal windowed test functions for bounding analytic ranks of elliptic curves.**
Sabrina Hackfort · Independent researcher · ORCID [0009-0007-7458-2487](https://orcid.org/0009-0007-7458-2487)

> **Built on the Radix Riemanniana** (Hackfort, 2026): the reproducing kernel of the Weil form on a window [−L, L],
> constructed from the Gamma factor, the conductor and the primes p < e^{2L} alone, without the zeros:
>
> **ℛ_{E,L}(w, z) = Φ(w)\* G⁻¹ Φ(z)**,  and the optimal rank bound is **r_an(E) ≤ 1 / ℛ_{E,L}(0, 0)**  (under GRH).
>
> G is the Gram matrix of the Weil form of E, and Φ(z) = (b̂₀(z), …, b̂_{K−1}(z)) is the Fourier basis on the window.

Assuming GRH, the explicit formula bounds the analytic rank of an elliptic curve E/ℚ from the traces
a_p for p < e^{2L}, for any admissible test function supported on [−L, L]. The usual choice is a fixed
Fejér-type kernel (Bober 2013; Klagsbrun–Sherman–Weigandt 2019). `radix-rank` instead chooses, **for each
individual curve**, the test function that minimises the bound in a finite-dimensional space on the window.
This choice amounts to one linear solve with the Gram matrix of the Weil form of E.
The optimal bound is 1 / ℛ_{E,L}(0,0), the reciprocal of the radix Riemanniana at the central point.
Every bound is certified with ball arithmetic.

Paper: [`paper/radix_rank.pdf`](paper/radix_rank.pdf) · Archived on Zenodo: [doi:10.5281/zenodo.23299461](https://doi.org/10.5281/zenodo.23299461)

## Reproduce

Requirements: Rust (stable), Python ≥ 3.11, `pip install -r requirements.txt` (PARI via `cypari2`, Arb via `python-flint`).

```bash
bash reproduce.sh quick     # ~1 min: 7 reference curves (ranks 0–6): optimise + certify
bash reproduce.sh cremona   # ~10 min: full Cremona benchmark (1642 curves) from scratch
```

The `quick` mode should end with every line reading `BEWIESEN` (= proven), i.e. the certified upper
end of the bound lies below rank + 1.

## What is in here

| path | content |
|---|---|
| `src/` | Rust optimiser: window Gram matrices, Weil kernel, scan over L, streaming mode (`RADIX_STREAM=1`) |
| `scripts/certify.py` | rigorous certification (Arb ball arithmetic) of a given test function |
| `scripts/export*.py` | build curve data (a_p via PARI) from Cremona's tables / the Stein–Watkins database |
| `scripts/*certify*.py`, `analyze*.py` | benchmark certification and evaluation |
| `certificates/` | all certified results reported in the paper |
| `data/reference_curves.jsonl` | 7 reference curves 11a1 … 5187563742a1 with their a_p |

## Results (all certified, see paper)

| sample | curves | primes saved (aggregate) | worse than fixed kernel |
|---|---|---|---|
| Cremona, N < 4.1·10⁵ | 1642 | 8.3–15.2 % | 0 |
| Stein–Watkins, 4.4·10⁶ ≤ N ≤ 7.5·10⁹ | 2500 | 17.0–22.7 % | 0 |
| single curves, N up to 5·10¹⁵ | 9 | 0–43 % | 0 |

"Never worse" holds by construction: the fixed kernel lies in the optimised space.
The prime-side values agree with zero-side sums computed from the zeros of L(E,s). The gap halves with each doubling of the zero height.

## Record curves: streamed pipeline

For curves with huge conductors (the ICARM leaderboard) the bound needs a_p for 10⁹–10¹¹ primes. One command
does everything, uses all cores and can be resumed at any point:

```bash
scripts/rank_pipeline.sh 4239 2e10      # ICARM curve id, prime bound
```

| step | tool | what it does |
|---|---|---|
| 1 | `scripts/prep_curve.py` | conductor from local data (PARI), root number, target from parity |
| 2 | `ap` (Rust) | a_p for all p < PMAX: baby-step/giant-step over the whole Hasse interval, a_p mod 2 from the 2-torsion test, affine arithmetic with batched inversions; never guesses (aborts if a_p is not uniquely determined) |
| 3 | `radix_stream` (Rust) | streams the compact `.rap` blocks (6 bytes per prime), bisection on L for Radix and Bober |
| 4 | `scripts/cert_rap.py` | Arb certification of every hit and both endpoints, in parallel (exact ball transfer between processes) |

Checks: on #4239 `ap` reproduces all 50 847 534 PARI values of a_p for p < 10⁹ bit for bit, 2.4–2.8× faster than PARI per core;
`radix_stream` reproduces the bounds of `highrank/record_4239.md` to all printed digits.

## Limitations

The bounds are conditional on GRH and concern the *analytic* rank. The basis size is K = 10 (on #4239, K = 20 changes the
bound only in the fifth decimal); the L-grid step is 0.01–0.05. Certification runs in Python on top of Arb (parallel).

## License

MIT. If you use this, please cite via [`CITATION.cff`](CITATION.cff).

An AI assistant (Claude, Anthropic) was used in preparing this work.
