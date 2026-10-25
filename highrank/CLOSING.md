# Closing the high-rank board

The board run (cap p < 10^10) left 16 of 395 curves where the Radix kernel had not yet reached the target.
On all 16 the Radix bound was already below Bober's at the same primes; only primes were missing.

## 1. Close the 16 open curves

```bash
python3 scripts/close_board.py
```

Takes the open curves from `highrank/board/summary.csv`, smallest remaining gap first, and carries each up the
ladder p < 2·10^10, 3·10^10, 5·10^10, 8·10^10, 1.2·10^11 until the Radix bound is below the target. a_p are
extended block by block (nothing is recomputed between stages). Every hit is certified (Rust interval prime sum and
Arb), and Bober is certified on the same primes. Resumable; results in `work/close/summary.md`.
Disk: about 6 bytes per prime (25 GB at 10^11), freed after each curve unless `--keep-ap`.

## 2. Curves with rational 2-torsion

```bash
python3 scripts/batch_icarm.py --ids <ids from compare.md> --cap 1e9 --work work
```

The 44 leaderboard curves of rank ≥ 20 with a rational point of order 2 are outside the scope of the
2-Selmer/class-group method of rankbound. The explicit formula does not care about torsion.

Result (25 Oct 2026, p < 10^8, certified with Arb): the exact analytic rank is certified on all 44
(rank 20: 34 curves, rank 21: 8, rank 22: 2), with parity, under GRH. Bober's fixed kernel also reaches all 44 and
needs 1.8–2.9 times as many primes (median 2.4). Table and per-curve results: `highrank/torsion2/`.

## 3. Comparison with rankbound

```bash
python3 scripts/compare_rankbound.py
```

Curve-by-curve table against the rankbound database (R. Jain, commit pinned in the script), with the curves
where the two methods differ and the 2-torsion curves outside the scope of rankbound. Output: `work/compare/`.
The two methods certify different statements: radix-rank bounds the analytic rank under GRH for L(E,s);
rankbound bounds the 2-Selmer rank, hence the Mordell–Weil rank, under GRH for Dedekind zeta functions and
without BSD.
