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
2-Selmer/class-group method of rankbound. The explicit formula does not care about torsion. (For curves with
2-torsion, descent via a 2-isogeny, e.g. mwrank, is the classical alternative and can give unconditional bounds;
curve#40 is stated to have rank exactly 20.)

Result (25 Oct 2026, p < 10^8, certified with Arb): the analytic rank is certified to be at most the rank lower
bound on all 44 curves (rank 20: 34 curves, rank 21: 8, rank 22: 2), with parity, under GRH. The 44 curves form 22
pairs of 2-isogenous curves (isogeny classes proven with PARI `ellisomat`, recorded as `isogenous_to`); isogenous
curves have the same L-function, so each pair shares one certificate. Bober's fixed kernel also reaches all of them
and needs 1.8–2.9 times as many primes (median 2.4). Table and per-curve results: `highrank/torsion2/`.

## 3. Comparison with rankbound

```bash
python3 scripts/compare_rankbound.py
```

Curve-by-curve table against the rankbound database (R. Jain, commit pinned in the script), with the curves
where the two methods differ and the 2-torsion curves outside the scope of rankbound. Output: `work/compare/`.
The two methods certify different statements: radix-rank bounds the analytic rank under GRH for L(E,s);
rankbound bounds the 2-Selmer rank, hence the Mordell–Weil rank, under GRH for Dedekind zeta functions and
without BSD.

## 4. Upper bounds on the leaderboard

```bash
python3 scripts/icarm_notes.py --collect          # board certificates -> highrank/board/curves/, then commit + push
python3 scripts/icarm_notes.py                    # dry run: work/icarm_notes/preview.md, nothing is posted
python3 scripts/icarm_notes.py --apply --limit 3  # first three, check them on the website
python3 scripts/icarm_notes.py --apply            # the rest
```

Appends one note per curve to the ICARM commentary, below the existing text, in the form
"r_an <= r subject to GRH ... With the r independent points the rank is exactly r subject to GRH and BSD", with a
link to the curve's certificate file. Skipped automatically: curves whose leaderboard rank differs from the
certified rank, curves that already carry a radix-rank note, and curves (or their isogenous partner) whose
commentary already states an upper bound or the exact rank. Every old commentary is saved before posting and checked
after posting; `--restore <id>` puts it back. `--fix-4239` corrects "r_an = 32" to "r_an <= 32" in the note on #4239.
