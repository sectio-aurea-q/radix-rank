# Kapsel — search rooms for high-rank specializations

A small toolkit for searching elliptic fibrations over Q(T) for specializations of high rank,
with the Radix bound as the upper side.

| Room | Tool | Status |
|---|---|---|
| 1 · Scan | `kapsel-scan` (Rust) | working, validated |
| 2 · Points (lower bound) | needs the generic sections of the fibration | open |
| 3 · Radix (upper bound) | `radix-rank`, `radix-stream`, `sieve-score` | working |

## Room 1: scan

Families are read from the ICARM leaderboard comments. Each family is checked: its stated
specialization must reproduce the listed curve (same j-invariant).

```bash
curl -s -o kapsel/database.json https://elliptic-rank.icarm.cloud/database.json
python3 kapsel/families_from_icarm.py kapsel/database.json kapsel/families
cargo build --release
./target/release/kapsel_scan kapsel/families/icarm4239.json out.jsonl <umax> <vmin> <vmax> [B1=1000] [B2=3000] [N=20000] [K=500] [threads]
```

For every prime 3 <= p < B2 the scanner tabulates a_p(t) for all t in F_p once (point count of
y^2 + (Lx+B)y = x(x+D)(x+E) over F_p). A parameter T = u/v then costs one table lookup per prime.
Stage 1 ranks all T with |u| <= umax, vmin <= v <= vmax, gcd(u,v) = 1 by the Mestre–Nagao sum
S1 over p < B1 and keeps the best N; stage 2 rescores them with S1, S2 over p < B2 and writes the best K.

Validation: S1 and S2 agree with PARI `ellap` on all 1959 parameters of a test box (|u|, v <= 40,
p < 500) to float precision. Throughput: about 0.9 million parameters per second per core.

## Ranking test (25 Oct 2026)

Seven leaderboard families; in each, the listed record specialization (rank 24–32) mixed with 120
random specializations of similar height, a_p for p < 10^4. Position of the record among 121
(sum over the seven families, 7 = perfect):

| score | sum of positions |
|---|---|
| Nagao S1 | 7 |
| Nagao S2 | 7 |
| fixed Fejér kernel, conductor term removed, L = 4.6 | 17 |
| Radix kernel, conductor term removed, L = 4.6 | 77 |
| Radix bound itself, L = 4.6 | 686 |

The Nagao sums find; the Radix kernel does not rank better. Its use is the upper side: excluding a
candidate with a certified bound before an expensive descent, and closing the exact rank of a hit.

## Room 2: what is missing

PARI `ellrank` does not finish within minutes on these curves (|Delta| ~ 10^118). Record searches
start from the specialized generic sections (17 for the X1152 fibration) and add points by descent.
The section basis of the fibration is needed before room 2 can run.
