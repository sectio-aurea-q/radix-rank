"""Benchmark-Datensatz: eine Kurve pro Isogenieklasse aus Cremonas Tabellen,
geschichtet nach Rang; a_p per PARI. Ausgabe: JSON-Zeilen für den Rust-Kern."""
import json, random, sys
import numpy as np
import cypari2

pari = cypari2.Pari(); pari.allocatemem(2 * 10**9, silent=True)
PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8200          # e^9 ≈ 8103 -> L <= 4.5
PER = int(sys.argv[2]) if len(sys.argv) > 2 else 150            # Kurven pro (Bereich, Rang)
random.seed(20261010)
files = ["data/allcurves.00000-09999", "data/allcurves.100000-109999", "data/allcurves.400000-409999"]
primes = [int(p) for p in pari.primes([2, PMAX])]
out = open("curves.jsonl", "w"); count = {}
for fn in files:
    byrank = {}
    for line in open(fn):
        N, iso, num, ai, r, t = line.split()
        if num != "1": continue
        byrank.setdefault(int(r), []).append((f"{N}{iso}{num}", int(N), json.loads(ai), int(r)))
    for r, lst in sorted(byrank.items()):
        pick = lst if len(lst) <= PER else random.sample(lst, PER)
        for label, N, ai, rank in pick:
            E = pari.ellinit(ai)
            an = pari.ellan(E, PMAX)
            ap = [int(an[p - 1]) for p in primes]
            out.write(json.dumps(dict(label=label, N=N, rank=rank, a=ai, primes=primes, ap=ap)) + "\n")
            count[(fn[-11:], r)] = count.get((fn[-11:], r), 0) + 1
out.close()
for k, v in sorted(count.items()): print(k, v)
print("gesamt", sum(count.values()))
