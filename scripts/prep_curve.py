"""
Prepare an ICARM leaderboard curve for the streamed pipeline (ap -> radix-stream -> cert_rap).

Usage: python3 prep_curve.py <icarm.json> <out_dir>
Writes <out_dir>/curve.json (label, N, logN, bad, rank, c4, c6, root_number, target) and
<out_dir>/small_ap.txt (a_p for p < 1000 and all bad primes below 2^63, via PARI).
Checks: the conductor equals the product of local conductors (PARI elllocalred), the model's c4/c6 match.
Target: with root number (-1)^rank, a bound < rank + 2 already gives r_an = rank (parity). Otherwise r_an has the
other parity, so r_an >= rank + 1, and a bound < rank + 3 gives r_an = rank + 1 (reported as a warning).
"""
import json, math, os, sys
import cypari2

pari = cypari2.Pari()
pari.allocatemem(2 * 10**9, silent=True)
c = json.load(open(sys.argv[1]))
out = sys.argv[2]
os.makedirs(out, exist_ok=True)
E = pari.ellinit([int(a) for a in c['ainvs']])
N = int(c['conductor'])
bad = [int(p) for p in c['bad_primes']]
prod, w = 1, -1
for p in bad:
    prod *= p ** int(pari.elllocalred(E, p)[0])
    w *= int(pari.ellrootno(E, p))
assert prod == N, "conductor does not match the local data"
c4, c6 = int(E[9]), int(E[10])
if 'curve_key' in c:
    k4, k6 = map(int, c['curve_key'].split(':'))
    assert (k4, k6) == (c4, c6), "c4/c6 of the model differ from the leaderboard key"
r = int(c['rank_lower_bound'])
parity = (w == (-1) ** r)
target = r + 2 if parity else r + 3
label = f"icarm{c['id']}"
cv = dict(label=label, N=str(N), logN=math.log(N), bad=[p for p in bad if p < 2**63], rank=r,
          c4=str(c4), c6=str(c6), root_number=w, target=target)
json.dump(cv, open(os.path.join(out, 'curve.json'), 'w'), indent=1)
ps = sorted(set([int(p) for p in pari.primes([2, 999])] + [p for p in bad if p < 2**63]))
with open(os.path.join(out, 'small_ap.txt'), 'w') as f:
    for p in ps:
        f.write(f"{p} {int(pari.ellap(E, p))}\n")
print(f"{label}: rank >= {r}, ln N = {math.log(N):.4f}, root number {w:+d} "
      f"({'parity: bound < ' + str(target) + ' gives r_an = ' + str(r) if parity else 'WARNING: parity says r_an >= ' + str(r + 1) + '; bound < ' + str(target) + ' gives r_an = ' + str(r + 1)})")
