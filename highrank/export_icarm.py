"""a_p for an ICARM leaderboard curve (global minimal model) for p < PMAX, JSON line for radix-rank."""
import json, math, sys, time, cypari2
pari = cypari2.Pari(); pari.allocatemem(2 * 10**9, silent=True)
cid, PMAX = int(sys.argv[1]), int(float(sys.argv[2]))
d = json.load(open('database.json'))
c = {x['id']: x for x in d['curves']}[cid]
E = pari.ellinit([int(a) for a in c['ainvs']])
N = int(c['conductor']); bad = sorted(int(p) for p in c['bad_primes'])
assert math.prod(bad) and all(N % p == 0 for p in bad)
t = time.time()
primes = [int(p) for p in pari.primes([2, PMAX])]
ap = [int(pari.ellap(E, p)) for p in primes]
out = dict(label=f"icarm{cid}", logN=math.log(N), N_str=str(N), rank=c['rank_lower_bound'], a=c['ainvs'],
           bad=[p for p in bad if p < 2**63], primes=primes, ap=ap)
json.dump(out, open(f"icarm{cid}_{sys.argv[2]}.json", "w"))
print(cid, c['rank_lower_bound'], len(primes), f"{time.time()-t:.0f}s", flush=True)
