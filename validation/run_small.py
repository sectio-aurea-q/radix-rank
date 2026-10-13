import json, math, cypari2, sys
from bober_formula import bober_bound
pari = cypari2.Pari(); pari.allocatemem(2*10**9, silent=True)
table = {'E20':(2.0,21.70),'E21':(2.5,22.68),'E22':(2.0,23.71),'E23':(2.5,24.49),'E24':(2.5,25.57)}
PMAX = 6_700_000
primes = [int(p) for p in pari.primes([2, PMAX])]
for n,(D,pub) in table.items():
    c = json.load(open(n+'.json')); E = pari.ellinit(c['a'])
    ap = [int(pari.ellap(E,p)) for p in primes]; bad = set(int(q) for q in c['bad'] if int(q) < 2**63)
    N = int(c['N']); logN = math.log(N)
    b, npr = bober_bound(logN, primes, ap, bad, D)
    print(f"{n}  Delta={D}  published {pub:.2f}   independent formula {b:.6f}  (rounded up {math.ceil(b*100)/100:.2f})  primes={npr}", flush=True)
    json.dump(dict(label=n, logN=logN, N_str=c['N'], rank=int(n[1:]), a=[str(x) for x in c['a']], bad=sorted(bad), primes=primes, ap=ap), open(n+'_ap.json','w'))
