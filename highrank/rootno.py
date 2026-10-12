import json, cypari2
pari = cypari2.Pari(); pari.allocatemem(2*10**9, silent=True)
cs = {x['id']: x for x in json.load(open('database.json'))['curves']}
for cid in (3762, 3835, 3897, 3900):
    c = cs[cid]; E = pari.ellinit([int(a) for a in c['ainvs']])
    N = int(c['conductor']); bad = [int(p) for p in c['bad_primes']]
    # conductor check: product of p^f_p from local data
    prod, w = 1, -1
    for p in bad:
        loc = pari.elllocalred(E, p); prod *= p ** int(loc[0]); w *= int(pari.ellrootno(E, p))
    print(cid, 'rank>=', c['rank_lower_bound'], 'conductor matches:', prod == N, 'root number', w,
          'parity ok:', w == (-1) ** c['rank_lower_bound'])
