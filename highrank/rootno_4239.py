import json, math, cypari2
pari = cypari2.Pari(); pari.allocatemem(2*10**9, silent=True)
c = json.load(open('c4239.json'))
E = pari.ellinit([int(a) for a in c['ainvs']])
N = int(c['conductor']); bad = [int(p) for p in c['bad_primes']]
prod, w = 1, -1
for p in bad:
    loc = pari.elllocalred(E, p); prod *= p ** int(loc[0]); w *= int(pari.ellrootno(E, p))
D = int(pari.member_disc(E)) if hasattr(pari,'member_disc') else int(E[11])
print('minimal model check: disc matches page:', abs(D)==int(c.get('discriminant',D)) )
print('ln N =', math.log(N), 'conductor matches local data:', prod == N, 'root number', w, 'parity ok:', w == (-1)**c['rank_lower_bound'])
print('bad primes:', [p if p<10**6 else f'<{len(str(p))}-digit>' for p in bad])
