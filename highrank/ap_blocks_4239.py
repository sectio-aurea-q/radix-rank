"""a_p for curve 4239 on blocks [k*STEP,(k+1)*STEP), block files ap/b{k}.p/.a (u64/i64 LE). args: worker nworkers kmax"""
import json, os, sys, numpy as np, cypari2
pari = cypari2.Pari(); pari.allocatemem(2*10**9, silent=True)
STEP = 10_000_000
w, nw, kmax = map(int, sys.argv[1:4])
E = pari.ellinit([int(a) for a in json.load(open('c4239.json'))['ainvs']])
for k in range(w, kmax, nw):
    fp = f'ap/b{k:03d}.p'
    if os.path.exists(fp + '.ok'): continue
    lo, hi = max(2, k*STEP), (k+1)*STEP - 1
    ps = [int(p) for p in pari.primes([lo, hi])]
    ap = [int(pari.ellap(E, p)) for p in ps]
    np.array(ps, dtype='<u8').tofile(fp); np.array(ap, dtype='<i8').tofile(f'ap/b{k:03d}.a')
    open(fp + '.ok', 'w').write(str(len(ps))); print(k, len(ps), flush=True)
