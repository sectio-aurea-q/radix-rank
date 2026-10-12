"""Append a_p for p in [start, PMAX) to raw little-endian files icarm{id}.p.u64 / .ap.i64, in chunks (resumable)."""
import json, os, sys, time, numpy as np, cypari2
pari = cypari2.Pari(); pari.allocatemem(3 * 10**9, silent=True)
cid, PMAX = int(sys.argv[1]), int(float(sys.argv[2])); STEP = 25_000_000
c = {x['id']: x for x in json.load(open('database.json'))['curves']}[cid]
E = pari.ellinit([int(a) for a in c['ainvs']])
fp, fa, done = f"icarm{cid}.p.u64", f"icarm{cid}.ap.i64", f"icarm{cid}.done"
lo = int(open(done).read()) if os.path.exists(done) else 2
while lo < PMAX:
    hi = min(lo + STEP, PMAX); t = time.time()
    ps = [int(p) for p in pari.primes([lo, hi - 1])]
    ap = [int(pari.ellap(E, p)) for p in ps]
    with open(fp, 'ab') as f: f.write(np.array(ps, dtype='<u8').tobytes())
    with open(fa, 'ab') as f: f.write(np.array(ap, dtype='<i8').tobytes())
    open(done, 'w').write(str(hi)); lo = hi
    print(cid, hi, len(ps), f"{time.time()-t:.0f}s", flush=True)
