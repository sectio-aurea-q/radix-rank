import json, sys, time
sys.path.insert(0, 'radix-rank/scripts')
import certify as C
cid = int(sys.argv[1]); keys = sys.argv[2].split(',')
cv = json.load(open(f'icarm{cid}_3e7.json')); cv['N'] = int(cv['N_str'])
d = json.loads(open(f'out{cid}_3e7.jsonl').readline())
for key in keys:
    h = d[key]
    if h is None: print(cid, key, 'no hit'); continue
    t = time.time(); b = C.certify(cv, repr(h['l']), h['coeffs'])
    target = d['rank'] + (1 if key.endswith('r1') else 2)
    print(f"#{cid} {key:9s} L={h['l']:<5} primes={h['primes']:>7}  float {h['bound']:.6f}  Arb {b.str(10, radius=True)}  "
          f"< {target}: {'PROVEN' if b.upper() < target else 'NOT proven'}  ({time.time()-t:.0f}s)", flush=True)
