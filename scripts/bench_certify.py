"""Alle Treffer des Rust-Kerns beweisfest (Arb) prüfen und auswerten."""
import json, sys
from multiprocessing import Pool
from certify import certify

curves = {}
for l in open("curves.jsonl"):
    d = json.loads(l); curves[d["label"]] = d
results = [json.loads(l) for l in open("results.jsonl")]
KEYS = [("std_r1", 1), ("radix_r1", 1), ("std_r2", 2), ("radix_r2", 2)]

def work(d):
    cv = curves[d["label"]]; out = dict(label=d["label"], N=d["N"], rank=d["rank"])
    for key, off in KEYS:
        h = d[key]
        if h is None:
            out[key] = None; continue
        b = certify(cv, repr(h["l"]), h["coeffs"])
        out[key] = dict(l=h["l"], primes=h["primes"], lo=float(b.lower()), hi=float(b.upper()),
                        proven=bool(b.upper() < d["rank"] + off),
                        sane=bool(b.upper() >= d["rank"]))   # nie unter dem wahren Rang
    return out

if __name__ == "__main__":
    with Pool() as pool:
        cert = pool.map(work, results, chunksize=8)
    json.dump(cert, open("certified.json", "w"))
    bad = [c for c in cert for k, _ in KEYS if c[k] and not c[k]["sane"]]
    print("Kurven:", len(cert), "| Zertifikate unter dem wahren Rang (Fehler!):", len(bad))
    for k, _ in KEYS:
        n_hit = sum(1 for c in cert if c[k]); n_pr = sum(1 for c in cert if c[k] and c[k]["proven"])
        print(f"{k:9s}: Treffer {n_hit}, davon beweisfest {n_pr}")
