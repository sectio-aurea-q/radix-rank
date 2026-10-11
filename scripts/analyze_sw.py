import json, numpy as np
cert = json.load(open("sw_certified.json"))
grp = {json.loads(l)["label"]: json.loads(l)["group"] for l in open("sw.jsonl")} if False else None
groups = {}
for l in open("sw_out.jsonl"):
    d = json.loads(l); groups[d["label"]] = "prim" if d["label"].startswith("SW") and d["N"] > 10**9 else "allg"
need = lambda c, k: c[k]["primes"] if c[k] and c[k]["proven"] else np.inf
res = {}
for crit, ks, kr in [("r+1", "std_r1", "radix_r1"), ("r+2", "std_r2", "radix_r2")]:
    S = np.array([need(c, ks) for c in cert]); R = np.array([need(c, kr) for c in cert])
    both = np.isfinite(S) & np.isfinite(R)
    print(f"\n== {crit} == bewiesen (Budget 9592 Primzahlen): Standard {np.isfinite(S).sum()}  Radix {np.isfinite(R).sum()} von {len(cert)}")
    print(f"   Radix weniger {np.sum(R<S)}, gleich {np.sum(R==S)}, mehr {np.sum(R>S)} | Ersparnis gesamt {100*(1-R[both].sum()/S[both].sum()):.1f} %")
    tab = []
    for g, lab in [("allg", "N 4·10⁶–10⁸"), ("prim", "N prim 3,6–7,5·10⁹")]:
        for r in range(5):
            m = np.array([groups[c["label"]] == g and min(c["rank"], 4) == r for c in cert])
            mb = m & both
            row = dict(group=lab, rank=r, n=int(m.sum()), std_ok=int(np.isfinite(S[m]).sum()), radix_ok=int(np.isfinite(R[m]).sum()),
                       med_std=float(np.median(S[mb])) if mb.any() else None, med_radix=float(np.median(R[mb])) if mb.any() else None,
                       saving=float(100*(1-R[mb].sum()/S[mb].sum())) if mb.any() else None)
            tab.append(row)
            print(f"   {lab:20s} Rang {r}{'+' if r==4 else ' '}: bewiesen Std {row['std_ok']:3d} / Radix {row['radix_ok']:3d} von {row['n']}  "
                  f"Median {row['med_std']:.0f} vs {row['med_radix']:.0f}  Ersparnis {row['saving']:.1f} %")
    res[crit] = dict(table=tab, total=len(cert), std=int(np.isfinite(S).sum()), radix=int(np.isfinite(R).sum()),
                     saving=float(100*(1-R[both].sum()/S[both].sum())), fewer=int(np.sum(R<S)), more=int(np.sum(R>S)))
json.dump(res, open("sw_summary.json", "w"), indent=1)
