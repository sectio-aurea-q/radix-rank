import json, numpy as np
cert = json.load(open("certified.json"))
def need(c, k): return c[k]["primes"] if c[k] and c[k]["proven"] else np.inf
rng = lambda N: "N<10^4" if N < 10**4 else ("N≈1,0·10^5" if N < 2*10**5 else "N≈4,0·10^5")
out = {}
for crit, ks, kr in [("r+1", "std_r1", "radix_r1"), ("r+2 (Parität)", "std_r2", "radix_r2")]:
    S = np.array([need(c, ks) for c in cert]); R = np.array([need(c, kr) for c in cert])
    both = np.isfinite(S) & np.isfinite(R)
    ratio = R[both] / S[both]
    print(f"\n== Kriterium {crit} ==")
    print(f"bewiesen bis 1028 Primzahlen (L<=4.5): Standard {np.isfinite(S).sum()}  Radix {np.isfinite(R).sum()}  von {len(cert)}")
    # WICHTIG: nur ueber die beidseitig entschiedenen Kurven zaehlen.
    # need() liefert np.inf fuer unentschieden, und np.inf == np.inf ist True --
    # eine Zaehlung ueber alle Kurven wuerde die beidseitig UNentschiedenen als
    # "gleich viele" mitzaehlen und die nur-von-einem-Kern entschiedenen verstecken.
    print(f"Radix braucht weniger: {np.sum(R[both]<S[both])}, gleich viele: "
          f"{np.sum(R[both]==S[both])}, mehr: {np.sum(R[both]>S[both])}  "
          f"(von {both.sum()} beidseitig entschiedenen)")
    print(f"nur Radix entscheidet: {np.sum(np.isfinite(R) & ~np.isfinite(S))}, "
          f"nur Standard: {np.sum(np.isfinite(S) & ~np.isfinite(R))}, "
          f"keiner: {np.sum(~np.isfinite(S) & ~np.isfinite(R))}")
    print(f"Primzahlen gesamt (beide bewiesen): Standard {S[both].sum():.0f}  Radix {R[both].sum():.0f}  -> Ersparnis {100*(1-R[both].sum()/S[both].sum()):.1f} %")
    print(f"Median Verhältnis Radix/Standard: {np.median(ratio):.3f}")
    for lab in ["N<10^4", "N≈1,0·10^5", "N≈4,0·10^5"]:
        for r in range(4):
            m = np.array([rng(c["N"]) == lab and c["rank"] == r for c in cert]) & both
            if m.sum() < 5: continue
            print(f"   {lab:11s} Rang {r}: n={m.sum():3d}  Primzahlen Median Standard {np.median(S[m]):6.0f}  Radix {np.median(R[m]):6.0f}  Ersparnis {100*(1-R[m].sum()/S[m].sum()):5.1f} %")
    budgets = [5,10,20,50,100,200,500,1028]
    out[crit] = dict(budgets=budgets, std=[int(np.sum(S<=b)) for b in budgets], radix=[int(np.sum(R<=b)) for b in budgets], total=len(cert))
    print("Budget -> bewiesen (Standard / Radix):", [(b, int(np.sum(S<=b)), int(np.sum(R<=b))) for b in budgets])
json.dump(out, open("summary.json","w"), indent=1)
