"""Stichprobe aus der Stein-Watkins-Datenbank (Watkins' analytische Ränge, numerisch, unbewiesen).
Pro Kurve: Führer (PARI-Kontrolle), Wurzelzahl, a_p bis PMAX. Ausgabe: JSON-Zeilen für radix-rank."""
import bz2, glob, json, random, sys
import cypari2

pari = cypari2.Pari(); pari.allocatemem(3 * 10**9, silent=True)
PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
PER = int(sys.argv[2]) if len(sys.argv) > 2 else 250
random.seed(20261010)
pool = {}
for fn in sorted(glob.glob("database_stein_watkins-20110713/swdb/*.bz2")):
    grp = "prim" if fn.split("/")[-1].startswith("p.") else "allg"
    lines = bz2.open(fn, "rt").read().splitlines()
    pairs, hdr = [], None
    for ln in lines:
        if ln[:1].isdigit(): hdr = ln
        elif ln.startswith("[") and hdr is not None:
            pairs.append((hdr, ln)); hdr = None          # erste Kurve je Klasse
    for h, c in pairs:
        N, bad, r = h.split()[:3]
        a = json.loads(c.split()[0])
        if int(N) < 10**6: continue                       # nur große Führer
        rr = min(int(r), 4)
        pool.setdefault((grp, rr), []).append((int(N), int(r), a))
primes = [int(p) for p in pari.primes([2, PMAX])]
out = open("sw.jsonl", "w"); stats = {}
for key in sorted(pool):
    lst = pool[key]
    pick = lst if len(lst) <= PER else random.sample(lst, PER)
    for N, r, a in pick:
        E = pari.ellinit(a)
        Np = int(pari.ellglobalred(E)[0])
        if Np != N:
            print("Führer-Abweichung", N, Np, a); continue
        w = int(pari.ellrootno(E))
        if w != (-1) ** r:
            print("Parität passt nicht zum SW-Rang", N, r, w); continue
        an = pari.ellan(E, PMAX)
        out.write(json.dumps(dict(label=f"SW{N}_{a}", N=N, rank=r, group=key[0], rootno=w, a=a,
                                  primes=primes, ap=[int(an[p - 1]) for p in primes])) + "\n")
        stats[key] = stats.get(key, 0) + 1
out.close()
for k in sorted(pool): print(k, "Pool", len(pool[k]), "gewählt", stats.get(k, 0))
print("gesamt", sum(stats.values()))
