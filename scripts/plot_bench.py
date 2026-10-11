import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
SURF, TXT, TXT2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"; BLUE, ORANGE = "#2a78d6", "#eb6834"
plt.rcParams.update({"font.size": 10.5, "axes.edgecolor": GRID, "axes.labelcolor": TXT2, "xtick.color": TXT2,
                     "ytick.color": TXT2, "figure.facecolor": SURF, "axes.facecolor": SURF})
cert = json.load(open("certified.json"))
need = lambda c, k: c[k]["primes"] if c[k] and c[k]["proven"] else np.inf
S = np.array([need(c, "std_r1") for c in cert]); R = np.array([need(c, "radix_r1") for c in cert])
fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.6, 7.2), dpi=170)
b = np.unique(np.concatenate([S[np.isfinite(S)], R[np.isfinite(R)]]))
a1.step(b, [100*np.mean(S <= x) for x in b], where="post", color=ORANGE, lw=2, label="Standard (Fejér/Bober)")
a1.step(b, [100*np.mean(R <= x) for x in b], where="post", color=BLUE, lw=2, label="Radix-Kern")
a1.set_xscale("log"); a1.set_xlabel("Primzahl-Budget"); a1.set_ylabel("Ränge beweisfest bestimmt (%)")
a1.set_title(f"{len(cert)} Kurven (Cremona, N bis 4·10⁵): Anteil bewiesener Ränge", loc="left", color=TXT, fontsize=11)
a1.legend(frameon=False, loc="lower right"); a1.set_ylim(0, 100)
labs, vals = [], []
for r in range(4):
    m = np.array([c["rank"] == r for c in cert]) & np.isfinite(S) & np.isfinite(R)
    labs.append(f"Rang {r}"); vals.append(100 * (1 - R[m].sum() / S[m].sum()))
bars = a2.bar(labs, vals, color=BLUE, width=0.55)
for bb, v in zip(bars, vals):
    a2.text(bb.get_x() + bb.get_width()/2, v + 0.4, f"{v:.1f} %", ha="center", fontsize=9.5, color=TXT)
a2.set_ylabel("eingesparte Primzahlen"); a2.set_ylim(0, max(vals) * 1.25)
a2.set_title("Ersparnis des Radix-Kerns je Rang (Kriterium Schranke < r+1)", loc="left", color=TXT, fontsize=11)
for a in (a1, a2):
    a.spines[["top", "right"]].set_visible(False); a.grid(axis="y", color=GRID, lw=0.6)
fig.tight_layout(); fig.savefig("benchmark.png"); print("ok")
