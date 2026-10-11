import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
SURF, TXT, TXT2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"; BLUE, AQUA = "#2a78d6", "#1baf7a"
plt.rcParams.update({"font.size": 10.5, "axes.edgecolor": GRID, "axes.labelcolor": TXT2, "xtick.color": TXT2,
                     "ytick.color": TXT2, "figure.facecolor": SURF, "axes.facecolor": SURF})
s = json.load(open("sw_summary.json"))
fig, axes = plt.subplots(2, 1, figsize=(6.6, 7.6), dpi=170)
for ax, crit, title in [(axes[0], "r+2", "mit Wurzelzahl (Schranke < r+2) – der Praxisfall"),
                        (axes[1], "r+1", "ohne Wurzelzahl (Schranke < r+1)")]:
    tab = s[crit]["table"]; x = np.arange(5); w = 0.38
    for j, (g, col) in enumerate([("N 4·10⁶–10⁸", BLUE), ("N prim 3,6–7,5·10⁹", AQUA)]):
        vals = [next(t["saving"] for t in tab if t["group"] == g and t["rank"] == r) for r in range(5)]
        bars = ax.bar(x + (j - 0.5) * (w + 0.02), vals, w, color=col, label=f"Führer {g[2:]}")
        for bb, v in zip(bars, vals):
            ax.text(bb.get_x() + bb.get_width() / 2, v + 0.6, f"{v:.0f}", ha="center", fontsize=8.5, color=TXT)
    ax.set_xticks(x); ax.set_xticklabels(["Rang 0", "Rang 1", "Rang 2", "Rang 3", "Rang ≥4"])
    ax.set_ylabel("eingesparte Primzahlen (%)")
    ax.set_title(f"Stein-Watkins, 2500 Kurven, {title}\nRadix vs. Standard: gesamt −{s[crit]['saving']:.0f} %, "
                 f"nie schlechter; bewiesen {s[crit]['radix']} statt {s[crit]['std']}", loc="left", color=TXT, fontsize=10)
    ax.set_ylim(0, 45); ax.spines[["top", "right"]].set_visible(False); ax.grid(axis="y", color=GRID, lw=0.6)
axes[0].legend(frameon=False, loc="upper right", fontsize=9)
fig.tight_layout(); fig.savefig("stein_watkins.png"); print("ok")
