"""Figure: the 16 hard curves of the ICARM high-rank board, Radix vs Bober on the same primes,
and the convergence of both bounds for #535. Reads highrank/close/summary_full.csv and highrank/close/535.json; writes paper/board16.pdf."""
import csv, json, re, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "serif", "font.serif": ["cmr10"], "mathtext.fontset": "cm", "axes.formatter.use_mathtext": True, "axes.unicode_minus": False,
                     "font.size": 9, "pdf.fonttype": 42, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6})
mid = lambda s: float(re.match(r"\[([-\d.]+)", s).group(1))
rows = list(csv.DictReader(open(sys.argv[1] if len(sys.argv) > 1 else "highrank/close/summary_full.csv")))
st = json.load(open(sys.argv[2] if len(sys.argv) > 2 else "highrank/close/535.json"))
rows.sort(key=lambda r: (int(r["target"]), float(r["logN"])))
INK, MUTED, GRID = "#111111", "#6b6b6b", "#d9d9d9"
fig, (a, b) = plt.subplots(1, 2, figsize=(6.3, 3.5), gridspec_kw={"width_ratios": [1.35, 1]})
for y, r in enumerate(rows):
    t = int(r["target"]); rd, bd = mid(r["radix_arb"]) - t, mid(r["bober_arb"]) - t
    done = r["radix_proven"] == "True"
    a.plot([rd, bd], [y, y], color=GRID, lw=1.2, zorder=1)
    a.plot(bd, y, "o", ms=5, mfc="white", mec=INK, mew=0.9, zorder=3)
    a.plot(rd, y, "o" if done else "o", ms=5, mfc=INK if done else "#9a9a9a", mec=INK, mew=0.9, zorder=4,
           **({} if done else {"ls": "none"}))
    if not done: a.text(bd + 0.02, y, "open at $p<5\\cdot10^{11}$", ha="left", va="center", fontsize=8, color=MUTED, style="italic")
a.axvline(0, color=INK, lw=0.8)
a.set_yticks(range(len(rows)), [f"#{r['id']}  (r = {r['rank']})" for r in rows], fontsize=8)
a.set_xlim(-0.17, 0.62)
a.set_xlabel("bound $-$ target  (target $= r+2$)")
a.text(-0.008, len(rows) - 0.3, "certified", ha="right", va="top", fontsize=8, color=MUTED); a.text(0.008, len(rows) - 0.3, "not reached", ha="left", va="top", fontsize=8, color=MUTED)
a.set_ylim(len(rows) + 0.4, -0.7)
for s in ("top", "right"): a.spines[s].set_visible(False)
a.set_title("(a) The 16 hardest board curves", fontsize=9, loc="left")
P = [s["pmax"] for s in st["stages"]]
R = [s["radix_end"] for s in st["stages"]]; B = [s["std_end"] for s in st["stages"]]
b.plot(P, B, "-o", color=INK, ms=4, mfc="white", mew=0.9, lw=1, label="fixed kernel of Bober [2]")
b.plot(P, R, "-o", color=INK, ms=4, mfc=INK, lw=1, label="radix Riemanniana (this paper)")
b.set_ylim(29.94, None); b.axhline(st["target"], color=INK, lw=0.8, ls="--"); b.text(P[0], st["target"] + 0.008, "target 30", va="bottom", fontsize=8, color=MUTED)
b.set_xscale("log"); import matplotlib.ticker as mt; b.xaxis.set_minor_formatter(mt.NullFormatter()); b.set_xticks([2e11, 4e11, 8e11], ["$2\\cdot10^{11}$", "$4\\cdot10^{11}$", "$8\\cdot10^{11}$"]); b.set_xlabel("prime bound $P$ (stage $p<P$)"); b.set_ylabel("bound (floating point; last stage certified)")
b.annotate("radix $<30$ at $8\\cdot10^{11}$", xy=(P[-1], R[-1]), xytext=(P[1], 29.965), fontsize=8,
           arrowprops=dict(arrowstyle="-", lw=0.6, color=MUTED))
for s in ("top", "right"): b.spines[s].set_visible(False)
b.set_title(f"(b) #{st['id']}, rank {st['rank']}, $\\log N = {st['logN']:.1f}$", fontsize=9, loc="left")
h, l = b.get_legend_handles_labels()
fig.legend(h[::-1], l[::-1], loc="lower center", ncol=2, frameon=False, fontsize=8, bbox_to_anchor=(0.5, 0.0))
fig.tight_layout(w_pad=1.5, rect=(0, 0.06, 1, 1))
fig.savefig("paper/board16.pdf")
