"""Figure 1: certified bounds for ICARM #4239 as a function of the window L.

Data: the twelve certified evaluations recorded in work/icarm4239/stream65.json
(fields "evals": [L, primes, fixed-kernel bound, optimal-kernel bound]).
Run from the repository root:  python3 scripts/plot_highrank.py
"""
import json, os, sys
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC  = os.path.join(ROOT, 'work', 'icarm4239', 'stream65.json')
FALLBACK = [[5.0,2466,48.090489237359705,47.830162695150705],
            [8.72,2290685,37.41959756954147,37.06619671231559],
            [10.58,76979786,35.499406932593686,35.12355891509995],
            [11.51,452518634,34.867086866304014,34.481985634569696],
            [11.97,1089752410,34.61284095718861,34.22211063888523],
            [12.2,1692168645,34.49044390692539,34.10113395327363],
            [12.32,2129234429,34.43103130127891,34.04135555401495],
            [12.38,2388539593,34.40283987865482,34.01244922622535],
            [12.39,2434733560,34.39821614717292,34.00768671417168],
            [12.4,2481819193,34.3936028146321,34.002936073013814],
            [12.41,2529822599,34.389020290187936,33.9982044781576],
            [12.44,2679479553,34.37520986459815,33.98403649247535]]
ev = json.load(open(SRC))['evals'] if os.path.exists(SRC) else FALLBACK
L = [e[0] for e in ev]; B = [e[2] for e in ev]; R = [e[3] for e in ev]

GREY, INK, RED = '#8a8a8a', '#141414', '#a8322a'
fig, ax = plt.subplots(1, 2, figsize=(9.4, 3.2))
for a in ax:
    a.axhline(34, color=RED, lw=1.2, ls=(0, (5, 3)), zorder=1)
    a.plot(L, B, 'o-', color=GREY, ms=3.2, lw=1.3, zorder=2, label='fixed kernel')
    a.plot(L, R, 'o-', color=INK,  ms=3.2, lw=1.7, zorder=3, label='optimal kernel')
    a.set_xlabel('window $L$'); a.grid(alpha=.22, lw=.5); a.tick_params(labelsize=8)
    for sp in ('top', 'right'):
        a.spines[sp].set_visible(False)
ax[0].set_xlim(4.6, 12.8); ax[0].set_ylim(33.5, 48.8)
ax[0].set_ylabel('bound')
ax[0].set_title('all twelve evaluated windows', fontsize=9)
ax[0].text(12.0, 35.4, 'target $34$', color=RED, fontsize=8, ha='right')
ax[0].legend(fontsize=8, frameon=False, loc='upper right', bbox_to_anchor=(1.0, 1.02))
ax[1].set_xlim(12.33, 12.47); ax[1].set_ylim(33.955, 34.45)
ax[1].set_title('the crossing', fontsize=9)
ax[1].plot([12.41], [33.9982044781576], '*', color=RED, ms=15, zorder=6, clip_on=False)
ax[1].annotate('$[33.99820448\\pm1.8\\cdot10^{-9}]<34$\n$L=12.41$,  $2\\,529\\,822\\,599$ primes',
               (12.411, 33.9982), textcoords='offset points', xytext=(10, 30), fontsize=7.6,
               ha='left', va='bottom',
               arrowprops=dict(arrowstyle='-', lw=.7, color=RED, shrinkA=0, shrinkB=6))
ax[1].text(12.468, 34.012, 'target $34$', color=RED, fontsize=8, ha='right', va='bottom')
plt.tight_layout(w_pad=2.2)
out = os.path.join(ROOT, 'paper', 'highrank.pdf')
plt.savefig(out, bbox_inches='tight')
print(out)
