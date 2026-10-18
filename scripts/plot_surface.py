"""Abbildung: die Radix Riemanniana als Flaeche. Quelle: paper/radix_surface.json
(erzeugt von scripts/radix_surface.py). Aufruf aus dem Wurzelverzeichnis."""
import json, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
# Type 42 (TrueType) instead of matplotlib's default Type 3: Type 3 fonts are not
# scalable, which breaks pdfTeX font expansion (microtype) and is rejected by several
# journals.  The serif family matches the body font of the paper.
matplotlib.rcParams.update({
    'pdf.fonttype': 42, 'ps.fonttype': 42,
    'font.family': 'serif', 'font.serif': ['DejaVu Serif'],
    'mathtext.fontset': 'dejavuserif',
})
import matplotlib.pyplot as plt
from matplotlib import cm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
d = json.load(open(os.path.join(ROOT, 'paper', 'radix_surface.json')))
L = np.array(d['L']); t = np.array(d['t']); H = np.array(d['h']); zs = d['zeros']
TMAX = 9.0                      # jenseits davon ist die Flaeche numerisch null
m = t <= TMAX; t = t[m]; H = H[:, m]
zs = [z for z in zs if z <= TMAX]

TT, LL = np.meshgrid(t, L)
fig = plt.figure(figsize=(9.6, 3.6))

ax = fig.add_subplot(1, 2, 1, projection='3d')
ax.plot_surface(TT, LL, H, cmap=cm.bone_r, rstride=1, cstride=1,
                linewidth=0, antialiased=True, alpha=.95)
for z in zs:
    ax.plot([z, z], [L[0], L[-1]], [0, 0], color='#a8322a', lw=1.0, zorder=10)
ax.set_xlabel('height $t$', labelpad=-4); ax.set_ylabel('window $L$', labelpad=-4)
ax.set_zlabel('$\\hat f(t)^2$', labelpad=-6)
ax.tick_params(labelsize=7, pad=-2); ax.view_init(elev=24, azim=-128)
ax.set_box_aspect((2.0, 1.5, 1.1))
ax.set_title('the radix as a surface', fontsize=9, pad=-2)

ax2 = fig.add_subplot(1, 2, 2)
for Lv, col, lw in ((0.6, '#c0c0c0', 1.0), (1.5, '#8a8a8a', 1.2), (3.0, '#141414', 1.6)):
    i = int(np.argmin(abs(L - Lv)))
    ax2.plot(t, H[i], color=col, lw=lw, label=f'$L={L[i]:g}$')
for z in zs:
    ax2.axvline(z, color='#a8322a', lw=.8, alpha=.8, zorder=0)
ax2.set_xlim(0, TMAX); ax2.set_ylim(0, None)
ax2.set_xlabel('height $t$'); ax2.set_ylabel('$\\hat f(t)^2$')
ax2.set_title('cross-sections; red: the zeros of $L(E,s)$', fontsize=9)
ax2.legend(fontsize=8, frameon=False); ax2.grid(alpha=.2, lw=.5)
ax2.tick_params(labelsize=8)
for sp in ('top', 'right'): ax2.spines[sp].set_visible(False)

plt.tight_layout(w_pad=1.4)
out = os.path.join(ROOT, 'paper', 'surface.pdf')
plt.savefig(out, bbox_inches='tight'); print(out)
