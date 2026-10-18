"""
Generate paper/board.tex (Table "board") and paper/board_macros.tex (macros \\BoardAbstract, \\BoardBrief) from work/batch/summary.csv.

Usage: python3 scripts/make_board.py [work/batch/summary.csv] [paper/board.tex]
Also copies summary.csv/summary.md to highrank/board/ for the repository.
"""
import csv, os, shutil, statistics, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, 'work', 'batch', 'summary.csv')
dst = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, 'paper', 'board.tex')

rows = list(csv.DictReader(open(src)))
T = lambda x: x == 'True'
by = {}
for r in rows:
    by.setdefault(int(r['rank']), []).append(r)


def med(xs):
    return statistics.median(xs) if xs else None


lines = [r'\begin{table}[!ht]', r'\centering\small',
         r'\begin{tabular}{rrrrrl}', r'\toprule',
         r'rank & curves & optimal & fixed & only optimal & median ratio\\', r'\midrule']
tot = dict(n=0, rp=0, bp=0, only=0, both=0)
ratios_all = []
meds = []
for k in sorted(by):
    rs = by[k]
    rp = sum(T(r['radix_proven']) for r in rs)
    bp = sum(T(r['bober_proven']) for r in rs)
    only = sum(T(r['radix_proven']) and not T(r['bober_proven']) for r in rs)
    both = sum(T(r['radix_proven']) and T(r['bober_proven']) for r in rs)
    rat = [float(r['prime_ratio']) for r in rs if r['prime_ratio']]
    ratios_all += rat
    m = med(rat)
    meds.append((k, m))
    lines.append(f"{k} & {len(rs)} & {rp} & {bp} & {only} & {m:.1f}\\\\" if m else
                 f"{k} & {len(rs)} & {rp} & {bp} & {only} & ---\\\\")
    for a, b in (('n', len(rs)), ('rp', rp), ('bp', bp), ('only', only), ('both', both)):
        tot[a] += b
lines += [r'\midrule',
          f"all & {tot['n']} & {tot['rp']} & {tot['bp']} & {tot['only']} & {med(ratios_all):.1f}\\\\",
          r'\bottomrule', r'\end{tabular}',
          r'\caption{ICARM curves of rank $24$--$31$: number of curves on which each kernel certifies $\ran\le r$ under GRH '
          r"with $p<CAP$; ``only optimal'' = certified by the optimal kernel but not by the fixed one; ratio = primes "
          r'needed by the fixed kernel divided by primes needed by the optimal kernel, over curves where both succeed.}',
          r'\label{tab:board}', r'\end{table}']
pmax = max(float(r['pmax']) for r in rows)
exp = f"{pmax:.0e}".replace('e+', 'e').split('e')
cap = f"10^{{{int(exp[1])}}}" if exp[0] == '1' else f"{exp[0]}\\cdot10^{{{int(exp[1])}}}"
mac = [f"\\newcommand{{\\BoardAbstract}}{{On all {tot['n']} leaderboard curves of rank $24$--$31$, with primes below ${cap}$, "
       f"the optimised kernel certifies $\\ran\\le r$ on {tot['rp']} curves and the fixed kernel on {tot['bp']}; "
       f"where both succeed the fixed kernel needs a median of {med(ratios_all):.1f} times as many primes.}}",
       f"\\newcommand{{\\BoardBrief}}{{On the {tot['n']} leaderboard curves of rank $24$--$31$ ($p<{cap}$) the optimised kernel "
       f"certifies $\\ran\\le r$ on {tot['rp']} curves, the fixed kernel on {tot['bp']}; {tot['only']} only by the optimised kernel, none only by "
       f"the fixed one (Section~\\ref{{sec:board}}).}}"]
mid = lambda a: float(a.strip('[').split()[0])
neither = [r for r in rows if not T(r['radix_proven']) and not T(r['bober_proven'])]
lower = sum(mid(r['radix_arb']) < mid(r['bober_arb']) for r in neither)
mac.append(f"\\newcommand{{\\BoardDetail}}{{On {len(neither)} curves neither kernel reaches the target with $p<{cap}$; "
           f"on {'all' if lower == len(neither) else lower} of them the certified bound of the optimised kernel is the smaller one; their certified bounds at the "
           f"cap are listed in the repository.}}")
assert not any(T(r['bober_proven']) and not T(r['radix_proven']) for r in rows), 'Bober-only case: edit BoardBrief'
ms = [m for _, m in meds if m]
mac.append(f"\\newcommand{{\\BoardTrend}}{{The median ratio grows with the rank, though not monotonically, from {ms[0]:.1f} at rank "
           f"${meds[0][0]}$ to {ms[-1]:.1f} at rank ${meds[-1][0]}$ (Table~\\ref{{tab:board}}).}}")
open(dst, 'w').write('\n'.join(lines).replace('CAP', cap) + '\n')
open(dst.replace('board.tex', 'board_macros.tex'), 'w').write('\n'.join(mac) + '\n')
bd = os.path.join(ROOT, 'highrank', 'board')
os.makedirs(bd, exist_ok=True)
for f in ('summary.csv', 'summary.md'):
    p = os.path.join(os.path.dirname(src), f)
    if os.path.exists(p):
        shutil.copy(p, bd)
print(f"{dst}: {tot}")
