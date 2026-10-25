"""
Curve-by-curve comparison of radix-rank with rankbound (R. Jain, 2026) on the ICARM high-rank board.

The two methods bound different things under different hypotheses:
  radix-rank  analytic rank r_an(E), GRH for L(E,s); exact rank needs BSD in addition
  rankbound   2-Selmer rank, hence rank E(Q), GRH for finitely many Dedekind zeta functions; no BSD;
              only for curves without rational 2-torsion; cannot close a curve whose 2-Selmer rank exceeds its rank
The table counts, for every board curve, whether each method certifies the exact value.

Usage: python3 scripts/compare_rankbound.py [--board highrank/board/summary.csv] [--close work/close/summary.csv]
                                            [--rankbound certificates.tsv] [--db work/database.json] [--out work/compare]
Without --rankbound the database of rankbound is downloaded (commit pinned below, so the comparison is reproducible).
"""
import argparse, csv, json, os, subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RB_COMMIT = 'f324d1bfff14fad0bcdbb73c75a53f0729e187e3'
RB_URL = f'https://raw.githubusercontent.com/rjain2470/rankbound/{RB_COMMIT}/certificates.tsv'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--board', default=os.path.join(ROOT, 'highrank', 'board', 'summary.csv'))
    ap.add_argument('--close', default=os.path.join(ROOT, 'work', 'close', 'summary.csv'))
    ap.add_argument('--rankbound', default='')
    ap.add_argument('--db', default=os.path.join(ROOT, 'work', 'database.json'))
    ap.add_argument('--out', default=os.path.join(ROOT, 'work', 'compare'))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rbp = a.rankbound or os.path.join(a.out, 'rankbound_certificates.tsv')
    if not os.path.exists(rbp):
        subprocess.run(['curl', '-sSfL', RB_URL, '-o', rbp], check=True)

    radix = {}
    for r in csv.DictReader(open(a.board)):
        radix[r['id']] = dict(rank=int(r['rank']), logN=float(r['logN']), proven=r['radix_proven'] == 'True',
                              primes=int(r['radix_primes']), arb=r['radix_arb'], src='board, p < 1e10')
    if os.path.exists(a.close):
        for r in csv.DictReader(open(a.close)):
            if r['radix_proven'] == 'True':
                radix[r['id']].update(proven=True, primes=int(r['radix_primes']), arb=r['radix_arb'],
                                      src=f"close, p < {float(r['pmax']):.1e}")
    rb = {r['icarm_id']: r for r in csv.DictReader(open(rbp), delimiter='\t')}

    rows, cnt = [], {'both': 0, 'radix_only': 0, 'rb_only': 0, 'neither': 0, 'rb_missing': 0}
    for cid, x in sorted(radix.items(), key=lambda kv: (kv[1]['rank'], kv[1]['logN'])):
        y = rb.get(cid)
        rb_ok = bool(y) and y['rank = r?'] == 'yes'
        if not y:
            cnt['rb_missing'] += 1
        key = ('both' if x['proven'] and rb_ok else 'radix_only' if x['proven'] else 'rb_only' if rb_ok else 'neither')
        cnt[key] += 1
        rows.append(dict(id=cid, rank=x['rank'], logN=x['logN'], radix=x['proven'], radix_primes=x['primes'],
                         radix_arb=x['arb'], radix_src=x['src'],
                         rankbound=(y['rank upper bound under GRH'] if y and y['rank upper bound under GRH'] else
                                    (y['verified'] if y else 'not in database')),
                         rankbound_exact=rb_ok))
    with open(os.path.join(a.out, 'compare.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    n = len(rows)
    lines = ['# radix-rank and rankbound on the ICARM high-rank board', '',
             f'{n} curves (rank 24-31). Exact value certified (under the respective hypotheses):', '',
             '| | count |', '|---|---|',
             f"| both methods | {cnt['both']} |", f"| radix-rank only | {cnt['radix_only']} |",
             f"| rankbound only | {cnt['rb_only']} |", f"| neither | {cnt['neither']} |", '',
             f"radix-rank: {cnt['both'] + cnt['radix_only']} of {n}; rankbound: {cnt['both'] + cnt['rb_only']} of {n}.", '',
             '## Curves where the two differ', '',
             '| ICARM | rank | ln N | radix-rank | rankbound upper bound |', '|---|---|---|---|---|']
    for r in rows:
        if r['radix'] != r['rankbound_exact']:
            rr = f"{r['radix_arb']} ({r['radix_primes']:,} primes, {r['radix_src']})" if r['radix'] else 'open'
            lines.append(f"| #{r['id']} | {r['rank']} | {r['logN']:.1f} | {rr} | {r['rankbound']} |")

    if os.path.exists(a.db):
        db = json.load(open(a.db))['curves']
        t2 = []
        for c in db:
            tors = c.get('torsion') or []
            if isinstance(tors, str):
                tors = json.loads(tors.replace("'", '"')) if tors.strip() not in ('', 'None') else []
            if int(c['rank_lower_bound']) >= 20 and any(int(t) % 2 == 0 for t in tors):
                t2.append((int(c['rank_lower_bound']), c['id'], tors))
        t2.sort(reverse=True)
        lines += ['', '## Outside the scope of rankbound: rational 2-torsion', '',
                  f"{len(t2)} leaderboard curves of rank >= 20 have a rational point of order 2. "
                  'rankbound does not apply to them; radix-rank does (scripts/batch_icarm.py --ids ...).', '',
                  'ids: ' + ','.join(str(i) for _, i, _ in t2)]
    open(os.path.join(a.out, 'compare.md'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines[:11]))
    print(f"-> {os.path.join(a.out, 'compare.md')}, compare.csv")


if __name__ == '__main__':
    main()
