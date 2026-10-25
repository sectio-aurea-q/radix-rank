"""
Close the open curves of the ICARM high-rank board: certified Radix bounds with more primes.

The board run (scripts/batch_icarm.py, cap 1e10) left curves where the Radix kernel had not yet reached the
target. This script carries each of them up a ladder of prime bounds until the target is reached:

  1. prep_curve.py   conductor, root number, target (as in the board run)
  2. ap              a_p for all p < P, extended stage by stage (blocks of AP_BLOCK = 1e9 are kept, so a larger P
                     only computes the new range)
  3. radix_stream    first window below the target for Radix and Bober, and both bounds at the largest window
  4. cert_fast.py    rigorous certification (Rust interval prime sum + Arb) of the Radix hit, and of Bober's hit or
                     its bound at the same data, so every table entry compares both kernels on identical primes
Curves are processed in order of the remaining gap (smallest first). Results: <work>/close/<id>.json, resumable
(finished curves are skipped), and <work>/close/summary.md + summary.csv.

Usage: python3 scripts/close_board.py [--ids 616,620] [--ladder 2e10,3e10,5e10,8e10,1.2e11] [--keep-ap]
                                      [--board work/batch/summary.csv] [--work work]
Disk: a_p data take about 6 bytes per prime (25 GB at P = 1e11); they are deleted after each curve unless --keep-ap.
"""
import argparse, csv, json, math, os, re, shutil, subprocess, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, 'target', 'release')
SCRIPTS = os.path.join(ROOT, 'scripts')
BLOCK = 10**9


def run(cmd, **kw):
    return subprocess.run(cmd, check=True, text=True, capture_output=True, **kw)


def mid(arb_txt):
    m = re.match(r'\[\s*([-+0-9.eE]+)', str(arb_txt))
    return float(m.group(1)) if m else float('nan')


def prime_count(x):
    return x / (math.log(x) - 1.1) if x > 10 else 4


def open_curves(board):
    rows = [r for r in csv.DictReader(open(board)) if r['radix_proven'] != 'True']
    for r in rows:
        r['gap'] = mid(r['radix_arb']) - float(r['target'])
    return sorted(rows, key=lambda r: r['gap'])


def certify(W, key):
    o = run([sys.executable, os.path.join(SCRIPTS, 'cert_fast.py'), os.path.join(W, 'curve.json'),
             os.path.join(W, 'rap', 'c'), os.path.join(W, 'stream.json'), key])
    return json.loads(o.stdout.strip().splitlines()[-1])


def one_curve(c, ladder, work, keep_ap, log):
    cid = c['id']
    W = os.path.join(work, 'close', f'icarm{cid}')
    rap = os.path.join(W, 'rap')
    os.makedirs(rap, exist_ok=True)
    json.dump(c, open(os.path.join(W, 'icarm.json'), 'w'))
    t0 = time.time()
    log(run([sys.executable, os.path.join(SCRIPTS, 'prep_curve.py'), os.path.join(W, 'icarm.json'), W]).stdout.strip())
    cv = json.load(open(os.path.join(W, 'curve.json')))
    res = dict(id=cid, rank=cv['rank'], logN=cv['logN'], root_number=cv['root_number'], target=cv['target'], stages=[])
    stream = None
    for P in ladder:
        need = 6 * prime_count(P) * 1.05
        free = shutil.disk_usage(W).free
        have = sum(os.path.getsize(os.path.join(rap, f)) for f in os.listdir(rap))
        if need - have > free:
            log(f"  #{cid}: not enough disk for p < {P:.0e} (need {(need - have) / 1e9:.0f} GB more, "
                f"free {free / 1e9:.0f} GB) -- stopping this curve")
            break
        if P % BLOCK:  # not block-aligned (small test ladders): start this stage from scratch
            shutil.rmtree(rap, ignore_errors=True)
            os.makedirs(rap)
        ts = time.time()
        run([os.path.join(BIN, 'ap'), cv['c4'], cv['c6'], os.path.join(W, 'small_ap.txt'), '0', repr(float(P)),
             os.path.join(rap, 'c')], env=dict(os.environ, AP_BLOCK=str(BLOCK)))
        ta = time.time() - ts
        out = run([os.path.join(BIN, 'radix_stream'), os.path.join(W, 'curve.json'), os.path.join(rap, 'c'),
                   str(cv['target']), '10', '0.01', '2'])
        stream = json.loads(out.stdout)
        json.dump(stream, open(os.path.join(W, 'stream.json'), 'w'))
        st = dict(pmax=P, ap_secs=round(ta, 1), radix_hit=bool(stream['radix_hit']), std_hit=bool(stream['std_hit']),
                  radix_end=stream['radix_end']['bound'], std_end=stream['std_end']['bound'])
        res['stages'].append(st)
        rtxt = 'HIT' if st['radix_hit'] else f"{st['radix_end']:.4f}"
        btxt = 'hit' if st['std_hit'] else f"{st['std_end']:.4f}"
        log(f"  #{cid} p < {P:.1e}: Radix {rtxt}  Bober {btxt}  (a_p {ta:.0f} s)")
        if stream['radix_hit']:
            break
    if stream is None:
        return None
    res['pmax_used'] = stream['pmax']
    for k in ('radix', 'std'):
        h = stream[f'{k}_hit']
        key = f'{k}_hit' if h else f'{k}_end'
        src = h or stream[f'{k}_end']
        cert = certify(W, key)
        res[k] = dict(reached=bool(h), L=src['l'], primes=src['primes'], bound=src['bound'], arb=cert['arb'],
                      proven=cert['proven'], cert_secs=cert['secs'], certified=key)
        log(f"  #{cid} {k:5s} {key}: Arb {cert['arb']} < {cv['target']}: {cert['proven']}")
    res['secs'] = round(time.time() - t0, 1)
    if not keep_ap:
        shutil.rmtree(rap, ignore_errors=True)
    return res


def summary(outdir):
    rows = []
    for f in sorted(os.listdir(outdir)):
        if f.endswith('.json'):
            r = json.load(open(os.path.join(outdir, f)))
            if 'radix' in r:
                rows.append(r)
    rows.sort(key=lambda r: (r['rank'], r['logN']))
    with open(os.path.join(outdir, 'summary.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['id', 'rank', 'logN', 'target', 'pmax', 'radix_reached', 'radix_primes', 'radix_arb', 'radix_proven',
                    'bober_reached', 'bober_primes', 'bober_arb', 'secs'])
        for r in rows:
            w.writerow([r['id'], r['rank'], f"{r['logN']:.2f}", r['target'], f"{r['pmax_used']:.2e}",
                        r['radix']['reached'], r['radix']['primes'], r['radix']['arb'], r['radix']['proven'],
                        r['std']['reached'], r['std']['primes'], r['std']['arb'], r['secs']])
    lines = ['# Closing the open curves of the high-rank board', '',
             'Certified (Arb) analytic-rank bounds under GRH, target rank+2 (parity). Both kernels on the same primes.', '',
             '| ICARM | rank | ln N | target | primes (Radix) | Radix bound (Arb) | Bober on the same primes |',
             '|---|---|---|---|---|---|---|']
    for r in rows:
        rp = f"{r['radix']['primes']:,}" if r['radix']['reached'] else f"> {r['radix']['primes']:,}"
        bb = 'reached' if r['std']['reached'] else f"{r['std']['arb']} (not reached)"
        lines.append(f"| #{r['id']} | {r['rank']} | {r['logN']:.1f} | < {r['target']} | {rp} | {r['radix']['arb']} | {bb} |")
    n_ok = sum(1 for r in rows if r['radix'].get('proven'))
    lines += ['', f"{len(rows)} curves processed; Radix proven on {n_ok}."]
    open(os.path.join(outdir, 'summary.md'), 'w').write('\n'.join(lines) + '\n')
    return len(rows), n_ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ids', default='')
    ap.add_argument('--ladder', default='2e10,3e10,5e10,8e10,1.2e11')
    ap.add_argument('--keep-ap', action='store_true')
    ap.add_argument('--work', default=os.path.join(ROOT, 'work'))
    ap.add_argument('--board', default='')
    ap.add_argument('--db', default='')
    a = ap.parse_args()
    ladder = [float(x) for x in a.ladder.split(',')]
    for P in ladder:
        assert P % BLOCK == 0 or P < BLOCK, f"ladder values must be multiples of {BLOCK:.0e} (got {P:.3e})"
    subprocess.run(['cargo', 'build', '--release', '--quiet', '--manifest-path', os.path.join(ROOT, 'Cargo.toml')], check=True)
    dbp = a.db or os.path.join(a.work, 'database.json')
    if not os.path.exists(dbp):
        os.makedirs(a.work, exist_ok=True)
        subprocess.run(['curl', '-sSfL', 'https://elliptic-rank.icarm.cloud/database.json', '-o', dbp], check=True)
    db = {str(c['id']): c for c in json.load(open(dbp))['curves']}
    if a.ids:
        order = [x.strip() for x in a.ids.split(',')]
    else:
        board = a.board or os.path.join(ROOT, 'highrank', 'board', 'summary.csv')
        order = [r['id'] for r in open_curves(board)]
    outdir = os.path.join(a.work, 'close')
    os.makedirs(outdir, exist_ok=True)
    logf = open(os.path.join(outdir, 'close.log'), 'a')

    def log(s):
        print(s, flush=True)
        logf.write(s + '\n')
        logf.flush()

    log(f"== {len(order)} curves, ladder {', '.join(f'{p:.1e}' for p in ladder)}")
    for i, cid in enumerate(order):
        f = os.path.join(outdir, f"{cid}.json")
        if os.path.exists(f) and json.load(open(f)).get('radix', {}).get('proven'):
            continue
        log(f"[{i + 1}/{len(order)}] #{cid} rank >= {db[cid]['rank_lower_bound']}")
        try:
            res = one_curve(db[cid], ladder, a.work, a.keep_ap, log)
        except subprocess.CalledProcessError as e:
            log(f"  #{cid} FAILED: {' '.join(map(str, e.cmd[:2]))}: {(e.stderr or '').strip()[-300:]}")
            continue
        if res:
            json.dump(res, open(f, 'w'), indent=1)
            n, ok = summary(outdir)
            log(f"  done in {res['secs']:.0f} s; {ok} of {n} closed so far")
    summary(outdir)
    log(f"== results: {os.path.join(outdir, 'summary.md')}")


if __name__ == '__main__':
    main()
