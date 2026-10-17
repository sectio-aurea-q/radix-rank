"""
Radix vs Bober on the ICARM high-rank board: certified analytic-rank upper bounds for many curves, one command.

For every selected curve (default: all with rank lower bound 24..31):
  1. prep_curve.py   conductor from local data, root number, target (rank+2 with parity, else rank+3)
  2. ap              a_p for all p < P, P = 1e6, 1e7, ... up to the cap, until Radix AND Bober reach the target
                     (the stage is recomputed from scratch each time; total cost <= 1.11x the final stage)
  3. radix_stream    first window below the target for Radix and for Bober (bisection, same data)
  4. cert_rap.py     Arb certification of every hit; if a kernel never reaches the target below the cap,
                     its bound at the largest window is certified instead and reported as "not reached"
Results: <work>/batch/<id>.json per curve (resumable: finished curves are skipped), and
<work>/batch/summary.md + summary.csv over all finished curves.

Usage: python3 scripts/batch_icarm.py [--min-rank 24] [--max-rank 31] [--cap 1e10] [--ids 3762,3835]
                                      [--no-cert] [--keep-ap] [--db database.json] [--redo-unreached]
  --redo-unreached: recompute curves where a kernel did not reach the target and the stored cap is below --cap
  (typical use: a fast pass with --cap 1e9 over everything, then --cap 1e10 --redo-unreached overnight)
"""
import argparse, csv, glob, json, math, os, shutil, subprocess, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, 'target', 'release')
SCRIPTS = os.path.join(ROOT, 'scripts')


def run(cmd, **kw):
    return subprocess.run(cmd, check=True, text=True, capture_output=True, **kw)


def stages(cap):
    p = 1e6
    while p < cap * 1.0001:
        yield p
        p *= 10
    if p / 10 < cap * 0.9999:
        yield cap


def one_curve(c, work, cap, certify, keep_ap, log):
    cid = c['id']
    W = os.path.join(work, f'icarm{cid}')
    os.makedirs(W, exist_ok=True)
    json.dump(c, open(os.path.join(W, 'icarm.json'), 'w'))
    t0 = time.time()
    r = run([sys.executable, os.path.join(SCRIPTS, 'prep_curve.py'), os.path.join(W, 'icarm.json'), W])
    log(r.stdout.strip())
    cv = json.load(open(os.path.join(W, 'curve.json')))
    target = cv['target']
    res = dict(id=cid, rank=cv['rank'], logN=cv['logN'], root_number=cv['root_number'], target=target,
               parity=(cv['root_number'] == (-1) ** cv['rank']), stages=[])
    rap = os.path.join(W, 'rap')
    stream = None
    for P in stages(cap):
        shutil.rmtree(rap, ignore_errors=True)
        os.makedirs(rap)
        ts = time.time()
        env = dict(os.environ, AP_BLOCK=str(max(1e6, P / 28)))
        run([os.path.join(BIN, 'ap'), cv['c4'], cv['c6'], os.path.join(W, 'small_ap.txt'), '0', repr(P),
             os.path.join(rap, 'c')], env=env)
        ta = time.time() - ts
        out = run([os.path.join(BIN, 'radix_stream'), os.path.join(W, 'curve.json'), os.path.join(rap, 'c'),
                   str(target), '10', '0.01', '2'])
        stream = json.loads(out.stdout)
        st = dict(pmax=P, ap_secs=round(ta, 1), radix_hit=bool(stream['radix_hit']), std_hit=bool(stream['std_hit']),
                  radix_end=stream['radix_end']['bound'], std_end=stream['std_end']['bound'])
        res['stages'].append(st)
        rtxt = 'hit' if st['radix_hit'] else f"{st['radix_end']:.4f}"
        btxt = 'hit' if st['std_hit'] else f"{st['std_end']:.4f}"
        log(f"  #{cid} p<{P:.0e}: Radix {rtxt}  Bober {btxt}  (a_p {ta:.0f}s)")
        if stream['radix_hit'] and stream['std_hit']:
            break
    json.dump(stream, open(os.path.join(W, 'stream.json'), 'w'))
    res['pmax_used'] = stream['pmax']
    for k in ('radix', 'std'):
        h = stream[f'{k}_hit']
        res[k] = dict(reached=bool(h), L=(h or stream[f'{k}_end'])['l'], primes=(h or stream[f'{k}_end'])['primes'],
                      bound=(h or stream[f'{k}_end'])['bound'])
        if certify:
            key = f'{k}_hit' if h else f'{k}_end'
            o = run([sys.executable, os.path.join(SCRIPTS, 'cert_rap.py'), os.path.join(W, 'curve.json'),
                     os.path.join(rap, 'c'), os.path.join(W, 'stream.json'), key])
            cert = json.loads(o.stdout.strip().splitlines()[-1])
            res[k]['arb'] = cert['arb']
            res[k]['proven'] = cert['proven']
            res[k]['cert_secs'] = cert['secs']
            log(f"  #{cid} {k:5s} {key}: Arb {cert['arb']} < {target}: {cert['proven']}")
    if res['radix']['reached'] and res['std']['reached']:
        res['prime_ratio'] = res['std']['primes'] / res['radix']['primes']
    res['secs'] = round(time.time() - t0, 1)
    if not keep_ap:
        shutil.rmtree(rap, ignore_errors=True)
    return res


def summary(outdir):
    rows = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(outdir, '*.json')))]
    rows = [r for r in rows if 'radix' in r]
    rows.sort(key=lambda r: (r['rank'], r['logN']))
    with open(os.path.join(outdir, 'summary.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['id', 'rank', 'logN', 'root_number', 'target', 'radix_reached', 'radix_primes', 'radix_arb',
                    'radix_proven', 'bober_reached', 'bober_primes', 'bober_arb', 'bober_proven', 'prime_ratio', 'pmax'])
        for r in rows:
            w.writerow([r['id'], r['rank'], f"{r['logN']:.2f}", r['root_number'], r['target'],
                        r['radix']['reached'], r['radix']['primes'], r['radix'].get('arb', ''), r['radix'].get('proven', ''),
                        r['std']['reached'], r['std']['primes'], r['std'].get('arb', ''), r['std'].get('proven', ''),
                        f"{r['prime_ratio']:.2f}" if 'prime_ratio' in r else '', f"{r['pmax_used']:.0e}"])
    lines = ['# Radix vs Bober on the ICARM high-rank board', '',
             'Primes needed for a certified (Arb) analytic-rank bound below the target (rank+2 when the root number',
             'gives parity). "> n" = target not reached with all primes below the cap; then the certified bound at',
             'the largest window is given. Under GRH.', '',
             '| ICARM | rank | ln N | target | Radix primes | Bober primes | ratio | Radix bound (Arb) |',
             '|---|---|---|---|---|---|---|---|']
    for r in rows:
        rp = f"{r['radix']['primes']:,}" if r['radix']['reached'] else f"> {r['radix']['primes']:,}"
        bp = f"{r['std']['primes']:,}" if r['std']['reached'] else f"> {r['std']['primes']:,}"
        ratio = f"{r['prime_ratio']:.1f}" if 'prime_ratio' in r else ('≥ ' + f"{r['std']['primes'] / r['radix']['primes']:.1f}"
                                                                       if r['radix']['reached'] else '—')
        lines.append(f"| #{r['id']} | {r['rank']} | {r['logN']:.1f} | < {r['target']:g} | {rp} | {bp} | {ratio} | "
                     f"{r['radix'].get('arb', round(r['radix']['bound'], 6))} |")
    done = [r for r in rows if 'prime_ratio' in r]
    if done:
        rs = sorted(r['prime_ratio'] for r in done)
        lines += ['', f"{len(rows)} curves; both kernels reached the target on {len(done)}. "
                  f"Prime ratio Bober/Radix: median {rs[len(rs)//2]:.1f}, min {rs[0]:.1f}, max {rs[-1]:.1f}. "
                  f"Radix proven on {sum(1 for r in rows if r['radix'].get('proven'))} curves."]
    open(os.path.join(outdir, 'summary.md'), 'w').write('\n'.join(lines) + '\n')
    return len(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--min-rank', type=int, default=24)
    ap.add_argument('--max-rank', type=int, default=31)
    ap.add_argument('--cap', type=float, default=1e10)
    ap.add_argument('--ids', default='')
    ap.add_argument('--no-cert', action='store_true')
    ap.add_argument('--keep-ap', action='store_true')
    ap.add_argument('--db', default='')
    ap.add_argument('--redo-unreached', action='store_true')
    ap.add_argument('--work', default=os.path.join(ROOT, 'work'))
    a = ap.parse_args()
    subprocess.run(['cargo', 'build', '--release', '--quiet', '--manifest-path', os.path.join(ROOT, 'Cargo.toml')], check=True)
    dbp = a.db or os.path.join(a.work, 'database.json')
    if not os.path.exists(dbp):
        os.makedirs(a.work, exist_ok=True)
        subprocess.run(['curl', '-sSfL', 'https://elliptic-rank.icarm.cloud/database.json', '-o', dbp], check=True)
    db = json.load(open(dbp))['curves']
    if a.ids:
        want = {int(x) for x in a.ids.split(',')}
        sel = [c for c in db if c['id'] in want]
    else:
        sel = [c for c in db if a.min_rank <= int(c['rank_lower_bound']) <= a.max_rank and c.get('conductor')]
    sel.sort(key=lambda c: (int(c['rank_lower_bound']), len(c['conductor'] or '')))
    outdir = os.path.join(a.work, 'batch')
    os.makedirs(outdir, exist_ok=True)
    logf = open(os.path.join(outdir, 'batch.log'), 'a')

    def log(s):
        print(s, flush=True)
        logf.write(s + '\n')
        logf.flush()

    log(f"== {len(sel)} curves, cap {a.cap:.0e}, certification {'off' if a.no_cert else 'on'}")
    for i, c in enumerate(sel):
        f = os.path.join(outdir, f"{c['id']}.json")
        if os.path.exists(f):
            old = json.load(open(f))
            unreached = not (old['radix']['reached'] and old['std']['reached'])
            if not (a.redo_unreached and unreached and old['pmax_used'] < a.cap * 0.999):
                continue
        log(f"[{i + 1}/{len(sel)}] #{c['id']} rank >= {c['rank_lower_bound']}")
        try:
            res = one_curve(c, a.work, a.cap, not a.no_cert, a.keep_ap, log)
        except subprocess.CalledProcessError as e:
            log(f"  #{c['id']} FAILED: {' '.join(e.cmd[:2])}: {(e.stderr or '').strip()[-300:]}")
            json.dump(dict(id=c['id'], error=(e.stderr or '')[-2000:]), open(f + '.err', 'w'))
            continue
        json.dump(res, open(f, 'w'), indent=1)
        n = summary(outdir)
        log(f"  done in {res['secs']:.0f}s; summary over {n} curves updated")
    summary(outdir)
    log(f"== results: {os.path.join(outdir, 'summary.md')}")


if __name__ == '__main__':
    main()
