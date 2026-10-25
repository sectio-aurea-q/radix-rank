"""
Append certified Radix upper bounds to the ICARM leaderboard commentary of each curve.

Safety first. ICARM keeps every commentary edit in the curve's history, but these are public entries on other
people's curves, so this script:
  * runs as a DRY RUN by default: it writes every note and every new commentary to <out>/preview.md and posts nothing;
  * with --apply asks once for confirmation, then for each curve
      - re-reads the live entry, skips it if the rank lower bound differs from the certified rank, if a radix-rank
        note is already present, or (unless --include-existing) if the commentary already states an explicit-formula
        bound under GRH,
      - saves the old commentary to <out>/backup/<id>.txt,
      - posts old commentary + blank line + note (the existing text is never changed),
      - reads the entry back and checks that the old text is an exact prefix and the note is present;
    it stops at the first failed check;
  * --restore <id,...> puts the saved commentary back.

Data (only certified results, nothing is recomputed or guessed):
  highrank/board/curves/<id>.json    board run, one file per curve (collect them with --collect)
  highrank/close/<id>.json           curves closed by scripts/close_board.py
  highrank/torsion2/<id>.json        curves with rational 2-torsion
Each file must contain rank, root_number, target and radix = {reached, proven, L, primes, arb}.
Every note links to its own certificate file on GitHub; the link is checked before anything is posted.

Usage:
  python3 scripts/icarm_notes.py --collect             copy work/batch and work/close certificates into highrank/
  python3 scripts/icarm_notes.py                       dry run -> work/icarm_notes/preview.md
  python3 scripts/icarm_notes.py --apply [--ids 630,628] [--limit 5]
  python3 scripts/icarm_notes.py --fix-4239 [--apply]  "r_an = 32 subject to GRH" -> "r_an <= 32 subject to GRH"
  python3 scripts/icarm_notes.py --restore 630
Token: environment variable ICARM_TOKEN, otherwise asked for (not echoed).
"""
import argparse, datetime, getpass, glob, json, math, os, re, shutil, sys, time, urllib.parse, urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SITE = os.environ.get('ICARM_SITE', 'https://elliptic-rank.icarm.cloud')
REPO_RAW = 'https://raw.githubusercontent.com/sectio-aurea-q/radix-rank/main/'
REPO_BLOB = 'https://github.com/sectio-aurea-q/radix-rank/blob/main/'
SOURCES = [('highrank/close', 'close'), ('highrank/torsion2', 'torsion2'), ('highrank/board/curves', 'board')]
UA = 'radix-rank-icarm-notes/1.0 (+https://github.com/sectio-aurea-q/radix-rank)'
EXISTING = re.compile(r'radix-rank|explicit[- ]formula|1112\.1503|Bober|rank exactly \d|exact rank \d|exact rank of \d|unconditional', re.I)


# ------------------------------------------------------------------ http
def http(method, url, token=None, data=None, tries=4):
    body = urllib.parse.urlencode(data).encode() if data is not None else None
    for k in range(tries):
        req = urllib.request.Request(url, data=body, method=method, headers={'User-Agent': UA})
        if token:
            req.add_header('authorization', f'Bearer {token}')
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.status, r.read().decode()
        except urllib.error.HTTPError as e:
            if e.code == 429 and k < tries - 1:
                time.sleep(10 * (k + 1)); continue
            return e.code, e.read().decode(errors='replace')
        except urllib.error.URLError:
            if k < tries - 1:
                time.sleep(5 * (k + 1)); continue
            raise


def live(cid):
    st, txt = http('GET', f'{SITE}/curve/{cid}.json')
    if st != 200:
        raise RuntimeError(f'curve {cid}: GET returned {st}')
    return json.loads(txt)


# ------------------------------------------------------------------ certificates
def arb_upper(arb):
    m = re.fullmatch(r'\[\s*([-+0-9.eE]+)\s*\+/-\s*([0-9.eE+-]+)\s*\]', arb.strip())
    if not m:
        raise ValueError(f'cannot read Arb enclosure {arb!r}')
    return float(m.group(1)) + float(m.group(2))


def load_certs():
    certs = {}
    for d, kind in SOURCES:
        for f in sorted(glob.glob(os.path.join(ROOT, d, '[0-9]*.json'))):
            r = json.load(open(f))
            cid = str(r.get('id', os.path.basename(f)[:-5]))
            x = r.get('radix') or {}
            if not (x.get('reached') and x.get('proven')):
                continue
            if cid in certs:          # first source wins: close > torsion2 > board
                continue
            certs[cid] = dict(id=cid, kind=kind, rank=int(r['rank']), root_number=int(r['root_number']),
                              target=float(r['target']), L=float(x['L']), primes=int(x['primes']), arb=x['arb'],
                              date=r.get('certified_on') or file_date(f), iso=r.get('isogenous_to'),
                              path=os.path.relpath(f, ROOT).replace(os.sep, '/'))
    return certs


def file_date(f):
    return datetime.date.fromtimestamp(os.path.getmtime(f)).isoformat()


def check_cert(c):
    """Internal consistency of one certificate; returns a list of problems (empty = fine)."""
    p = []
    r = c['rank']
    if c['root_number'] != (-1) ** r:
        p.append(f"root number {c['root_number']} does not give parity (-1)^{r}")
    if c['target'] != r + 2:
        p.append(f"target {c['target']} is not rank+2")
    try:
        if not arb_upper(c['arb']) < c['target']:
            p.append(f"Arb enclosure {c['arb']} is not below the target {c['target']}")
    except ValueError as e:
        p.append(str(e))
    if not (0 < c['L'] < 20 and c['primes'] > 0):
        p.append(f"implausible window L={c['L']} / primes={c['primes']}")
    return p


def note(c):
    r = c['rank']
    x = f"{math.exp(2 * c['L']):.3g}".replace('e+0', 'e').replace('e+', 'e')
    sign = '+1' if c['root_number'] == 1 else '-1'
    iso = (f" The 2-isogenous curve#{c['iso']} has the same L-function, hence the same bound."
           if c.get('iso') else '')
    return (f"[edit MEG]: Upper bound ({c['date']}, Sabrina Hackfort (MEG)): r_an <= {r} subject to GRH. "
            f"Explicit-formula bound with a per-curve optimised test function (radix-rank), L = {c['L']:.2f}, "
            f"all {c['primes']:,} primes p < e^(2L) = {x}, certified in interval arithmetic (Arb): "
            f"{c['arb']} < {r + 2:g}; root number {sign}, so r_an has the parity of {r}. "
            f"With the {r} independent points the rank is exactly {r} subject to GRH and BSD. "
            f"{iso.strip()}{' ' if iso else ''}Certificate: {REPO_BLOB}{c['path']}")


def link_ok(path):
    st, _ = http('GET', REPO_RAW + path)
    return st == 200


# ------------------------------------------------------------------ actions
def collect():
    """Copy certificates from the work folder into the repository: board run (work/batch) and closing run
    (work/close; certificates written by an older close_board.py that did not copy them itself)."""
    for src_dir, dst_rel in (('batch', 'highrank/board/curves'), ('close', 'highrank/close')):
        src = sorted(glob.glob(os.path.join(ROOT, 'work', src_dir, '[0-9]*.json')))
        dst = os.path.join(ROOT, dst_rel)
        n = 0
        for f in src:
            r = json.load(open(f))
            if 'radix' not in r:
                continue
            os.makedirs(dst, exist_ok=True)
            r.setdefault('certified_on', file_date(f))
            json.dump(r, open(os.path.join(dst, os.path.basename(f)), 'w'), indent=1)
            n += 1
        print(f"{n} certificates from work/{src_dir} copied to {dst_rel}/")
    print("commit and push them before --apply")


def plan(certs, ids, include_existing, out, check_links):
    rows, skipped = [], []
    for cid in ids:
        c = certs[cid]
        probs = check_cert(c)
        if probs:
            skipped.append((cid, 'certificate: ' + '; '.join(probs))); continue
        cur = live(cid)
        com = cur.get('commentary') or ''
        if int(cur['rank_lower_bound']) != c['rank']:
            skipped.append((cid, f"leaderboard rank is now {cur['rank_lower_bound']}, certificate is for {c['rank']}")); continue
        if 'radix-rank' in com:
            skipped.append((cid, 'radix-rank note already present')); continue
        if EXISTING.search(com) and not include_existing:
            skipped.append((cid, 'commentary already states an upper bound or the exact rank (use --include-existing)')); continue
        if c.get('iso') and not include_existing:
            pc = live(c['iso']).get('commentary') or ''
            if EXISTING.search(pc) and 'radix-rank' not in pc:
                skipped.append((cid, f"isogenous curve#{c['iso']} already states an upper bound or the exact rank")); continue
        if check_links and not link_ok(c['path']):
            skipped.append((cid, f"certificate file not on GitHub yet: {c['path']} (commit and push first)")); continue
        n = note(c)
        new = (com + '\r\n\r\n' + n) if com else n
        rows.append(dict(id=cid, cert=c, old=com, note=n, new=new))
        time.sleep(0.2)
    lines = ['# ICARM notes (preview, nothing posted)', '',
             f'{len(rows)} notes to post, {len(skipped)} curves skipped.', '']
    for r in rows:
        lines += [f"## curve#{r['id']} (rank {r['cert']['rank']}, {r['cert']['kind']})", '', r['note'], '']
    if skipped:
        lines += ['## skipped', ''] + [f'- curve#{i}: {why}' for i, why in skipped]
    return rows, skipped, '\n'.join(lines) + '\n'


def post(cid, content, token):
    st, txt = http('POST', f'{SITE}/curve/{cid}/commentary', token=token, data={'content': content})
    if st >= 300:
        raise RuntimeError(f'curve {cid}: POST returned {st}: {txt[:200]}')


def apply(rows, token, out, sleep):
    bdir = os.path.join(out, 'backup'); os.makedirs(bdir, exist_ok=True)
    log = open(os.path.join(out, 'applied.jsonl'), 'a')
    for k, r in enumerate(rows):
        cid = r['id']
        cur = live(cid).get('commentary') or ''
        if cur != r['old']:
            print(f"curve#{cid}: commentary changed since the dry run -- skipped, run again"); continue
        open(os.path.join(bdir, f'{cid}.txt'), 'w', newline='').write(cur)
        post(cid, r['new'], token)
        time.sleep(1)
        after = live(cid).get('commentary') or ''
        ok_prefix = after.startswith(cur) if cur else True
        ok_note = r['note'] in after
        log.write(json.dumps(dict(id=cid, time=datetime.datetime.now().isoformat(timespec='seconds'),
                                  prefix_kept=ok_prefix, note_present=ok_note, note=r['note'])) + '\n'); log.flush()
        print(f"[{k + 1}/{len(rows)}] curve#{cid}: original kept {ok_prefix}, note present {ok_note}")
        if not (ok_prefix and ok_note):
            print("STOP: check failed. Old text is in", os.path.join(bdir, f'{cid}.txt'),
                  f"-- restore with: python3 scripts/icarm_notes.py --restore {cid}")
            sys.exit(2)
        time.sleep(sleep)


def restore(ids, token, out):
    for cid in ids:
        f = os.path.join(out, 'backup', f'{cid}.txt')
        old = open(f, newline='').read()
        post(cid, old, token)
        ok = (live(cid).get('commentary') or '') == old
        print(f"curve#{cid}: restored {ok}")


def fix_4239(token, do_apply, out):
    cur = live(4239).get('commentary') or ''
    wrong = 'r_an = 32 subject to GRH'
    right = 'r_an <= 32 subject to GRH'
    if cur.count(wrong) != 1:
        print('4239: expected phrase not found exactly once -- nothing changed'); return
    new = cur.replace(wrong, right)
    print('4239: ' + wrong + '  ->  ' + right)
    if not do_apply:
        print('(dry run; add --apply to post)'); return
    os.makedirs(os.path.join(out, 'backup'), exist_ok=True)
    open(os.path.join(out, 'backup', '4239.txt'), 'w', newline='').write(cur)
    post(4239, new, token)
    after = live(4239).get('commentary') or ''
    print('4239: fixed', after == new)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--apply', action='store_true')
    ap.add_argument('--ids', default='')
    ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--include-existing', action='store_true')
    ap.add_argument('--collect', action='store_true')
    ap.add_argument('--restore', default='')
    ap.add_argument('--fix-4239', action='store_true')
    ap.add_argument('--no-link-check', action='store_true')
    ap.add_argument('--sleep', type=float, default=3.0)
    ap.add_argument('--out', default=os.path.join(ROOT, 'work', 'icarm_notes'))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if a.collect:
        collect(); return
    need_token = a.apply or a.restore
    token = os.environ.get('ICARM_TOKEN') or (getpass.getpass('ICARM token (not shown): ') if need_token else None)
    if a.restore:
        restore([x.strip() for x in a.restore.split(',')], token, a.out); return
    if a.fix_4239:
        fix_4239(token if a.apply else None, a.apply, a.out); return
    certs = load_certs()
    ids = [x.strip() for x in a.ids.split(',')] if a.ids else sorted(certs, key=lambda i: (-certs[i]['rank'], int(i)))
    ids = [i for i in ids if i in certs and i != '4239']
    if a.limit:
        ids = ids[:a.limit]
    print(f"{len(certs)} certificates found; checking {len(ids)} curves against the live leaderboard ...")
    rows, skipped, md = plan(certs, ids, a.include_existing, a.out, not a.no_link_check)
    open(os.path.join(a.out, 'preview.md'), 'w').write(md)
    print(f"{len(rows)} notes ready, {len(skipped)} skipped -> {os.path.join(a.out, 'preview.md')}")
    if not a.apply:
        print('dry run: nothing was posted. Read the preview, then run again with --apply.'); return
    if not rows:
        return
    ans = input(f"Post {len(rows)} notes to ICARM under your account? Type JA to continue: ")
    if ans.strip() != 'JA':
        print('aborted, nothing posted'); return
    apply(rows, token, a.out, a.sleep)
    print('done. Log:', os.path.join(a.out, 'applied.jsonl'))


if __name__ == '__main__':
    main()
