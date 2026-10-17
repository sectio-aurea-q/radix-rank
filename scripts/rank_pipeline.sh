#!/usr/bin/env bash
# Certify an analytic-rank bound (under GRH) for an ICARM leaderboard curve, streamed and parallel:
#   1. prep_curve.py   conductor + root number check, target from parity, a_p table for p < 1000 and bad primes
#   2. ap              a_p for all primes < PMAX (Rust BSGS, all cores, resumable, compact .rap blocks)
#   3. radix_stream    first window with bound < target (Radix and Bober), and both bounds at the largest window
#   4. cert_fast.py    rigorous certification of every hit and of both endpoints (Rust intervals + Arb)
#
# usage:  scripts/rank_pipeline.sh <icarm_id> <pmax> [K=10]
#   e.g.  scripts/rank_pipeline.sh 4239 2e10
# Re-running continues where it stopped (finished a_p blocks are kept); a larger PMAX extends the data.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
ID="$1"; PMAX="$2"; K="${3:-10}"
W="$ROOT/work/icarm$ID"
mkdir -p "$W/rap"
cargo build --release --quiet --manifest-path "$ROOT/Cargo.toml"
BIN="$ROOT/target/release"
[ -s "$W/icarm.json" ] || curl -sSfL "https://elliptic-rank.icarm.cloud/curve/$ID.json" -o "$W/icarm.json"
python3 "$ROOT/scripts/prep_curve.py" "$W/icarm.json" "$W"
read -r C4 C6 TARGET < <(python3 -c "import json; c = json.load(open('$W/curve.json')); print(c['c4'], c['c6'], c['target'])")

echo "== a_p for p < $PMAX"
t0=$(date +%s)
"$BIN/ap" "$C4" "$C6" "$W/small_ap.txt" 0 "$PMAX" "$W/rap/c" 2> >(grep -v "^points" >&2)
echo "   done in $(( $(date +%s) - t0 )) s"

echo "== windows (target: bound < $TARGET)"
"$BIN/radix_stream" "$W/curve.json" "$W/rap/c" "$TARGET" "$K" > "$W/stream.json"
python3 - "$W/stream.json" <<'EOF'
import json, sys
d = json.load(open(sys.argv[1]))
for k in ('radix_hit', 'std_hit', 'radix_end', 'std_end'):
    h = d[k]
    print(f"   {k:10s}", "not reached" if h is None else f"L = {h['l']:.2f}  primes {h['primes']:>12,}  bound {h['bound']:.6f}")
EOF

echo "== Arb certification"
for key in radix_hit std_hit radix_end std_end; do
  if python3 -c "import json,sys; sys.exit(0 if json.load(open('$W/stream.json'))['$key'] else 1)"; then
    python3 "$ROOT/scripts/cert_fast.py" "$W/curve.json" "$W/rap/c" "$W/stream.json" "$key" | tee -a "$W/certificates.jsonl"
  fi
done
echo "== results: $W/stream.json, $W/certificates.jsonl"
