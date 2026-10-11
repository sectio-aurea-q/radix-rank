#!/usr/bin/env bash
# radix-rank — reproduce the results of the paper.
#   ./reproduce.sh quick     # ~1 min: 7 reference curves (ranks 0-6), build + optimise + certify
#   ./reproduce.sh cremona   # ~10 min: full Cremona benchmark (1642 curves) from scratch
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
MODE="${1:-quick}"
cargo build --release --quiet --manifest-path "$ROOT/Cargo.toml"
BIN="$ROOT/target/release/radix-rank"
export PYTHONPATH="$ROOT/scripts"
mkdir -p "$ROOT/work" && cd "$ROOT/work"
case "$MODE" in
  quick)
    "$BIN" "$ROOT/data/reference_curves.jsonl" reference_out.jsonl 0.3 4.5 0.05 10
    python3 "$ROOT/scripts/certify.py" "$ROOT/data/reference_curves.jsonl" reference_out.jsonl ;;
  cremona)
    mkdir -p data
    for f in allcurves.00000-09999 allcurves.100000-109999 allcurves.400000-409999; do
      [ -f "data/$f" ] || curl -sSfL -o "data/$f" "https://raw.githubusercontent.com/JohnCremona/ecdata/master/allcurves/$f"
    done
    python3 "$ROOT/scripts/export.py" 8200 150
    "$BIN" curves.jsonl results.jsonl 0.3 4.5 0.05 10
    python3 "$ROOT/scripts/bench_certify.py"
    python3 "$ROOT/scripts/analyze.py" ;;
  *) echo "usage: $0 [quick|cremona]"; exit 1 ;;
esac
