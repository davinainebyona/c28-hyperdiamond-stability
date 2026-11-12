#!/usr/bin/env bash
# Helpers for the GitHub Actions workflow.
#   ci.sh setup            install Quantum ESPRESSO, Python packages and the pseudopotential
#   ci.sh publish MSG      commit dft/results to main (retrying if another job pushed first)
#   ci.sh logs NAME        save input files and the tail of each output to the dft-logs branch
set -u
cd "$(dirname "$0")"
REMOTE="https://x-access-token:${GITHUB_TOKEN:-}@github.com/${GITHUB_REPOSITORY:-}.git"
ident() { git -C "$1" config user.name "github-actions[bot]"; git -C "$1" config user.email "41898282+github-actions[bot]@users.noreply.github.com"; }

case "$1" in
  setup)
    set -e
    sudo apt-get update -qq
    sudo apt-get install -y -qq quantum-espresso
    pip install -q ase==3.29.0 phonopy==4.10.0 spglib==2.8.0 numpy
    mkdir -p pseudo results
    curl -fsSL -o pseudo/C.pbe-n-kjpaw_psl.1.0.0.UPF https://pseudopotentials.quantum-espresso.org/upf_files/C.pbe-n-kjpaw_psl.1.0.0.UPF
    test "$(wc -c < pseudo/C.pbe-n-kjpaw_psl.1.0.0.UPF)" -gt 100000
    which pw.x; nproc
    ;;
  publish)
    D=$(mktemp -d)
    git clone -q --depth 1 --branch main "$REMOTE" "$D" || exit 1
    ident "$D"
    for i in 1 2 3 4 5 6; do
      mkdir -p "$D/dft/results"; cp -f results/*.json results/*.cif results/*.npy "$D/dft/results/" 2>/dev/null
      git -C "$D" add -A dft/results
      git -C "$D" diff --cached --quiet && { echo "nothing to publish"; exit 0; }
      git -C "$D" commit -q -m "$2"
      git -C "$D" push -q origin main && { echo "published"; exit 0; }
      sleep $((RANDOM % 15 + 3)); git -C "$D" fetch -q origin main; git -C "$D" reset -q --hard origin/main
    done
    exit 1
    ;;
  logs)
    D=$(mktemp -d)
    git clone -q --depth 1 --branch dft-logs "$REMOTE" "$D" 2>/dev/null || { git init -q -b dft-logs "$D"; git -C "$D" remote add origin "$REMOTE"; }
    ident "$D"
    for i in 1 2 3 4 5 6; do
      mkdir -p "$D/$2"
      for f in *.out; do [ -f "$f" ] && tail -n 250 "$f" > "$D/$2/$f.tail.txt"; done
      for f in *.out; do [ -f "$f" ] && grep -E "^!|P=|new unit-cell volume|Total force|convergence has been achieved|PWSCF.*WALL|Sym. Ops|number of k points" "$f" > "$D/$2/$f.summary.txt"; done
      cp -f *.in "$D/$2/" 2>/dev/null
      date -u > "$D/$2/updated.txt"
      git -C "$D" add -A; git -C "$D" commit -q -m "logs: $2" || exit 0
      git -C "$D" push -q origin dft-logs && exit 0
      sleep $((RANDOM % 15 + 3)); git -C "$D" pull -q --rebase origin dft-logs || true
    done
    exit 0
    ;;
esac
