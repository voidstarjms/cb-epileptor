#!/bin/bash
set -e

for run in run1 run2 run3 run4; do
    echo "===== $run ====="
    python3 run.py -m t --no-cb --params parameters/params4.yaml --out-dir output/$run
    echo ""
done
