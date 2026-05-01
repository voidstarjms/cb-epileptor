#!/bin/bash

CONDOR_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJ_DIR="$(cd "$CONDOR_DIR/.." && pwd)"
SWEEP_DIR="$PROJ_DIR/src/sweep"

RUN_NUM=1
while [ -d "$SWEEP_DIR/figures/${RUN_NUM}_sweep_debug" ]; do
    RUN_NUM=$((RUN_NUM + 1))
done

LOG_DIR="$CONDOR_DIR/logs/${RUN_NUM}_synchrony"
mkdir -p "$LOG_DIR"
echo "$RUN_NUM" > "$SWEEP_DIR/current_run.txt"

echo "Run number: $RUN_NUM"
echo "Created $LOG_DIR"
echo "Wrote $SWEEP_DIR/current_run.txt"

cat > "$CONDOR_DIR/condor.sub" << EOF
Universe   = vanilla
Executable = $CONDOR_DIR/wrapper.sh
Arguments  = run_single_sim.py --params \$(params_file)

Initialdir = $SWEEP_DIR

Output = $LOG_DIR/out_\$(Process).log
Error  = $LOG_DIR/err_\$(Process).log
Log    = $LOG_DIR/condor.log

Request_Cpus   = 1
Request_Memory = 8GB

+CSCI_GrpDesktop = True

Queue params_file from $SWEEP_DIR/params_list.txt
EOF
