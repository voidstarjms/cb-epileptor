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
mkdir -p "$SWEEP_DIR/data/results"   # where transferred-back metrics pkls land
echo "$RUN_NUM" > "$SWEEP_DIR/current_run.txt"

echo "Run number: $RUN_NUM"
echo "Created $LOG_DIR"
echo "Wrote $SWEEP_DIR/current_run.txt"

cat > "$CONDOR_DIR/condor.sub" << EOF
Universe   = vanilla
Executable = $CONDOR_DIR/wrapper.sh
Arguments  = run_single_sim.py --params \$(params_file)

Initialdir = $PROJ_DIR/src/sweep

should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
transfer_input_files    = $PROJ_DIR/src, $CONDOR_DIR/cn_venv.tar.gz
transfer_output_files   = metrics_\$Fn(params_file).pkl
transfer_output_remaps  = "metrics_\$Fn(params_file).pkl = data/results/metrics_\$Fn(params_file).pkl"

Output = $LOG_DIR/out_\$(Process).log
Error  = $LOG_DIR/err_\$(Process).log
Log    = $LOG_DIR/condor.log

Request_Cpus   = 1
Request_Memory = 16GB

+CSCI_GrpDesktop = True

Queue params_file from $SWEEP_DIR/params_list.txt
EOF
