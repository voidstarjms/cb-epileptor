#!/bin/bash

# Ensure correct number of arguments
if [ $# != 1 ]; then
    echo "Usage ./setup_condor.sh [run_name]"
    exit 1
fi

CONDOR_DIR="$(realpath .)"
PROJ_DIR="$(realpath ..)"
SWEEP_DIR="$PROJ_DIR/src/sweep"
PARAM_DIR="$SWEEP_DIR/params/$1/"

# Make sure specified run exists
if [ ! -d $PARAM_DIR ]; then
    echo "Parameter directory ${$PARAM_DIR} does not exist"
    exit 2
fi

LOG_DIR="$CONDOR_DIR/logs/${$1}_synchrony"
mkdir -p "$LOG_DIR"
mkdir -p "$SWEEP_DIR/data/results"   # where transferred metrics pkls land

echo "Run name: $1"
echo "Created $LOG_DIR"
echo "Wrote $SWEEP_DIR/current_run.txt"

cat > "$CONDOR_DIR/condor.sub" << EOF
Universe   = vanilla
Executable = $CONDOR_DIR/wrapper.sh
Arguments  = run_single_sim.py --params \$(params_file)

Initialdir = $SWEEP_DIR

should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
transfer_input_files    = $PROJ_DIR/src, $CONDOR_DIR/cn_venv.tar.gz
transfer_output_files   = metrics_\$Fn(params_file).pkl
transfer_output_remaps  = "metrics_\$Fn(params_file).pkl = data/results/metrics_\$Fn(params_file).pkl"

Output = $LOG_DIR/out_\$(Process).log
Error  = $LOG_DIR/err_\$(Process).log
Log    = $LOG_DIR/condor.log

Request_Cpus   = 1
Request_Memory = 32GB

+CSCI_GrpDesktop = True

Queue params_file from $PARAM_DIR/params_list.txt
EOF
