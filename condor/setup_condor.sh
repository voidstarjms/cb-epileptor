#!/bin/bash
# Write condor/condor.sub for one sweep run and create its runs/ folder.
#
#   condor/setup_condor.sh <run_name>
#
# Expects data/sweep/<run_name>/ from generate_params.py. Each job's raw output
# comes back to data/sweep/<run_name>/runs/<job_id>.pkl, ready for build_table.py.
# Request_Memory defaults to 32GB; override with REQUEST_MEMORY=8GB once a job's
# real usage is known (see the condor.log "Memory (MB)" usage line).
set -e

if [ $# != 1 ]; then
    echo "Usage: $0 <run_name>"
    exit 1
fi
RUN="$1"

CONDOR_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJ_DIR="$(cd "$CONDOR_DIR/.." && pwd)"
RUN_DIR="$PROJ_DIR/data/sweep/$RUN"
JOB_LIST="$RUN_DIR/job_list.txt"
RUNS_DIR="$RUN_DIR/runs"
LOG_DIR="$CONDOR_DIR/logs/$RUN"
REQUEST_MEMORY="${REQUEST_MEMORY:-32GB}"

if [ ! -f "$JOB_LIST" ]; then
    echo "No job list at $JOB_LIST; run generate_params.py --name $RUN first"
    exit 2
fi

# Refuse to mix a new submission into an earlier run's results.
if [ -d "$RUNS_DIR" ] && [ -n "$(ls -A "$RUNS_DIR")" ]; then
    echo "$RUNS_DIR already has results; use a new run name"
    exit 3
fi

# Condor doesn't create directories for transfer_output_remaps.
mkdir -p "$RUNS_DIR" "$LOG_DIR"

cat > "$CONDOR_DIR/condor.sub" << EOF
Universe   = vanilla
Executable = $CONDOR_DIR/wrapper.sh
# wrapper.sh runs from src/sweep; the job YAML is transferred to the sandbox root.
Arguments  = run_single_sim.py --params ../../\$Fnx(params_file)

# job_list.txt paths are relative to the project root.
Initialdir = $PROJ_DIR

should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
transfer_input_files    = $PROJ_DIR/src, $CONDOR_DIR/cn_venv.tar.gz, \$(params_file)
transfer_output_files   = output_\$Fn(params_file).pkl
transfer_output_remaps  = "output_\$Fn(params_file).pkl = $RUNS_DIR/\$Fn(params_file).pkl"

Output = $LOG_DIR/out_\$(Process).log
Error  = $LOG_DIR/err_\$(Process).log
Log    = $LOG_DIR/condor.log

Request_Cpus   = 1
Request_Memory = $REQUEST_MEMORY

+CSCI_GrpDesktop = True

Queue params_file from $JOB_LIST
EOF

echo "Run: $RUN ($(grep -c . "$JOB_LIST") jobs)"
echo "Results will land in $RUNS_DIR/"
echo "Logs in $LOG_DIR"
echo "Wrote $CONDOR_DIR/condor.sub; submit with: condor_submit $CONDOR_DIR/condor.sub"
