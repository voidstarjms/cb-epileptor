#!/bin/bash
# Write condor/condor.sub for one sweep run and create its output folders.
#
#   condor/setup_condor.sh <run_name>
#
# Expects src/sweep/params/<run_name>/ from generate_params.py. Each job's raw
# output comes back to output/sweep/<run_name>/jobs/<job_id>/output.pkl, next to
# a copy of the run's sweep_config.yaml, ready for build_table.py.
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
SWEEP_DIR="$PROJ_DIR/src/sweep"
PARAM_DIR="$SWEEP_DIR/params/$RUN"
JOB_LIST="$PARAM_DIR/job_list.txt"
OUT_DIR="$PROJ_DIR/output/sweep/$RUN"
LOG_DIR="$CONDOR_DIR/logs/$RUN"
REQUEST_MEMORY="${REQUEST_MEMORY:-32GB}"

if [ ! -f "$JOB_LIST" ]; then
    echo "No job list at $JOB_LIST; run generate_params.py --name $RUN first"
    exit 2
fi

# Refuse to mix a new submission into an earlier run's results.
if [ -d "$OUT_DIR/jobs" ] && [ -n "$(find "$OUT_DIR/jobs" -name output.pkl -print -quit)" ]; then
    echo "$OUT_DIR/jobs already has results; use a new run name"
    exit 3
fi

# Condor doesn't create directories for transfer_output_remaps, so make one per job.
while read -r job_path; do
    [ -z "$job_path" ] && continue
    mkdir -p "$OUT_DIR/jobs/$(basename "$job_path" .yaml)"
done < "$JOB_LIST"
cp "$PARAM_DIR/sweep_config.yaml" "$OUT_DIR/sweep_config.yaml"
mkdir -p "$LOG_DIR"

cat > "$CONDOR_DIR/condor.sub" << EOF
Universe   = vanilla
Executable = $CONDOR_DIR/wrapper.sh
Arguments  = run_single_sim.py --params \$(params_file)

Initialdir = $SWEEP_DIR

should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
transfer_input_files    = $PROJ_DIR/src, $CONDOR_DIR/cn_venv.tar.gz
transfer_output_files   = output_\$Fn(params_file).pkl
transfer_output_remaps  = "output_\$Fn(params_file).pkl = $OUT_DIR/jobs/\$Fn(params_file)/output.pkl"

Output = $LOG_DIR/out_\$(Process).log
Error  = $LOG_DIR/err_\$(Process).log
Log    = $LOG_DIR/condor.log

Request_Cpus   = 1
Request_Memory = $REQUEST_MEMORY

+CSCI_GrpDesktop = True

Queue params_file from $JOB_LIST
EOF

echo "Run: $RUN ($(grep -c . "$JOB_LIST") jobs)"
echo "Results will land in $OUT_DIR/jobs/"
echo "Logs in $LOG_DIR"
echo "Wrote $CONDOR_DIR/condor.sub; submit with: condor_submit $CONDOR_DIR/condor.sub"
