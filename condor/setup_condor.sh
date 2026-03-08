#!/bin/bash
# Run this once before each sweep submission, from anywhere:
#   bash condor/setup_condor.sh
#
# It will:
#   1. Find the next run number (matches figures/{run_num}_sweep_debug convention)
#   2. Create condor/logs/{run_num}_synchrony/
#   3. Write src/current_run.txt so aggregate.py uses the same run number
#   4. Generate condor/condor.sub with the correct absolute paths for this machine

CONDOR_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJ_DIR="$(cd "$CONDOR_DIR/.." && pwd)"
FIGURES_DIR="$PROJ_DIR/src/figures"

# Find next run number — same logic as aggregate.py
RUN_NUM=1
while [ -d "$FIGURES_DIR/${RUN_NUM}_sweep_debug" ]; do
    RUN_NUM=$((RUN_NUM + 1))
done

LOG_DIR="$CONDOR_DIR/logs/${RUN_NUM}_synchrony"
mkdir -p "$LOG_DIR"
echo "$RUN_NUM" > "$PROJ_DIR/src/current_run.txt"

echo "Run number: $RUN_NUM"
echo "Created $LOG_DIR"
echo "Wrote src/current_run.txt"

cat > "$CONDOR_DIR/condor.sub" << EOF
Universe   = vanilla
Executable = $CONDOR_DIR/wrapper.sh
Arguments  = run_single_sim.py --ce \$(ce) --x0 \$(x0)

Initialdir = $PROJ_DIR/src

Output = $LOG_DIR/out_\$(Process).log
Error  = $LOG_DIR/err_\$(Process).log
Log    = $LOG_DIR/condor.log

Request_Cpus   = 1
Request_Memory = 2GB

Queue ce, x0 from $PROJ_DIR/src/params_list.txt
EOF

echo "Generated condor/condor.sub for run $RUN_NUM"
echo ""
echo "Next steps:"
echo "  1. cd src && python generate_params.py"
echo "  2. python run_single_sim.py --ce 0.25 --x0 -3.5   # test one job"
echo "  3. condor_submit condor/condor.sub"
echo "  4. condor_q"
echo "  5. cd src && python aggregate.py"
