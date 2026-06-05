#!/bin/sh
# Build a relocatable venv tarball for HTCondor file transfer.
#
# Run this ONCE on a pool machine (so the packed interpreter matches the execute
# computers), and re-run whenever dependencies change. Output: condor/cn_venv.tar.gz,
# which condor.sub ships to each job's scratch sandbox.
set -e

CONDOR_DIR="$(cd "$(dirname "$0")" && pwd)"

# Point this at your actual venv. The repo has conflicting references — README says
# SR-CB/.cn_venv, the old wrapper.sh resolved to CompNeuro/.cn_venv, and the local copy
# is at the repo root — so make it overridable (CN_VENV=... ) and verify before packing.
VENV="${CN_VENV:-$CONDOR_DIR/../.cn_venv}"

if [ ! -x "$VENV/bin/python" ]; then
    echo "No venv found at $VENV (set CN_VENV=/path/to/.cn_venv to override)" >&2
    exit 1
fi

pip show venv-pack >/dev/null 2>&1 || pip install venv-pack
venv-pack -p "$VENV" -o "$CONDOR_DIR/cn_venv.tar.gz" --force

echo "Wrote $CONDOR_DIR/cn_venv.tar.gz from $VENV"
