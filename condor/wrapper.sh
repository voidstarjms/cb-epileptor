#!/bin/sh
# HTCondor executable. Runs inside the per-job scratch sandbox with the transferred
# inputs present (the src/ tree and cn_venv.tar.gz). Unpacks the packed venv onto the
# execute machine's local disk, activates it, then runs the job from src/sweep.
set -e

mkdir -p venv
tar -xzf cn_venv.tar.gz -C venv
. venv/bin/activate
cd src/sweep
exec python "$@"
