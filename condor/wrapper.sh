#!/bin/sh
# Activates the project venv and runs python with any arguments passed in.
# Used as the HTCondor executable.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
. "$SCRIPT_DIR/../../.cn_venv/bin/activate"
exec python "$@"
