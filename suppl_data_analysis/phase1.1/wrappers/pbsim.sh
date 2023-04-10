#!/usr/bin/env bash
set -e

eval "$(conda 'shell.bash' 'hook' 2>/dev/null)"
conda activate yasim-pbsim1
exec pbsim "${@}"
