#!/usr/bin/env bash
eval "$(conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate yasim_c_elegans_as_depth_analysis
set -vueo pipefail

SHDIR="$(readlink -f "$(dirname "${0}")")"
pip install -r "${SHDIR}"/spladder_requirements.txt
pip install --no-deps spladder==3.0.2
