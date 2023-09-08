# shellcheck shell=bash
# BaseRC for hot update without rebuilding.
# Change your working directory to here before sourcing this file.

PYTHONPATH="$(pwd):$(pwd)/src:$(pwd)/deps/labw_utils/src:${PYTHONPATH:-}"
export PYTHONPATH

PATH="$(pwd)/singularity/pbsim_builder:${PATH:-}"
export PATH
