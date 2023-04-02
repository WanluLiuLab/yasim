#shell=bash
# BaseRC for hot update without rebuilding.
# Change your working directory to here before sourcing this file.

export PYTHONPATH="$(pwd):$(pwd)/src:$(pwd)/deps/labw_utils/src:$(pwd)/deps/pbcore:${PYTHONPATH:-}"
