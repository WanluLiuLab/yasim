#!/usr/bin/env bash
echo "Deprecated due to lack of pre-trained models"
exit 1
set -e
SHDIR="$(dirname "$(readlink -f "${0}")")"
if which simulatION &>> /dev/null; then
    exec simulatION "${@}"
fi
if ! which conda &>> /dev/null; then
    echo "conda not found!"
    exit 127
fi

if ! conda env list | grep ^nanopore_simulation &>> /dev/null; then
    conda create -y -n nanopore_simulation python=3
fi
eval "$(conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate nanopore_simulation
if which simulatION &>> /dev/null; then
    exec "${CONDA_PREFIX}"/bin/python "$(which simulatION)" "${@}"
fi
"${CONDA_PREFIX}"/bin/pip install "${SHDIR}"/../../../../dist/simulatION-0.3-py3-none-any.whl
exec "${CONDA_PREFIX}"/bin/python "$(which simulatION)" "${@}"
