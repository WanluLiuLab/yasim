#!/usr/bin/env bash
set -e
if which nanosim &>> /dev/null; then
    exec nanosim "${@}"
fi
if ! which conda &>> /dev/null; then
    echo "conda not found!"
    exit 127
fi

if ! conda env list | grep ^nanosim &>> /dev/null; then
    conda create -y -n nanosim -c bioconda nanosim==3.0.0
fi

eval "$(conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate nanosim
exec nanosim "${@}"
