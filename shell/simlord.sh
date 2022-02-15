#!/usr/bin/env bash
set -e
if ! which conda &>> /dev/null; then
    echo "conda not found!"
    exit 127
fi

if ! conda env list | grep ^simlord &>> /dev/null; then
    conda create -y -n simlord -c bioconda simlord
fi

eval "$(conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate simlord
simlord "${@}"
