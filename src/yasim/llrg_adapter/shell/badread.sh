#!/usr/bin/env bash
set -e
if which badread &>> /dev/null; then
    exec badread "${@}"
fi
if ! which conda &>> /dev/null; then
    echo "conda not found!"
    exit 127
fi

if ! conda env list | grep ^badread &>> /dev/null; then
    conda create -y -n badread -c bioconda badread==0.2.0 python-edlib
fi

eval "$(conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate badread
exec badread "${@}"
