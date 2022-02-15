#!/usr/bin/env bash
set -ue
if not which conda; then
    echo "conda not found!"
    exit 127
fi

if not conda env list | grep ^badread; then
    conda create -n badread -c bioconda badread
fi

conda activate badread
badread "${@}"
