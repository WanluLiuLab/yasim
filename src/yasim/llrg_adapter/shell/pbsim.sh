#!/usr/bin/env bash
set -e
if which pbsim &>> /dev/null; then
    exec pbsim "${@}"
fi
if ! which conda &>> /dev/null; then
    echo "conda not found!"
    exit 127
fi

if ! conda env list | grep ^pbsim &>> /dev/null; then
    conda create -y -n pbsim -c bioconda pbsim==1.0.3
fi

eval "$(conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate pbsim
exec pbsim "${@}"
