#!/usr/bin/env bash
set -

set -C
if ! echo -e "sync\t${$}" >sync.lock 2>/dev/null; then
    set +C
    echo -e "$(timestamp)\tSYNC\tOCCUPIED" >>act.log
    errh "Repository being synced by $(cat sync.lock)"
fi
set +C



if which pbsim2 &>> /dev/null; then
    exec pbsim2 "${@}"
fi
if ! which conda &>> /dev/null; then
    echo "conda not found!"
    exit 127
fi

if ! conda env list | grep ^pbsim2 &>> /dev/null; then
    conda create -y -n pbsim2 -c bioconda pbsim2==2.0.1
fi

eval "$(conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate pbsim2
exec pbsim "${@}"
