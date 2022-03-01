#!/usr/bin/env bash
echo "Deprecated due to error number of trials when simulating ONT."
exit 1
set -e
SHDIR="$(dirname "$(readlink -f "${0}")")"
if which SiLiCO &>> /dev/null; then
    exec SiLiCO "${@}"
fi
if ! which conda &>> /dev/null; then
    echo "conda not found!"
    exit 127
fi

if ! conda env list | grep ^silico &>> /dev/null; then
    conda create -y -n silico -c conda-forge -c bioconda python=2.7.15 numpy=1.16.5 pybedtools=0.8.2 pysam=0.8.4
fi
eval "$(conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate silico
if "${CONDA_PREFIX}"/bin/pip freeze | grep ^SiLiCO &>> /dev/null; then
    exec "${CONDA_PREFIX}"/bin/python -m SiLiCO "${@}"
fi
"${CONDA_PREFIX}"/bin/pip install "${SHDIR}"/../../../../dist/SiLiCO-1.0.1-py2-none-any.whl
exec "${CONDA_PREFIX}"/bin/python -m SiLiCO "${@}"
