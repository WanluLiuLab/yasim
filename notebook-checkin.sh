#!/usr/bin/env bash
set -eu
SHDIR="$(readlink -f "$(dirname "${0}")")"
cd "${SHDIR}" || exit 1

find . | grep '\.ipynb$' | grep -v ipynb_checkpoints | grep -v _build | while read -r fn; do
    target_fn="${fn/.ipynb/.ipynb.md}"
    if [ -e "${target_fn}" ] && [ "${target_fn}" -nt "${fn}" ]; then
        echo "CONVERT ${fn} -> ${target_fn}: REFUSE TO OVERWRITE NEWER FILE"
        continue
    fi
    echo "CONVERT ${fn} -> ${target_fn}: START"
    if jupytext --to md:myst --from notebook --output "${target_fn}" "${fn}" &>>nbconvert.log; then
        echo "CONVERT ${fn} -> ${target_fn}: FIN"
    else
        echo "CONVERT ${fn} -> ${target_fn}: ERR"
    fi
done
