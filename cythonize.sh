#!/usr/bin/env bash
while read -r line; do
    cythonize -i -a "${line}" &
done <<< "$(find . | grep -v venv | grep -e '\.pyx$')"
wait
