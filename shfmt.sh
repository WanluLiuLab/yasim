#!/usr/bin/env bash

git ls-files | grep '\.sh$' | while read line; do
    shfmt -i 4 --write "${line}"
done
