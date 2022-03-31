#!/usr/bin/env bash
# Shell script that counts how many lines of code was written by us.
# shellcheck disable=SC2086
SOURCES=$(
    git ls-files |\
    grep -v '.maint' |\
    grep -v '.idea/' |\
    xargs
)

if which scc &>/dev/null;then
    scc ${SOURCES}
elif which cloc &>/dev/null;then
    cloc ${SOURCES}
else
    echo "scc or cloc required!"
fi
