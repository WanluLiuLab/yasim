# shellcheck shell=bash
# shellcheck disable=SC2034
#===============================================================================
# Copyright (C) 2021-2022. tetgs authors
#
# This file is a part of tetgs, which is licensed under MIT,
# a copy of which can be obtained at <https://opensource.org/licenses/MIT>.
#
# NAME: libstr.sh -- Shell string library.
#
# VERSION HISTORY:
# 2021-07-14 1.1  : Migrated from LinuxMiniPrograms
# 2021-07-14 1.1  : TTY detection added.
#
#===============================================================================

# ANSI color is supported by most terminals, but should not appear in files.
# Check for color support
RED=""
GREEN=""
YELLOW=""
NOCOLOR=""
HAVE_COLOR=0
COLORS="$(tput colors 2>/dev/null)" || COLORS=0
# shellcheck disable=SC2086
[ ${COLORS} -gt 2 ] && HAVE_COLOR=1
[[ "${TERM:-}" =~ "256" ]] && HAVE_COLOR=1
if [ ${HAVE_COLOR} -eq 1 ]; then
    RED="\033[31m"
    GREEN="\033[32m"
    YELLOW="\033[33m"
    NOCOLOR="\033[0m"
fi

trimstr() {
    : "${1#"${1%%[![:space:]]*}"}" # Trim trailing space
    : "${_%"${_##*[![:space:]]}"}" # Trim leading space
    printf '%s\n' "${_}"
}

errh() {
    builtin echo -e "$(date) ${RED}[ERROR] ${*}${NOCOLOR}" >&2
    exit 1
}

warnh() {
    builtin echo -e "$(date) ${RED}[WARNING] ${*}${NOCOLOR}" >&2
}

infoh() {
    builtin echo -e "$(date) ${YELLOW}[INFO] ${*}${NOCOLOR}" >&2
}

debugh() {
    builtin echo -e "$(date) ${GREEN}[DEBUG] ${*}${NOCOLOR}" >&2
}
