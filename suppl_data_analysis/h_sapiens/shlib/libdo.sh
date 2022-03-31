# shellcheck shell=bash
# shellcheck disable=SC2015
# shellcheck disable=SC2086
#===============================================================================
# Copyright (C) 2021-2022. tetgs authors
#
# This file is a part of tetgs, which is licensed under MIT,
# a copy of which can be obtained at <https://opensource.org/licenses/MIT>.
#
# NAME: libdo.sh -- The logger and monitor wrapper.
#
# VERSION HISTORY:
# 2021-07-14 2.1  : Migrated from LinuxMiniPrograms.
# 2021-08-14 2.1  : Function prefixed with '__' to avoid misleading.
#
#===============================================================================

__DO_ECHO() {
    [ "${LIBDO_LOG_MODE:-}" = "S" ] && return || true
    [ -n "${LIBDO_LOG:-}" ] && echo "${*}" >>"${LIBDO_LOG}" || echo -e "\033[33m${*}\033[0m" >&2
}
__DO() {
    local LIBDO_CMD="${*}"
    DO_ECHO "LIBDO IS GOING TO EXECUTE ${LIBDO_CMD}"
    DO_ECHO "LIBDO STARTED AT $(date "+%Y-%m-%d %H:%M:%S")"
    local LIBDO_PID
    if [ -z "${LIBDO_LOG:-}" ]; then
        eval "${LIBDO_CMD}" &
        LIBDO_PID=${!}
    else
        case "${LIBDO_LOG_MODE:-}" in
        "2")
            eval ${LIBDO_CMD} 2>>"${LIBDO_LOG}" &
            LIBDO_PID=${!}
            ;;
        "3")
            eval ${LIBDO_CMD} >>"${LIBDO_LOG}" &
            LIBDO_PID=${!}
            ;;
        "4")
            eval ${LIBDO_CMD} &>>"${LIBDO_LOG}" &
            LIBDO_PID=${!}
            ;;
        *)
            eval ${LIBDO_CMD} &
            LIBDO_PID=${!}
            ;;
        esac
    fi
    [ -z "${LIBDO_TOP_PID:-}" ] && DO_ECHO "LIBDO PID ${LIBDO_PID}" || DO_ECHO "LIBDO PID ${LIBDO_PID} with top_pid ${LIBDO_TOP_PID}"
    wait ${LIBDO_PID} && LIBDO_PRIV=0 || LIBDO_PRIV=${?}
    DO_ECHO "LIBDO STOPPED AT $(date "+%Y-%m-%d %H:%M:%S")"
    if [ ${LIBDO_PRIV} -ne 0 ]; then
        DO_ECHO "LIBDO FAILED, GOT \$?=${LIBDO_PRIV}"
        if [ -n "${LIBDO_TOP_PID:-}" ]; then
            DO_ECHO "LIBDO WILL KILL ${LIBDO_TOP_PID}"
            kill -9 "${LIBDO_TOP_PID}"
        fi
    else
        DO_ECHO "LIBDO EXITED SUCCESSFULLY"
    fi
    return ${LIBDO_PRIV}
}
