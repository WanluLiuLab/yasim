# shellcheck shell=bash
#===============================================================================
# Copyright (C) 2021-2022. tetgs authors
#
# This file is a part of tetgs, which is licensed under MIT,
# a copy of which can be obtained at <https://opensource.org/licenses/MIT>.
#
# NAME: libisopt.sh -- Shell option-handling library.
#
# VERSION HISTORY:
# 2021-07-14 1.3  : Migrated from LinuxMiniPrograms.
# 2021-07-14 1.3  : Regular expression bugs fixed.
#
#===============================================================================

isopt() {
    case "${1:-}" in
    -? | --* | -?=*)
        return 0
        ;;
    *)
        return 1
        ;;
    esac
}
