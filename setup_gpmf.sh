#!/usr/bin/env bash
#===============================================================================
# Copyright (C) 2021-2022. gpmf authors
#
# This file is a part of tetgs, which is licensed under MIT,
# a copy of which can be obtained at <https://opensource.org/licenses/MIT>.
#
# NAME: setup.sh -- Setup gpmf environment
#
# VERSION HISTORY:
# 2021-10-26 0.1  : Purposed and added by YU Zhejian.
#
#===============================================================================
set -ue
SHDIR="$(dirname "$(readlink -f "${0}")")"
cd "${SHDIR}" || exit 1
[ -d .maint ] || submodule add https://gitee.com/yuzjlab/gpmf .maint
bash .maint/setup.sh
