#!/usr/bin/env make
#===============================================================================
# Copyright (C) 2021-2022. tetgs authors
#
# This file is a part of tetgs, which is licensed under MIT,
# a copy of which can be obtained at <https://opensource.org/licenses/MIT>.
#
# NAME: Makefile -- Makefile of building sphinx documentations
#
# VERSION HISTORY:
# 2021-08-20 0.1  : Purposed and added by YU Zhejian.
#
#===============================================================================

SPHINX_OPTS := -j auto

SPHINX_APIDOC_ARGS := --separate \
	--no-toc \
	--module-first \
	--maxdepth 2 \
	--force \
	--output-dir $(CURDIR)/_apidoc \
	--doc-project 'API Docs'

BUILD_DIR := $(CURDIR)/_build

.PHONY: all

all: pdf html

.PHONY: apidoc
apidoc:
	sphinx-apidoc $(SPHINX_APIDOC_ARGS) $(PWD)/src/bioutils
	sphinx-apidoc $(SPHINX_APIDOC_ARGS) $(PWD)/src/commonutils
	sphinx-apidoc $(SPHINX_APIDOC_ARGS) $(PWD)/src/yasim

html:
	sphinx-build -M html $(CURDIR) $(BUILD_DIR) $(SPHINX_OPTS)

pdf: apidoc
	sphinx-build -M latexpdf $(CURDIR) $(BUILD_DIR) $(SPHINX_OPTS)

clean:
	$(RM) $(BUILD_DIR)

distclean:
