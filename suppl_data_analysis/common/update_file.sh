#!/usr/bin/env bash

rsync -av . "$(cat ssh_target.conf)"
