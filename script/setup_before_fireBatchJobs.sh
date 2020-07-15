#!/bin/bash
set -e
#--------------- Setup ---------------#
echo "[MESSAGE] Check directories for plots..." && ./script/mkplotdir.sh
echo "[MESSAGE] Check executable..." && make
echo "[MESSAGE] Ready!"
