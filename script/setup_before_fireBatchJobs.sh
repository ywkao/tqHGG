#!/bin/bash
set -e
#--------------- Setup ---------------#
echo "[MESSAGE] Check directories for plots..." && ./script/mkplotdir.sh
echo "[MESSAGE] cp src/preselection.cpp src/preselection_exe.cpp" && cp src/preselection.cpp src/preselection_exe.cpp
echo "[MESSAGE] Check executable..." && make
echo "[MESSAGE] Ready!"
