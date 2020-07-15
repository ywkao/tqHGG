#!/bin/bash
set -e
dir="/home/ykao/qjob/qSubMessage"
file_latest=${dir}/`ls -rt ${dir} | grep "run" | tail -n1 | head -n1`
less ${file_latest}
