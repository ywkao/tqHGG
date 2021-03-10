#!/usr/bin/env python
import os, sys
import subprocess
import parallel_utils

def init():
    src_file = 'src/generalChiSquareStudy_%s.cpp' % (channel)
    exe_file = 'src/generalChiSquareStudy_%s_exe.cpp' % (channel)
    subprocess.call('cp %s %s' % (src_file, exe_file), shell=True)
    subprocess.call('make', shell=True)
    subprocess.call("mkdir -p result_top_reco_study/hist_factory_hadronic", shell=True)

channel = ["hadronic", "leptonic"]
channel = ["hadronic"]

tags = ["tt_hut", "tt_hct", "st_hut", "st_hct"]
tags = ["tt_hut"]

ext = "kapibara"
log = "result_top_reco_study/log_hadronic.txt"

command_list = []
for tag in tags:
    print tag
    exe = "./bin/generalChiSquareStudy_%s_exe" % channel
    output = "result_top_reco_study/chi2study_%s_%s_%s.root" % (channel, tag, ext)
    command = 'time %s %s %s %s' % (exe, tag, output, dataset)
    command_list.append(command)
    #command_list.append('time ./script/exe_generalChiSquareStudy.sh %s %s %s | tee %s' % (tag, channel, ext, log))

parallel_utils.submit_jobs(command_list, 24)


