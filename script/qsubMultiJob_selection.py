#!/usr/bin/env python2
import os
import sys

year = sys.argv[1]
channel = sys.argv[2]
# pickup list_samples of specified year{{{
list_samples=""
if year == "2016":
    list_samples = "lists/list_selection_2016"
elif year == "2017":
    list_samples = "lists/list_selection"
elif year == "2018":
    list_samples = "lists/list_selection"
elif year == "161718":
    list_samples = "lists/list_selection"
else: # year == 2017old
    list_samples="./lists/ListRootFiles"

#}}}
print year, list_samples

command='"cd /wk_cms2/ykao/CMSSW_9_4_10/src/2017/tqHGG; ./script/exe_selection_batch.sh {0} {1} {2} {3}"'

commentOut="#"
ext="root"
counter = 0
file = open(list_samples, "r")
for x in file:
    if commentOut in x: continue
    counter = counter + 1 
    tagName = 'selection_{0}_{1}'.format(year, counter)
    # identify the type of the file (isDirectory){{{
    isDirectory=""
    if ext in x:
        isDirectory="rootfile"
    else:
        isDirectory="directory"
    #}}}
    #--------------------------------------------------
    print tagName, command.format(x.rstrip(), year, isDirectory, channel)
    os.system('./script/submitJOB.py --command={} --name={}'.format(command.format(x.rstrip(), year, isDirectory, channel), tagName))
