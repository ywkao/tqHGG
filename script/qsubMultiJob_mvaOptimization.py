#!/usr/bin/env python2
command='"cd /wk_cms2/ykao/CMSSW_9_4_10/src/2017/tqHGG; time root -l -b -q {} | tee {}"'
import os

list = [
#"src/tmp/testAll_24_v1/TMVAClassification.C",
#"src/tmp/testAll_24_v2/TMVAClassification.C",
#"src/tmp/testAll_24_v3/TMVAClassification.C",
"src/tmp/testAll_24_v4/TMVAClassification.C",
"src/tmp/testAll_24_v5/TMVAClassification.C",
"src/tmp/testAll_24_v6/TMVAClassification.C",
"src/tmp/testAll_24_v7/TMVAClassification.C",
"src/tmp/testAll_24_v8/TMVAClassification.C",
"src/tmp/testAll_24_v9/TMVAClassification.C",
"src/tmp/testAll_24_v10/TMVAClassification.C",
"src/tmp/testAll_24_v11/TMVAClassification.C",
"src/tmp/testAll_24_v12/TMVAClassification.C",
"src/tmp/testAll_24_v13/TMVAClassification.C",
"src/tmp/testAll_24_v14/TMVAClassification.C",
"src/tmp/testAll_24_v15/TMVAClassification.C"
]

counter = 0
for dataset in list:
    counter = counter + 1
    tagName = 'mvaOpt_{}'.format(counter)
    log = 'log/log_{}.txt'.format(tagName)
    #--------------------------------------------------
    #os.system('echo --command={} --name={}'.format(command.format(dataset, log), tagName))
    os.system('./script/submitJOB.py --command={} --name={}'.format(command.format(dataset, log), tagName))
