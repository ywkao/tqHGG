#!/usr/bin/env python2

# test{{{
#---------- test ----------#
#command='"cd /wk_cms2/ykao/CMSSW_9_4_10/src/2017/tqHGG/script && ./hello.py"'
#import os
#os.system('./script/submitJOB.py --command={} --name={}'.format(command, 'runMVAjob'))
#}}}

command='"cd /wk_cms2/ykao/CMSSW_9_4_10/src/2017/tqHGG; ./script/mvaTasks/mvaTraining.sh {1}"'
import os
mvaNum = 7
for i in range(1, mvaNum+1):
    tagName = 'testAll_{:02d}'.format(i)

    # edit src/TMVAClassification_leptonic.C && cp src/TMVAClassification_leptonic.C src/tmp/task_${tag}/TMVAClassification.C
    os.system('./script/mvaTasks/mvaTask_{0:02d}.sh {1}'.format(i, tagName))

    # submit jobs; root -l -b -q src/TMVAClassification.C
    os.system('./script/submitJOB.py --command={} --name={}'.format(command.format(i, tagName), 'runMVAjob_{:02d}'.format(i)))

# legacy{{{
#command='"cd /wk_cms2/ykao/CMSSW_9_4_10/src/2017/tqHGG/script && ./hello.py"'
#command='"cd /wk_cms/ltsai/LbFrame/MC/lotofRun && cmsRun BPH-RunIISpring16DR80-00058_1_cfg.py runNum={} fileList=fLink{:02d} && cmsRun BPH-RunIISpring16DR80-00058_2_cfg.py runNum={} && /bin/rm BPH-RunIISpring16DR80-00058_step1_{:02d}.root"'
#runNumberFrom=20
#runNumberTo=30
#import os
#
## use for specificFileName
#spNum=[20]
#for i in spNum:
##for i in range(runNumberFrom, runNumberTo):
#    os.system('./submitJOB.py --command={} --name={}'.format(command.format(i,i,i,i), 'runMCjob_{:02d}'.format(i)))
#}}}
