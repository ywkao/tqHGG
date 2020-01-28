#!/usr/bin/env python2
# default only one command to be executed.
# if you want to do more than one, modify this file
# usage:
#  in command line
#    ./.submitJOB.py --command='cmsRun myJob.py && doSomeThing' --name="myQJobSent"
#
#  in python ( to use for loop to submit jobs )
#    os.system('./submitJOB.py --command={} --name={}'.format(command.format(i,i,i,i), 'runMCjob_{:02d}'.foramt(i))
#
# note, you need to check if you have such directory to write qjob summary!
defaultStorageFolder='/home/ykao/qjob/qSubResult/'
defaultMessageFolder='/home/ykao/qjob/qSubMessage/'
defaultErrorFolder  ='/home/ykao/qjob/qSubMessage/'
USER='ykao'

###############################
# template file to be submitted
###############################
submitSample='''
#!/usr/bin/env sh

#PBS -V
#PBS -j oe
#PBS -q cms
#PBS -d {0}
#PBS -o {1}
#PBS -e {2}

# setup CMSSW env
cd {3} && eval `scramv1 runtime -sh` && cd {4}
# commands to be run
{5}
'''


######################
# get CMSSW envirnment
######################
def getCMSSWVersion():
    import commands
    baseDir=commands.getstatusoutput('echo $CMSSW_BASE')
    if baseDir[0]:
        print "you haven't set CMSSW envirnment"
        exit()
    print 'current CMSSW env: {0}'.format(baseDir[1])
    return '{}/src'.format(baseDir[1])

def getCurrentPath():
    import commands
    return commands.getoutput('echo $PWD')



########################
# add parser to the code
########################
def addOption():
    import argparse
    parser = argparse.ArgumentParser(description='Hiiii')
    parser.add_argument(
            '--command', '-c', type=str, default=None,
            help='put your commands in shell'
            )
    parser.add_argument(
            '--name', '-N', type=str, default='job',
            help='job name of qSub'
            )
    parser.add_argument(
            '--user', '-u', type=str, default=USER,
            help='set user name'
            )
    parser.add_argument(
            '--lowPriority', action='store_true',
            help=''
            )

    return parser.parse_args()

if __name__ == "__main__":
    args=addOption()
    if not args.command:
        print 'you need add options, use [--help]'
        exit()

    cmsEnv=getCMSSWVersion()
    pwd=getCurrentPath()

    file=open( "/tmp/{}/tmpSh.sh".format(args.user), "w" )
    file.write( submitSample.format(defaultStorageFolder, defaultMessageFolder, defaultErrorFolder, cmsEnv, pwd, args.command) )
    file.close()
    import os
    if args.lowPriority:
        os.system( "qsub /tmp/{}/tmpSh.sh -N {} -p -500".format(args.user, args.name) )
    else:
        os.system( "qsub /tmp/{}/tmpSh.sh -N {} -p 1".format(args.user, args.name) )
