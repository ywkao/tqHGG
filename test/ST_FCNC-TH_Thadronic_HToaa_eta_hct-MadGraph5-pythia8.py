from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_('General')
config.General.requestName = FlashGG_'ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8-190401-v1'
config.General.workArea = 'yuwei_crab'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../MyFlashggPlugins/flashggAnalysisNtuplizer/test/flashggAnalysisNtuplizerWithSyst_cfg.py'
config.JobType.pyCfgParams = ['processType=sig', 'filename=ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8', 'doHTXS=False', 'doSystematics=False', 'runMiniAOD=True']
config.JobType.maxMemoryMB = 2000
config.JobType.maxJobRuntimeMin = 2000

config.section_('Data')
config.Data.inputDataset = '/ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/ykao/T2CH'
config.Data.outputDatasetTag = FlashGG_'ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8-190401-v1'

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN' #if you have eos
