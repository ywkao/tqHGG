# tqHGG
t to cH analysis

Searches for top FCNH with H decays to diphoton analysis.

Although the analysis code is performed under CMSSW_9_4_10,
it should be independent of the version of CMSSW because the only used tool is ROOT.

Usage:
1. cmsenv #In order to use ROOT. 
2. Modify the directory of input root files in the shell script, exePreselection.
3. ./exePreselection #Perform preselection ~7hrs. Additional option: -d or --dryRun for testing purpose. 
4. ./doSelection #Perform selection ~4min.
5. root -l -b -q src/stackHist.C #Make stack plots.
