#!/bin/env  python
#*************************************************************************
#
#  Filename    : makePUweight.py
#  Description : Making pileup weights using configuration
#  Author      : Yi-Mu "Enoch" Chen [ ensc@hep1.phys.ntu.edu.tw ]
#
#*************************************************************************
import argparse
import importlib
import os

import ROOT


def main():
    parser = argparse.ArgumentParser(
        description='Creating the Pileupweigt in plain text files')
    parser.add_argument('-x', '--crosssection', dest='xsec', type=int,
                        help='Minimum bias cross section to use',
                        default=None,
                        )
    parser.add_argument('-f', '--fileprefix', dest='prefix', type=str,
                        help='Storage prefix',
                        default='./data/pileupweights',
                        )
    parser.add_argument('-m', '--mcmodule', dest='mcmodule', type=str,
                        help='MC Module to use from SimGeneral.MixingModule',
                        default='mix_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU_cfi',
                        #default='mix_2017_25ns_UltraLegacy_PoissonOOTPU_cfi',
                        #default='mix_2016_25ns_Moriond17MC_PoissonOOTPU_cfi',
                        )
    parser.add_argument('-l', '--lumimask', dest='lumimask', type=str,
                        help='Lumimask json file to use',
                        default='data/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt',
                        #default='data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt',
                        )
    parser.add_argument('-p', '--pufile', dest='pufile', type=str,
                        help='Pileup json file to use',
                        default='data/pileup_latest.txt',
                        #default='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt',
                        )
    args = parser.parse_args()

    if not args.xsec:
        print "Minimum bias cross section input required"
        parser.print_help()
        return 1

    #-------------------------------------------------------------------------
    #   Getting MC Pileup from CMSSW python modules
    #-------------------------------------------------------------------------
    puset = importlib.import_module('SimGeneral.MixingModule.' + args.mcmodule)
    mcpileup = puset.mix.input.nbPileupEvents.probValue

    #-------------------------------------------------------------------------
    #   Getting Data Pile up
    #-------------------------------------------------------------------------
    os.system( """
        pileupCalc.py       \
        -i {0}              \
        --inputLumiJSON {1} \
        --calcMode true     \
        --minBiasXsec {2}   \
        --maxPileupBin  {3} \
        --numPileupBins {3} \
        ./data/PileUp.root
    """.format(args.lumimask, args.pufile, args.xsec, 75 ) ) #because the .txt bin number is 75
    datapufile = ROOT.TFile.Open('./data/PileUp.root')
    datapuhist = datapufile.Get('pileup')
    datapuhist.Scale(1. / datapuhist.Integral())

    #-------------------------------------------------------------------------
    #   Calculating Pile up weights
    #-------------------------------------------------------------------------
    mcweight = []
    orgweightsum = 0
    puweightsum = 0
    mchist   = ROOT.TH1D("mcpu","mcpu",len(mcpileup),0,len(mcpileup))
    for i in range(0, len(mcpileup)):
        mcweight.append(datapuhist.GetBinContent(i + 1) / mcpileup[i])
        mchist.SetBinContent(i+1,mcpileup[i])
        orgweightsum = orgweightsum + mcpileup[i]
        puweightsum = puweightsum + (mcpileup[i] * mcweight[i])
        print i, mcweight[i], mcpileup[i], datapuhist.GetBinContent(i + 1)

    mcfile = ROOT.TFile.Open("./data/MCPileUp.root","update")
    mchist.Write()
    print puweightsum
    print orgweightsum
    print len(mcweight)
    # os.system('rm -rf /tmp/PileUp.root')

    #-------------------------------------------------------------------------
    #   Saving to file
    #-------------------------------------------------------------------------
    outfileformat = "{0}_{1}.csv"
    outfile = open(outfileformat.format(args.prefix, args.xsec), "w")
    outstring = '\n'.join(str(x) for x in mcweight)
    outfile.write(outstring + '\n')
    outfile.close()


if __name__ == "__main__":
    main()
