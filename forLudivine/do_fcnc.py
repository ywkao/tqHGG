# vim: set fdm=marker:
# imports, parser{{{
import sys, os
sys.path.append("../")
import itertools
import parallel_utils
import workflow_utils
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--tag", help = "tag to denote with", type=str)
parser.add_argument("--baby_version", help = "which version of babies to use", type=str)
args = parser.parse_args()
#}}}
# bdt.xml files{{{
#bdt_lep = "../MVAs/Leptonic_4June2019_v1.7_FCNC_bdt.xml"
#bdt_had = "../MVAs/Hadronic_12June2019_v1.7_impute_FCNC_bdt.xml"

#bdt_lep_hut = "../MVAs/Leptonic_v1.7_27Jun2019_FCNC_SingleBDT_hut__bdt.xml"
#bdt_lep_hct = "../MVAs/Leptonic_v1.7_27Jun2019_FCNC_SingleBDT_hct__bdt.xml"
#bdt_had_hut = "../MVAs/Hadronic_v1.7_27Jun2019_FCNC_SingleBDT_impute_hut__bdt.xml"
#bdt_had_hct = "../MVAs/Hadronic_v1.7_27Jun2019_FCNC_SingleBDT_impute_hct__bdt.xml"

bdt_lep_hut = "../MVAs/Leptonic_v4.11_14Jan2020_hut__bdt.xml"
bdt_lep_hct = "../MVAs/Leptonic_v4.11_14Jan2020_hct__bdt.xml"
bdt_had_hut = "../MVAs/Hadronic_v4.11_14Jan2020_impute_hut__bdt.xml"
bdt_had_hct = "../MVAs/Hadronic_v4.11_14Jan2020_impute_hct__bdt.xml"
#}}}
from subprocess import *
os.chdir("../")
do_looping = True
if do_looping:
    #p = Popen("echo" + " Hello World! (from do_fcnc.py)", shell=True)
    p = Popen("pwd", shell=True)
    result = os.waitpid(p.pid, 0)

    # Loopers
    #parallel_utils.run('python looper_wrapper.py --fcnc --channel "Leptonic" --baby_version "%s" --tag "%s" --selection "ttHLeptonic_RunII_MVA_Presel" --bkg_options "none" --years "2017"' % (args.baby_version, args.tag + "_hut_BDT"))
    #parallel_utils.run('python looper_wrapper.py --fcnc --channel "Hadronic" --baby_version "%s" --tag "%s" --selection "ttHHadronic_RunII_MVA_Presel" --bkg_options "impute" --years "2017"' % (args.baby_version, args.tag + "_impute_hut_BDT"))
    #parallel_utils.run('python looper_wrapper.py --fcnc --channel "Leptonic" --baby_version "%s" --tag "%s" --selection "ttHLeptonic_RunII_MVA_Presel" --bkg_options "none" --years "2017"' % (args.baby_version, args.tag + "_hct_BDT"))
    #parallel_utils.run('python looper_wrapper.py --fcnc --channel "Hadronic" --baby_version "%s" --tag "%s" --selection "ttHHadronic_RunII_MVA_Presel" --bkg_options "impute" --years "2017"' % (args.baby_version, args.tag + "_impute_hct_BDT"))

    # Data/MC Plots
    os.chdir("Plots")
    parallel_utils.run('python plot_wrapper.py --input_file "../%s" --plot_type "std_2017" --plot_labels "FCNC Leptonic|Loose MVA Presel." --signals "TT_FCNC_hut|ST_FCNC_hut" --backgrounds "DiPhoton|GammaJets|TTGG|TTGJets|TTJets|VG"' % ("ttHLeptonic_RunII_MVA_Presel_%s_histogramsRunII.root" % (args.tag + "_hut_BDT_FCNC")))
    parallel_utils.run('python plot_wrapper.py --input_file "../%s" --plot_type "std_2017" --plot_labels "FCNC Leptonic|Loose MVA Presel." --signals "TT_FCNC_hct|ST_FCNC_hct" --backgrounds "DiPhoton|GammaJets|TTGG|TTGJets|TTJets|VG"' % ("ttHLeptonic_RunII_MVA_Presel_%s_histogramsRunII.root" % (args.tag + "_hct_BDT_FCNC")))
    parallel_utils.run('python plot_wrapper.py --input_file "../%s" --plot_type "std_2017" --plot_labels "FCNC Hadronic|Loose MVA Presel." --signals "TT_FCNC_hut|ST_FCNC_hut" --backgrounds "DiPhoton|QCD_GammaJets_imputed|TTGG|TTGJets|TTJets|VG"' % ("ttHHadronic_RunII_MVA_Presel_%s_histogramsRunII.root" % (args.tag + "_impute_hut_BDT_FCNC")))
    parallel_utils.run('python plot_wrapper.py --input_file "../%s" --plot_type "std_2017" --plot_labels "FCNC Hadronic|Loose MVA Presel." --signals "TT_FCNC_hct|ST_FCNC_hct" --backgrounds "DiPhoton|QCD_GammaJets_imputed|TTGG|TTGJets|TTJets|VG"' % ("ttHHadronic_RunII_MVA_Presel_%s_histogramsRunII.root" % (args.tag + "_impute_hct_BDT_FCNC")))
    os.chdir("../")

# Legacy code{{{
  # Loopers{{{
  #parallel_utils.run('python looper_wrapper.py --fcnc --channel "Leptonic" --baby_version "%s" --tag "%s" --selection "ttHLeptonic_RunII_MVA_Presel" --bkg_options "none" --years "2017" --bdt "%s"' % (args.baby_version, args.tag + "_hut_BDT", bdt_lep_hut))
  #parallel_utils.run('python looper_wrapper.py --fcnc --channel "Hadronic" --baby_version "%s" --tag "%s" --selection "ttHHadronic_RunII_MVA_Presel" --bkg_options "impute" --years "2017" --bdt "%s"' % (args.baby_version, args.tag + "_impute_hut_BDT", bdt_had_hut))
  #parallel_utils.run('python looper_wrapper.py --fcnc --channel "Leptonic" --baby_version "%s" --tag "%s" --selection "ttHLeptonic_RunII_MVA_Presel" --bkg_options "none" --years "2017" --bdt "%s"' % (args.baby_version, args.tag + "_hct_BDT", bdt_lep_hct))
  #parallel_utils.run('python looper_wrapper.py --fcnc --channel "Hadronic" --baby_version "%s" --tag "%s" --selection "ttHHadronic_RunII_MVA_Presel" --bkg_options "impute" --years "2017" --bdt "%s"' % (args.baby_version, args.tag + "_impute_hct_BDT", bdt_had_hct))
  #}}}
#  # Loopers, signal regions{{{
#  #parallel_utils.run('python looper_wrapper.py --channel "Leptonic" --baby_version "%s" --tag "%s" --selection "FCNC_Leptonic_Hut_RunII_SR_Inclusive" --bkg_options "none" --years "2017" --bdt "%s"' % (args.baby_version, args.tag + "_hut_BDT", bdt_lep_hut))
#  #parallel_utils.run('python looper_wrapper.py --channel "Hadronic" --baby_version "%s" --tag "%s" --selection "FCNC_Hadronic_Hut_RunII_SR_Inclusive" --bkg_options "impute" --years "2017" --bdt "%s"' % (args.baby_version, args.tag + "_impute_hut_BDT", bdt_had_hut))
#  #parallel_utils.run('python looper_wrapper.py --channel "Leptonic" --baby_version "%s" --tag "%s" --selection "FCNC_Leptonic_Hct_RunII_SR_Inclusive" --bkg_options "none" --years "2017" --bdt "%s"' % (args.baby_version, args.tag + "_hct_BDT", bdt_lep_hct))
#  #parallel_utils.run('python looper_wrapper.py --channel "Hadronic" --baby_version "%s" --tag "%s" --selection "FCNC_Hadronic_Hct_RunII_SR_Inclusive" --bkg_options "impute" --years "2017" --bdt "%s"' % (args.baby_version, args.tag + "_impute_hct_BDT", bdt_had_hct))
#  #}}}
#  # Data/MC Plots{{{
#  os.chdir("Plots")
#  #parallel_utils.run('python plot_wrapper.py --input_file "../%s" --plot_type "std_2017" --plot_labels "FCNC Leptonic|Loose MVA Presel." --signals "TT_FCNC_hut|ST_FCNC_hut" --backgrounds "DiPhoton|GammaJets|TTGG|TTGJets|TTJets|VG"' % ("ttHLeptonic_RunII_MVA_Presel_%s_histogramsRunII.root" % (args.tag + "_hut_BDT_FCNC")))
#  parallel_utils.run('python plot_wrapper.py --input_file "../%s" --plot_type "std_2017" --plot_labels "FCNC Leptonic|Loose MVA Presel." --signals "TT_FCNC_hct|ST_FCNC_hct" --backgrounds "DiPhoton|GammaJets|TTGG|TTGJets|TTJets|VG"' % ("ttHLeptonic_RunII_MVA_Presel_%s_histogramsRunII.root" % (args.tag + "_hct_BDT_FCNC")))
#  
#  parallel_utils.run('python plot_wrapper.py --input_file "../%s" --plot_type "std_2017" --plot_labels "FCNC Hadronic|Loose MVA Presel." --signals "TT_FCNC_hut|ST_FCNC_hut" --backgrounds "DiPhoton|QCD_GammaJets_imputed|TTGG|TTGJets|TTJets|VG"' % ("ttHHadronic_RunII_MVA_Presel_%s_histogramsRunII.root" % (args.tag + "_impute_hut_BDT_FCNC")))
#  parallel_utils.run('python plot_wrapper.py --input_file "../%s" --plot_type "std_2017" --plot_labels "FCNC Hadronic|Loose MVA Presel." --signals "TT_FCNC_hct|ST_FCNC_hct" --backgrounds "DiPhoton|QCD_GammaJets_imputed|TTGG|TTGJets|TTJets|VG"' % ("ttHHadronic_RunII_MVA_Presel_%s_histogramsRunII.root" % (args.tag + "_impute_hct_BDT_FCNC")))
#
#  os.chdir("../")
#  #}}}
#  # MVA BabyMaker{{{
#  #parallel_utils.run('python looper_wrapper.py --channel "Leptonic" --baby_version "%s" --tag "%s" --selection "ttHLeptonic_RunII_MVA_Presel" --bkg_options "none" --years "2017" --babymaker --fcnc' % (args.baby_version, args.tag))
#  #parallel_utils.run('python looper_wrapper.py --channel "Hadronic" --baby_version "%s" --tag "%s" --selection "ttHHadronic_RunII_MVA_Presel" --bkg_options "impute" --years "2017" --babymaker --fcnc' % (args.baby_version, args.tag + "_impute"))
#  #}}}
## others{{{
##do_mvas = True
##if do_mvas:
##  os.chdir("../MVAs/")
##  # MVA Prep
##  command_list = []
##  #command_list.append('python prep.py --input "../Loopers/MVABaby_ttHLeptonic_%s_FCNC.root" --channel "Leptonic" --fcnc_hut --tag "_hut"' % (args.tag))
##  #command_list.append('python prep.py --input "../Loopers/MVABaby_ttHHadronic_%s_FCNC.root" --channel "Hadronic" --fcnc_hut --tag "_hut"' % (args.tag + "_impute")) 
##  #command_list.append('python prep.py --input "../Loopers/MVABaby_ttHLeptonic_%s_FCNC.root" --channel "Leptonic" --fcnc_hct --tag "_hct"' % (args.tag))
##  #command_list.append('python prep.py --input "../Loopers/MVABaby_ttHHadronic_%s_FCNC.root" --channel "Hadronic" --fcnc_hct --tag "_hct"' % (args.tag + "_impute"))
##  parallel_utils.submit_jobs(command_list, 4)
##
##
##  # MVA Prep, no mass constraints
##  #parallel_utils.run('python prep.py --input "../Loopers/MVABaby_ttHHadronic_%s_FCNC.root" --channel "Hadronic" --fcnc_hut --tag "_hut_no_mass_constraint" --no_mass_constraint' % (args.tag + "_impute"))
##  #parallel_utils.run('python prep.py --input "../Loopers/MVABaby_ttHHadronic_%s_FCNC.root" --channel "Hadronic" --fcnc_hct --tag "_hct_no_mass_constraint" --no_mass_constraint' % (args.tag + "_impute"))
##
##  # MVA Training
##  #parallel_utils.run('python train.py --input "ttHLeptonic_%s_FCNC_features_hut.hdf5" --channel "Leptonic" --tag "%s" --ext ""' % (args.tag, args.tag + "_hut"))
##  #parallel_utils.run('python train.py --input "ttHHadronic_%s_FCNC_features_hut.hdf5" --channel "Hadronic" --tag "%s" --ext ""' % (args.tag + "_impute", args.tag + "_impute_hut"))
##  #parallel_utils.run('python train.py --input "ttHLeptonic_%s_FCNC_features_hct.hdf5" --channel "Leptonic" --tag "%s" --ext ""' % (args.tag, args.tag + "_hct"))
##  #parallel_utils.run('python train.py --input "ttHHadronic_%s_FCNC_features_hct.hdf5" --channel "Hadronic" --tag "%s" --ext ""' % (args.tag + "_impute", args.tag + "_impute_hct")) 
##
##  # MVA Training, no mass constraint
##  #parallel_utils.run('python train.py --input "ttHHadronic_%s_FCNC_features_hut_no_mass_constraint.hdf5" --channel "Hadronic" --tag "%s" --ext ""' % (args.tag + "_impute", args.tag + "_impute_hut_no_mass_constraint"))
##  #parallel_utils.run('python train.py --input "ttHHadronic_%s_FCNC_features_hct_no_mass_constraint.hdf5" --channel "Hadronic" --tag "%s" --ext ""' % (args.tag + "_impute", args.tag + "_impute_hct_no_mass_constraint")) 
##
##  # MVA Training (multi)
##  #parallel_utils.run('python train.py --input "ttHLeptonic_%s_FCNC_features_hut.hdf5" --channel "Leptonic" --tag "%s" --ext "" --multi' % (args.tag, args.tag + "_multi" + "_hut"))
##  #parallel_utils.run('python train.py --input "ttHHadronic_%s_FCNC_features_hut.hdf5" --channel "Hadronic" --tag "%s" --ext "" --multi' % (args.tag + "_impute", args.tag + "_multi" + "_impute_hut"))
##  #parallel_utils.run('python train.py --input "ttHLeptonic_%s_FCNC_features_hct.hdf5" --channel "Leptonic" --tag "%s" --ext "" --multi' % (args.tag, args.tag + "_multi" + "_hct"))
##  #parallel_utils.run('python train.py --input "ttHHadronic_%s_FCNC_features_hct.hdf5" --channel "Hadronic" --tag "%s" --ext "" --multi' % (args.tag + "_impute", args.tag + "_multi" + "_impute_hct"))
##  #parallel_utils.run('python train.py --input "ttHHadronic_%s_FCNC_features_hct.hdf5" --channel "Hadronic" --tag "%s_EqualClassWeights" --ext "" --multi --equal_weights' % (args.tag + "_impute", args.tag + "_multi_impute_hct"))
##
##
##  # FCNC vs SM Higgs BDT/DNN
##  command_list = []
##  #command_list.append('python prep_dnn.py --input "../Loopers/MVABaby_ttHLeptonic_%s_FCNC.root" --channel "Leptonic" --fcnc --z_score --signal "FCNC_hut" --backgrounds "ttH,tH" --tag "hut_FCNC_vs_SMHiggs"' % (args.tag))
##  #command_list.append('python prep.py --input "../Loopers/MVABaby_ttHLeptonic_%s_FCNC.root" --channel "Leptonic" --fcnc_hut --FCNC_vs_SMHiggs --tag "_hut_FCNC_vs_SMHiggs"' % (args.tag))
##  #command_list.append('python prep_dnn.py --input "../Loopers/MVABaby_ttHLeptonic_%s_FCNC.root" --channel "Leptonic" --fcnc --z_score --signal "FCNC_hct" --backgrounds "ttH,tH" --tag "hct_FCNC_vs_SMHiggs"' % (args.tag))
##  #command_list.append('python prep.py --input "../Loopers/MVABaby_ttHLeptonic_%s_FCNC.root" --channel "Leptonic" --fcnc_hct --FCNC_vs_SMHiggs --tag "_hct_FCNC_vs_SMHiggs"' % (args.tag))
##
##  #command_list.append('python prep_dnn.py --input "../Loopers/MVABaby_ttHHadronic_%s_FCNC.root" --channel "Hadronic" --fcnc --z_score --signal "FCNC_hut" --backgrounds "ttH,tH,ggH" --tag "hut_FCNC_vs_SMHiggs_addMassConstraint" --add_mass_constraints' % (args.tag + "_impute"))
##  ##command_list.append('python prep_dnn.py --input "../Loopers/MVABaby_ttHHadronic_%s_FCNC.root" --channel "Hadronic" --fcnc --z_score --signal "FCNC_hut" --backgrounds "ttH,tH,ggH" --tag "hut_FCNC_vs_SMHiggs_addExtMassConstraint" --add_mass_constraints --add_ext_mass_constraints' % (args.tag + "_impute"))
##  ###command_list.append('python prep_dnn.py --input "../Loopers/MVABaby_ttHHadronic_%s_FCNC.root" --channel "Hadronic" --fcnc --z_score --signal "FCNC_hut" --backgrounds "ttH,tH" --tag "hut_FCNC_vs_SMHiggs_addExtMassConstraint" --add_mass_constraints --add_ext_mass_constraints' % (args.tag + "_impute"))
##  #command_list.append('python prep_dnn.py --input "../Loopers/MVABaby_ttHHadronic_%s_FCNC.root" --channel "Hadronic" --fcnc --z_score --signal "FCNC_hut" --backgrounds "ttH,tH,ggH" --tag "hut_FCNC_vs_SMHiggs"' % (args.tag + "_impute"))
##  #  # command_list.append('python prep_dnn.py --input "../Loopers/MVABaby_ttHHadronic_%s_FCNC.root" --channel "Hadronic" --fcnc --z_score --signal "FCNC_hut" --backgrounds "ttH,tH,ggH" --tag "hut_FCNC_vs_SMHiggs_dontTrainWithGGH" --ggH_treatment "dont_train"' % (args.tag + "_impute"))
##    # command_list.append('python prep_dnn.py --input "../Loopers/MVABaby_ttHHadronic_%s_FCNC.root" --channel "Hadronic" --fcnc --z_score --signal "FCNC_hut" --backgrounds "ttH,tH,ggH" --tag "hut_FCNC_vs_SMHiggs_scaleGGH" --ggH_treatment "scale"' % (args.tag + "_impute"))
##  #command_list.append('python prep.py --input "../Loopers/MVABaby_ttHHadronic_%s_FCNC.root" --channel "Hadronic" --fcnc_hut --FCNC_vs_SMHiggs --tag "_hut_FCNC_vs_SMHiggs"' % (args.tag + "_impute"))
##  #command_list.append('python prep_dnn.py --input "../Loopers/MVABaby_ttHHadronic_%s_FCNC.root" --channel "Hadronic" --fcnc --z_score --signal "FCNC_hct" --backgrounds "ttH,tH,ggH" --tag "hct_FCNC_vs_SMHiggs_addMassConstraint" --add_mass_constraints' % (args.tag + "_impute"))
##  #command_list.append('python prep_dnn.py --input "../Loopers/MVABaby_ttHHadronic_%s_FCNC.root" --channel "Hadronic" --fcnc --z_score --signal "FCNC_hct" --backgrounds "ttH,tH,ggH" --tag "hct_FCNC_vs_SMHiggs_addExtMassConstraint" --add_mass_constraints --add_ext_mass_constraints' % (args.tag + "_impute"))
##  #command_list.append('python prep_dnn.py --input "../Loopers/MVABaby_ttHHadronic_%s_FCNC.root" --channel "Hadronic" --fcnc --z_score --signal "FCNC_hct" --backgrounds "ttH,tH,ggH" --tag "hct_FCNC_vs_SMHiggs"' % (args.tag + "_impute"))
##  #command_list.append('python prep.py --input "../Loopers/MVABaby_ttHHadronic_%s_FCNC.root" --channel "Hadronic" --fcnc_hct --FCNC_vs_SMHiggs --tag "_hct_FCNC_vs_SMHiggs"' % (args.tag + "_impute"))
##  parallel_utils.submit_jobs(command_list, 8)
##
##  # Training
##  #parallel_utils.run('python train.py --input "ttHLeptonic_%s_FCNC_features_hut_FCNC_vs_SMHiggs.hdf5" --channel "Leptonic" --tag "%s" --ext ""' % (args.tag, args.tag + "_hut_FCNC_vs_SMHiggs"))
##  #parallel_utils.run('python train_dnn.py --input "ttHLeptonic_%s_FCNC_dnn_features_hut_FCNC_vs_SMHiggs.hdf5" --channel "Leptonic" --tag "%s" --preprocess_scheme "preprocess_scheme_Leptonic_hut_FCNC_vs_SMHiggs.json"' % (args.tag, args.tag + "_hut_FCNC_vs_SMHiggs"), False)
##  #parallel_utils.run('python train.py --input "ttHLeptonic_%s_FCNC_features_hct_FCNC_vs_SMHiggs.hdf5" --channel "Leptonic" --tag "%s" --ext ""' % (args.tag, args.tag + "_hct_FCNC_vs_SMHiggs"))
##  #parallel_utils.run('python train_dnn.py --input "ttHLeptonic_%s_FCNC_dnn_features_hct_FCNC_vs_SMHiggs.hdf5" --channel "Leptonic" --tag "%s" --preprocess_scheme "preprocess_scheme_Leptonic_hct_FCNC_vs_SMHiggs.json"' % (args.tag, args.tag + "_hct_FCNC_vs_SMHiggs"), False)
##
##  #parallel_utils.run('python train.py --input "ttHHadronic_%s_FCNC_features_hut_FCNC_vs_SMHiggs.hdf5" --channel "Hadronic" --tag "%s" --ext ""' % (args.tag + "_impute", args.tag + "_impute_hut_FCNC_vs_SMHiggs"))
##  #parallel_utils.run('python train_dnn.py --input "ttHHadronic_%s_FCNC_dnn_features_hut_FCNC_vs_SMHiggs.hdf5" --channel "Hadronic" --tag "%s" --preprocess_scheme "preprocess_scheme_Hadronic_hut_FCNC_vs_SMHiggs.json"' % (args.tag + "_impute", args.tag + "_impute_hut_FCNC_vs_SMHiggs"), False)
##    #parallel_utils.run('python train_dnn.py --input "ttHHadronic_%s_FCNC_dnn_features_hut_FCNC_vs_SMHiggs_dontTrainWithGGH.hdf5" --channel "Hadronic" --tag "%s"' % (args.tag + "_impute", args.tag + "_impute_hut_FCNC_vs_SMHiggs_dontTrainWithGGH"))
##    #parallel_utils.run('python train_dnn.py --input "ttHHadronic_%s_FCNC_dnn_features_hut_FCNC_vs_SMHiggs_scaleGGH.hdf5" --channel "Hadronic" --tag "%s"' % (args.tag + "_impute", args.tag + "_impute_hut_FCNC_vs_SMHiggs_scaleGGH"))
##  #parallel_utils.run('python train_dnn.py --input "ttHHadronic_%s_FCNC_dnn_features_hut_FCNC_vs_SMHiggs_addMassConstraint.hdf5" --channel "Hadronic" --tag "%s" --preprocess_scheme "preprocess_scheme_Hadronic_hut_FCNC_vs_SMHiggs_addMassConstraint.json"' % (args.tag + "_impute", args.tag + "_impute_hut_FCNC_vs_SMHiggs_addMassConstraint"), False)
##  ##parallel_utils.run('python train_dnn.py --input "ttHHadronic_%s_FCNC_dnn_features_hut_FCNC_vs_SMHiggs_addExtMassConstraint.hdf5" --channel "Hadronic" --tag "%s" --preprocess_scheme "preprocess_scheme_Hadronic_hut_FCNC_vs_SMHiggs_addExtMassConstraint.json"' % (args.tag + "_impute", args.tag + "_impute_hut_FCNC_vs_SMHiggs_addExtMassConstraint"), False)
##
##  #parallel_utils.run('python train.py --input "ttHHadronic_%s_FCNC_features_hct_FCNC_vs_SMHiggs.hdf5" --channel "Hadronic" --tag "%s" --ext ""' % (args.tag + "_impute", args.tag + "_impute_hct_FCNC_vs_SMHiggs"))
##  #parallel_utils.run('python train_dnn.py --input "ttHHadronic_%s_FCNC_dnn_features_hct_FCNC_vs_SMHiggs_addMassConstraint.hdf5" --channel "Hadronic" --tag "%s" --preprocess_scheme "preprocess_scheme_Hadronic_hct_FCNC_vs_SMHiggs_addMassConstraint.json"' % (args.tag + "_impute", args.tag + "_impute_hct_FCNC_vs_SMHiggs_addMassConstraint"), False)
##  #parallel_utils.run('python train_dnn.py --input "ttHHadronic_%s_FCNC_dnn_features_hct_FCNC_vs_SMHiggs_addExtMassConstraint.hdf5" --channel "Hadronic" --tag "%s" --preprocess_scheme "preprocess_scheme_Hadronic_hct_FCNC_vs_SMHiggs_addExtMassConstraint.json"' % (args.tag + "_impute", args.tag + "_impute_hct_FCNC_vs_SMHiggs_addExtMassConstraint"), False)
##
##  # Zip DNN score into final fit ntuple
##  command_list = []
##  #command_list.append('python prep.py --input "../Loopers/MVABaby_ttHLeptonic_%s_FCNC.root" --channel "Leptonic" --fcnc_hut --tag "_hut_nonRes" --dnn_models "dnn_weights/metadata_Leptonic_%s_hut_FCNC_vs_SMHiggs.json" --dont_train_with_dnn --non_resonant_bkg' % (args.tag, args.tag))
##  #command_list.append('python prep.py --input "../Loopers/MVABaby_ttHLeptonic_%s_FCNC.root" --channel "Leptonic" --fcnc_hct --tag "_hct_nonRes" --dnn_models "dnn_weights/metadata_Leptonic_%s_hct_FCNC_vs_SMHiggs.json" --dont_train_with_dnn --non_resonant_bkg' % (args.tag, args.tag))
##  #command_list.append('python prep.py --input "../Loopers/MVABaby_ttHHadronic_%s_FCNC.root" --channel "Hadronic" --fcnc_hut --tag "_hut_nonRes" --dnn_models "dnn_weights/metadata_Hadronic_%s_hut_FCNC_vs_SMHiggs_addMassConstraint.json" --dont_train_with_dnn --non_resonant_bkg' % (args.tag + "_impute", args.tag + "_impute"))
##  #command_list.append('python prep.py --input "../Loopers/MVABaby_ttHHadronic_%s_FCNC.root" --channel "Hadronic" --fcnc_hct --tag "_hct_nonRes" --dnn_models "dnn_weights/metadata_Hadronic_%s_hct_FCNC_vs_SMHiggs_addMassConstraint.json" --dont_train_with_dnn --non_resonant_bkg' % (args.tag + "_impute", args.tag + "_impute"))
##
##
##  parallel_utils.submit_jobs(command_list, 8)
##
##  #parallel_utils.run('python train.py --input "ttHLeptonic_%s_FCNC_features_hut_nonRes.hdf5" --channel "Leptonic" --tag "%s" --ext "" --reference_mva "dnn_weights/metadata_Leptonic_%s_hut_FCNC_vs_SMHiggs.json" --reference_mva_name "fcnc_vs_smhiggs_dnn"' % (args.tag, args.tag + "_hut_FCNC_vs_NonRes", args.tag))
##  #parallel_utils.run('python train.py --input "ttHLeptonic_%s_FCNC_features_hct_nonRes.hdf5" --channel "Leptonic" --tag "%s" --ext "" --reference_mva "dnn_weights/metadata_Leptonic_%s_hct_FCNC_vs_SMHiggs.json" --reference_mva_name "fcnc_vs_smhiggs_dnn"' % (args.tag, args.tag + "_hct_FCNC_vs_NonRes", args.tag))
##  
##  #parallel_utils.run('python train.py --input "ttHHadronic_%s_FCNC_features_hut_nonRes.hdf5" --channel "Hadronic" --tag "%s" --ext "" --reference_mva "dnn_weights/metadata_Hadronic_%s_hut_FCNC_vs_SMHiggs_addMassConstraint.json" --reference_mva_name "fcnc_vs_smhiggs_dnn"' % (args.tag + "_impute", args.tag + "_impute_hut_FCNC_vs_NonRes", args.tag + "_impute"))
##  #parallel_utils.run('python train.py --input "ttHHadronic_%s_FCNC_features_hct_nonRes.hdf5" --channel "Hadronic" --tag "%s" --ext "" --reference_mva "dnn_weights/metadata_Hadronic_%s_hct_FCNC_vs_SMHiggs_addMassConstraint.json" --reference_mva_name "fcnc_vs_smhiggs_dnn"' % (args.tag + "_impute", args.tag + "_impute_hct_FCNC_vs_NonRes", args.tag + "_impute"))
##  # End FCNC vs SM Higgs BDT/DNN
##}}}
##}}}
