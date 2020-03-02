#!/bin/bash
set -e
tag=$1
file=src/mva_TMVAClassification_leptonic.C 
tmpDir=src/tmp/task_${tag}
# function{{{
function Edit(){
    line=$1; pre=$2; content=$3;
    sed -i -e "${line}c ${content}" $file
    sed -i -e ''"$line"'s/^'"$pre"'/   '"$pre"'/g' $file
}
#}}}
Edit 66 "const" "const char TAG[32] = \"${tag}\";"
#// Cut optimisation{{{
Edit 81 "Use" "Use[\"Cuts\"]            = 0;"
Edit 82 "Use" "Use[\"CutsD\"]           = 0;"
Edit 83 "Use" "Use[\"CutsPCA\"]         = 0;"
Edit 84 "Use" "Use[\"CutsGA\"]          = 0; // [ERROR]"
Edit 85 "Use" "Use[\"CutsSA\"]          = 0;"
#}}}
#// 1-dimensional likelihood ({{{
Edit 88 "Use" "Use[\"Likelihood\"]      = 0;"
Edit 89 "Use" "Use[\"LikelihoodD\"]     = 0; // the \"D\" extension indicates decorrelated input variables (see option strings)"
Edit 90 "Use" "Use[\"LikelihoodPCA\"]   = 0; // the \"PCA\" extension indicates PCA-transformed input variables (see option strings)"
Edit 91 "Use" "Use[\"LikelihoodKDE\"]   = 0; // [ERROR]"
Edit 92 "Use" "Use[\"LikelihoodMIX\"]   = 0;"
#}}}
#// Mutidimensional likelihood and Nearest-Neighbour methods{{{
Edit 95 "Use" "Use[\"PDERS\"]           = 0;"
Edit 96 "Use" "Use[\"PDERSD\"]          = 0;"
Edit 97 "Use" "Use[\"PDERSPCA\"]        = 0;"
Edit 98 "Use" "Use[\"PDEFoam\"]         = 0;"
Edit 99 "Use" "Use[\"PDEFoamBoost\"]    = 0; // uses generalised MVA method boosting"
Edit 100 "Use" "Use[\"KNN\"]             = 0; // k-nearest neighbour method"
#}}}
#// Linear Discriminant Analysis{{{
Edit 103 "Use" "Use[\"LD\"]              = 0; // Linear Discriminant identical to Fisher"
Edit 104 "Use" "Use[\"Fisher\"]          = 0;"
Edit 105 "Use" "Use[\"FisherG\"]         = 0;"
Edit 106 "Use" "Use[\"BoostedFisher\"]   = 0; // uses generalised MVA method boosting"
Edit 107 "Use" "Use[\"HMatrix\"]         = 0;"
#}}}
#// Function Discriminant analysis{{{
Edit 110 "Use" "Use[\"FDA_GA\"]          = 0; // minimisation of user-defined function using Genetics Algorithm"
Edit 111 "Use" "Use[\"FDA_SA\"]          = 0;"
Edit 112 "Use" "Use[\"FDA_MC\"]          = 0;"
Edit 113 "Use" "Use[\"FDA_MT\"]          = 0;"
Edit 114 "Use" "Use[\"FDA_GAMT\"]        = 0;"
Edit 115 "Use" "Use[\"FDA_MCMT\"]        = 0;"
#}}}
#// Neural Networks (all are feed-forward Multilayer Perceptrons){{{
Edit 118 "Use" "Use[\"MLP\"]             = 0; // Recommended ANN"
Edit 119 "Use" "Use[\"MLPBFGS\"]         = 0; // Recommended ANN with optional training method"
Edit 120 "Use" "Use[\"MLPBNN\"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator"
Edit 121 "Use" "Use[\"CFMlpANN\"]        = 0; // Depreciated ANN from ALEPH"
Edit 122 "Use" "Use[\"TMlpANN\"]         = 0; // ROOT's own ANN"
Edit 124 "Use" "Use[\"DNN_GPU\"]         = 0; // CUDA-accelerated DNN training."
Edit 126 "Use" "Use[\"DNN_GPU\"]         = 0;"
Edit 130 "Use" "Use[\"DNN_CPU\"]         = 0; // Multi-core accelerated DNN."
Edit 132 "Use" "Use[\"DNN_CPU\"]         = 0;"
#}}}
#// Support Vector Machine{{{
Edit 136 "Use" "Use[\"SVM\"]             = 0;"
Edit 137 "Use" "Use[\"BDT\"]             = 1; // uses Adaptive Boost"
Edit 138 "Use" "Use[\"BDTG\"]            = 1; // uses Gradient Boost"
Edit 139 "Use" "Use[\"BDTB\"]            = 0; // uses Bagging"
Edit 140 "Use" "Use[\"BDTD\"]            = 0; // decorrelation + Adaptive Boost"
Edit 141 "Use" "Use[\"BDTF\"]            = 0; // allow usage of fisher discriminant for node splitting"
#}}}
#// Friedman's RuleFit method, ie, an optimised series of cuts ({{{
Edit 144 "Use" "Use[\"RuleFit\"]         = 0;"
#}}}
#mkdir -p ${tmpDir}
#cp -p ${file} ${tmpDir}/TMVAClassification.C
