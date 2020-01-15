#!/bin/bash
set -e
file=src/generalChiSquareStudy.cpp

function EditFile(){
    
    loose_condition="bool bool_bjet_is_loose  = $1;"
    medium_condition="bool bool_bjet_is_medium = $2;"
    tight_condition="bool bool_bjet_is_tight  = $3;"
    num_condition="bool bool_num_bjets_is_exactly_one = $4;"
    
    sed -i "43c ${loose_condition}" $file
    sed -i "44c ${medium_condition}" $file
    sed -i "45c ${tight_condition}" $file
    sed -i "46c ${num_condition}" $file

    ./execution.sh | tee log_chi_had_$1_$2_$3_$4
}

EditFile false false true true
EditFile false true false true
EditFile true false false true
EditFile false false true false
EditFile false true false false
EditFile true false false false

#bool bool_bjet_is_loose  = true;
#bool bool_bjet_is_medium = true;
#bool bool_bjet_is_tight  = false;
#bool bool_num_bjets_is_exactly_one = true;
