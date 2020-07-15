#!/bin/bash
set -e

for file in \
log_chi_had_false_false_true_true \
log_chi_had_false_false_true_false \
log_chi_had_false_true_false_true \
log_chi_had_false_true_false_false \
log_chi_had_true_false_false_true \
log_chi_had_true_false_false_false

do
    #--- b-jet ---#
    #echo $file
    #cat $file | grep -i info-hist
    #cat $file | grep -i percentage
    #echo ""
    #echo ""

    echo $file
    cat $file | grep Correctly | grep sig
    #echo "--------------------"
    #cat $file | grep Correctly | grep tbw 
    #echo "--------------------"
    #cat $file | grep Correctly | grep tqh 
    echo ""
    echo ""
done

