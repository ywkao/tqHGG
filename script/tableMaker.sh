#!/bin/bash
# vim: set fdm=marker:
set -e

tag="$1"
#tag="_update"
#tag="_latest"

targetFile=tables/tables.tex
write_method=0 # for initiation
time_stamp=`date +"%Y-%m-%d"`

# functions{{{
function write(){
    message=$1;
    if [[ $write_method -eq 0 ]];
    then
        echo $message > $targetFile
    else
        echo $message >> $targetFile
    fi
}

function make_table(){
    log=$1
    channel=$2
    write ""
    write "{\renewcommand{\arraystretch}{1.0}"
    write "\begin{center}"
    write "    \captionof{table}{${time_stamp} ${channel}}"
    write "\begin{tabular}{lrrl}"
    write "    \hline\hline"
    write "    Processes & Entries &\multicolumn{2}{c}{Yields}\\\\"
    write "    \hline\hline"
    cat ${log} | grep "&" | grep -v -E "Processes|MC|Data" | grep -v -E "ggH|VH|VBF|ttH" | grep -v "res" >> ${targetFile}
    write "    \hline"
    cat ${log} | grep "&" | grep -E "MC|Data" >> ${targetFile}
    write "    \hline\hline\\\\"
    write "\end{tabular}"
    write "\end{center}"
    write "}"
}
#}}}

write "\clearpage"; write_method=1
write "\chapter{Table of yields}"
make_table plots_hadronic${tag}/161718/log/info_stack_plots_hadronic.txt "hadronic"
write "\clearpage"
make_table plots_leptonic${tag}/161718/log/info_stack_plots_leptonic.txt "leptonic"

sed -i 's/_/\\_/g' ${targetFile}

cd tables; ./Compile.sh
echo "mv tableMaker.pdf tableMaker${tag}.pdf"; mv tableMaker.pdf tableMaker${tag}.pdf;
