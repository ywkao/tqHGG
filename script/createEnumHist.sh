#!/bin/bash
set -e

# createEnumHist {{{
fin="include/list_histograms.h"
file="include/enumhist.h"
num=`cat ${fin} | grep -v "//" | wc -l`

echo "//### AUTOMATICALLY CREATED BY ./script/createEnumHist.sh ###//" > ${file}
echo "//### Original file: ${fin}" >> ${file}
echo "#ifndef __ENUMHIST_H__" >> ${file}
echo "#define __ENUMHIST_H__" >> ${file}
echo "#include <string>" >> ${file}
echo "using namespace std;" >> ${file}
echo "" >> ${file}


echo "enum histList {" >> ${file}
cat ${fin} | grep -v "//" | awk '{print $2}' | awk 'BEGIN{FS="*"}{print "    "$2","}' >> ${file}
echo "    totalHistNum" >> ${file}
echo "};" >> ${file}


echo "std::string histNames[totalHistNum]{" >> ${file}
cat ${fin} | grep -v "//" | awk '{print $5}' | awk 'BEGIN{FS="("}{print $2}' | awk 'BEGIN{FS=","}{if(NR=="'${num}'") print "    "$1; else print "    "$1","}' >> ${file}
echo "};" >> ${file}


echo "int histNbins[totalHistNum]{" >> ${file}
cat ${fin} | grep -v "//" | awk '{print $7 $2}' | awk 'BEGIN{FS=","}{if(NR=="'${num}'") print "    "$1"//"$2; else print "    "$1",//"$2}' | sed 's/\*//g' >> ${file}
echo "};" >> ${file}


echo "double histBinLow[totalHistNum]{" >> ${file}
cat ${fin} | grep -v "//" | awk '{print $8 $2}' | awk 'BEGIN{FS=","}{if(NR=="'${num}'") print "    "$1"//"$2; else print "    "$1",//"$2}' | sed 's/\*//g' >> ${file}
echo "};" >> ${file}


echo "double histBinHigh[totalHistNum]{" >> ${file}
cat ${fin} | grep -v "//" | awk '{print $9 $2}' | awk 'BEGIN{FS=");"}{if(NR=="'${num}'") print "    "$1"//"$2; else print "    "$1",//"$2}' | sed 's/\*//g' >> ${file}
echo "};" >> ${file}

echo "#endif" >> ${file}

echo "[INFO] the ${file} is updated!"
#}}}

dir_tmp="/tmp/ykao"
# automatic substitution stackHist.C {{{
time_stamp=`date +"%Y-%m-%d"`

stackFile=src/stackHist.C
targetFile=${dir_tmp}/tmp_stackHist.C
backupFile=${dir_tmp}/stackHist_${time_stamp}.C
#--------------------------------------------------
tmpFile_01="${dir_tmp}/tmp_stackHist_01"
num_tag=`grep -n "//### MakeStackHist" ${stackFile} | awk 'BEGIN{FS=":"}{print $1}'`
sed -n "1,${num_tag}p" ${stackFile} > ${tmpFile_01}
#--------------------------------------------------
tmpFile_02="${dir_tmp}/tmp_stackHist_02"
enumFile=include/enumhist.h
num_start=`grep -n "enum histList" ${enumFile} | awk 'BEGIN{FS=":"}{print $1}'`
num_end=`grep -n "totalHistNum" ${enumFile} | awk 'BEGIN{FS=":"}{print $1}' | head -n1`
num_start=$((num_start+1))
num_end=$((num_end-1))
sed -n "${num_start},${num_end}p" ${enumFile} | awk '{print $1}' | awk 'BEGIN{FS=","}{print "MakeStackHist(\""$1"\");"}' > ${tmpFile_02}
sed -i 's/Make/    Make/g' ${tmpFile_02}
#--------------------------------------------------
tmpFile_03="${dir_tmp}/tmp_stackHist_03"
num_tag=`grep -n "//end of MakeStackHist" ${stackFile} | awk 'BEGIN{FS=":"}{print $1}'`
sed -n "${num_tag},$ p" ${stackFile} > ${tmpFile_03}
#--------------------------------------------------
cat ${tmpFile_01} >  ${targetFile}
cat ${tmpFile_02} >> ${targetFile}
cat ${tmpFile_03} >> ${targetFile}
#--------------------------------------------------
cp -p ${stackFile} ${backupFile}
cp -p ${targetFile} ${stackFile}
echo "[INFO] the ${stackFile} is updated!"

rm ${dir_tmp}/tmp_*
#}}}
