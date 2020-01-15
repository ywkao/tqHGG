#!/bin/bash
set -e

fin="include/histogramList.h"
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
