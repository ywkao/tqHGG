#!/bin/bash
set -e

FILE=include/main.h
function prepare(){
    NAME=$1
    BIN=bin
    #BIN=bin_2
    #BIN=bin_3

    if [ -f ./${BIN}/${NAME} ]; then
        echo "[INFO] rm ./${BIN}/${NAME} ./build/${NAME}.o"; rm ./${BIN}/${NAME} ./build/${NAME}.o
    else
        echo "[INFO] No executable ./${BIN}/${NAME}."
    fi

    # executable with int NPu (default setting in the code)
    make
    echo "[INFO] mv ./${BIN}/${NAME} ./${BIN}/${NAME}_npu_int"; mv ./${BIN}/${NAME} ./${BIN}/${NAME}_npu_int
    echo "[INFO] mv ./build/${NAME}.o ./build/${NAME}_npu_int.o"; mv ./build/${NAME}.o ./build/${NAME}_npu_int.o

    # executable with float NPu
    sed -i 's/Int_t EvtInfo_NPu;/\/\/Int_t EvtInfo_NPu;/' ${FILE}
    sed -i 's/\/\/float EvtInfo_NPu;/float EvtInfo_NPu;/' ${FILE}
    if [ ${NAME} == "preselection" ] || [ ${NAME} == "preselection_npustudy" ]; then
        sed -i 's/EvtInfo_NPu\/I/EvtInfo_NPu\/F/' src/${NAME}_exe.cpp
    fi
    make
    echo "[INFO] mv ./${BIN}/${NAME} ./${BIN}/${NAME}_npu_float"; mv ./${BIN}/${NAME} ./${BIN}/${NAME}_npu_float
    echo "[INFO] rm ./build/${NAME}.o"; rm ./build/${NAME}.o

    # reset!
    if [ ${NAME} == "preselection" ] || [ ${NAME} == "preselection_npustudy" ]; then
        sed -i 's/EvtInfo_NPu\/F/EvtInfo_NPu\/I/' src/${NAME}_exe.cpp
    fi
    sed -i 's/\/\/Int_t EvtInfo_NPu;/Int_t EvtInfo_NPu;/' ${FILE}
    sed -i 's/float EvtInfo_NPu;/\/\/float EvtInfo_NPu;/' ${FILE}
    echo "[INFO] mv ./${BIN}/${NAME}_npu_int ./${BIN}/${NAME}"; mv ./${BIN}/${NAME}_npu_int ./${BIN}/${NAME}
    echo "[INFO] mv ./build/${NAME}_npu_int.o ./build/${NAME}.o"; mv ./build/${NAME}_npu_int.o ./build/${NAME}.o
}

stage=$1
prepare ${stage}


#prepare "preselection_npustudy"
#prepare "preselection"
#prepare "selection"


#echo "rm ./bin/selection ./build/selection.o"
#rm ./bin/selection ./build/selection.o
#make
#
#echo "mv ./bin/selection ./bin/selection_npu_float"
#mv ./bin/selection ./bin/selection_npu_float
#echo "rm ./build/selection.o"
#rm ./build/selection.o
