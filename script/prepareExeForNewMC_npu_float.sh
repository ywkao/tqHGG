#!/bin/bash
set -e

FILE=include/main.h
function prepare(){
    NAME=$1

    if [ -f ./bin/${NAME} ]; then
        echo "[INFO] rm ./bin/${NAME} ./build/${NAME}.o"; rm ./bin/${NAME} ./build/${NAME}.o
    else
        echo "[INFO] No executable ./bin/${NAME}."
    fi

    # executable with int NPu (default setting in the code)
    make
    echo "[INFO] mv ./bin/${NAME} ./bin/${NAME}_npu_int"; mv ./bin/${NAME} ./bin/${NAME}_npu_int
    echo "[INFO] mv ./build/${NAME}.o ./build/${NAME}_npu_int.o"; mv ./build/${NAME}.o ./build/${NAME}_npu_int.o

    # executable with float NPu
    sed -i 's/Int_t EvtInfo_NPu;/\/\/Int_t EvtInfo_NPu;/' ${FILE}
    sed -i 's/\/\/float EvtInfo_NPu;/float EvtInfo_NPu;/' ${FILE}
    if [ ${NAME} == "preselection" ] || [ ${NAME} == "preselection_npustudy" ]; then
        sed -i 's/EvtInfo_NPu\/I/EvtInfo_NPu\/F/' src/${NAME}.cpp
    fi
    make
    echo "[INFO] mv ./bin/${NAME} ./bin/${NAME}_npu_float"; mv ./bin/${NAME} ./bin/${NAME}_npu_float
    echo "[INFO] rm ./build/${NAME}.o"; rm ./build/${NAME}.o

    # reset!
    if [ ${NAME} == "preselection" ] || [ ${NAME} == "preselection_npustudy" ]; then
        sed -i 's/EvtInfo_NPu\/F/EvtInfo_NPu\/I/' src/${NAME}.cpp
    fi
    sed -i 's/\/\/Int_t EvtInfo_NPu;/Int_t EvtInfo_NPu;/' ${FILE}
    sed -i 's/float EvtInfo_NPu;/\/\/float EvtInfo_NPu;/' ${FILE}
    echo "[INFO] mv ./bin/${NAME}_npu_int ./bin/${NAME}"; mv ./bin/${NAME}_npu_int ./bin/${NAME}
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
