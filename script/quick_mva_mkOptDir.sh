#!/bin/bash
tag=$1
source="src/mva_TMVAClassification_leptonic.C"
tmpdir="src/tmp/${tag}"
target="${tmpdir}/TMVAClassification.C"

mkdir -p ${tmpdir}
cp ${source} ${target}

echo "ls ${target}"; ls ${target};
echo "diff ${source} ${target}"; diff ${source} ${target};
