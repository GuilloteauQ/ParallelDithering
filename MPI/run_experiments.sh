#! /bin/bash

FILE_TO_WRITE=$1
IMAGE=$2
ITERATIONS=$3
MAX_PROCESSES=$4

for ((i = 1; i <= ${ITERATIONS}; i++))
do
    echo ">> ${i}/${ITERATIONS}"
    for ((p = 1; p <= ${MAX_PROCESSES}; p++))
    do
        mpirun -np ${p} ./dithering_bw ${IMAGE} >> ${FILE_TO_WRITE}
    done
done

