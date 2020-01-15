#! /bin/bash

IMAGE=Images/mandrill.pgm
FACTOR=1
WITH_HOSTFILE=0

while [[ $# -gt 0 ]]
do
        key="$1"

        case $key in
            -p)
                PROCESSES=$2
                shift
                shift
                ;;
            -k)
                K=$2
                shift
                shift
                ;;
            -r)
                R=$2
                shift
                shift
                ;;
            -f)
                FACTOR=$2
                shift
                shift
                ;;
	    --hostfile)
		WITH_HOSTFILE=1
		HOSTFILE_PATH=$2
		shift
		shift
		;;
            --image)
                IMAGE=$2
                shift
                shift
                ;;
            *)
                echo "error while parsing the arguments"
                exit
                ;;
        esac
    done

if [[ $WITH_HOSTFILE -ne 1 ]]
then
	mpirun --allow-run-as-root -np ${PROCESSES} ./dithering_bw ${IMAGE} ${K} ${R} ${FACTOR}
else
	mpirun --allow-run-as-root -np ${PROCESSES} --hostfile ${HOSTFILE_PATH} ./dithering_bw ${IMAGE} ${K} ${R} ${FACTOR}
fi



