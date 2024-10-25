#!/bin/bash

TMP_FILE="attempted_list.log"
TMP_FAILED="failed_list.log"
GUPPY_BIN=""
REF=""
REFIDX=""
MODEL=""

# terminate script
die() {
	echo -e "$1" >&2
	echo
	exit 1
}

LOG=start_end_trace.log

# MAX_PROC=$(nproc)
# MAX_PROC=$(echo "${MAX_PROC}/2" | bc)
MAX_PROC=1

# Colour codes for printing
YELLOW="\e[33m"
RED="\e[31m"
NORMAL="\033[0;39m"

## Handle flags
while getopts "d:l:f:p:g:r:i:m:" o; do
    case "${o}" in
        d)
            TMP_FILE=${OPTARG}
            ;;
		l)
            LOG=${OPTARG}
			;;
        f)
            TMP_FAILED=${OPTARG}
            ;;
        p)
            MAX_PROC=${OPTARG}
            ;;
        g)
            GUPPY_BIN=${OPTARG}
            ;;
        r) 
            REF=${OPTARG}
            ;;
        i)  
            REFIDX=${OPTARG}
            ;;
        m)  
            MODEL=${OPTARG}
            ;;
        *)
            echo "[pipeline.sh] Incorrect args"
            exit 1
            ;;
    esac
done
shift $((OPTIND-1))

$SAMTOOLS --version &> /dev/null || die $RED"[pipeline.sh] samtools not found in path. Exiting."$NORMAL

echo "[pipeline.sh] Starting pipeline with $MAX_PROC max processes"
#test -e ${LOG}  && rm ${LOG}
counter=0
while read FILE
do
(
    MODBAM_FILEPATH=$FILE # first argument
    MODBAM_DIR=${MODBAM_FILEPATH%/*} # strip filename from .pod5 filepath
    PARENT_DIR=${MODBAM_DIR%/*} # get folder one heirarchy higher
    # name of the .bam file (strip the path and get only the name with extension)
    MODBAM_FILENAME=$(basename $MODBAM_FILEPATH)
    # name of the .bam file without the extension
    MODBAM_PREFIX=${MODBAM_FILENAME%.*}

    MODBAM_INDEX_FILEPATH=$MODBAM_DIR/$MODBAM_PREFIX.bai

    LOG_DIR=$PARENT_DIR/logs/
    test -d $LOG_DIR/ || { mkdir -p $LOG_DIR/; echo "[pipeline.sh] Created $LOG_DIR/. individual logs for each conversion will be here."; }

    LOG_FILEPATH=$LOG_DIR/$MODBAM_PREFIX.log

    START_TIME=$(date)
    echo "[pipeline.sh::${t6}]  Indexing $MODBAM_FILEPATH"
    /usr/bin/time -v ${SAMTOOLS} index -@ 8 $MODBAM_FILEPATH >> $LOG_FILEPATH 2>&1 || die $RED"Indexing $MODBAM_FILEPATH failed. Please check log at $LOG_FILEPATH"$NORMAL
    t7=$(date)
    echo "[pipeline.sh::${t7}]  Finished indexing $MODBAM_FILEPATH"
    echo -e "$MODBAM_FILEPATH\tsam-index\t${t6}\t${t7}" >> ${LOG}
    
    END_TIME=$(date)

    echo "[pipeline.sh::${END_TIME}]  Finished pipeline for $MODBAM_FILEPATH output: $MODBAM_FILEPATH"
   
    echo -e "${MODBAM_FILEPATH}\ttotal\t${START_TIME}\t${END_TIME}" >> ${LOG}
)&
    ((counter++))
    if [ $counter -ge $MAX_PROC ]; then
        echo "[pipeline.sh] Waiting for $counter jobs to finish."
        wait
        counter=0
    fi
done
wait
