#!/bin/bash

# Default values
TIME_LOG="realfreq_pipeline_time.log"
TMP_FAILED="realfreq_failed_list.log"
MAX_PROC=1
CHECK_TOOLS=false

# Colour codes
YELLOW="\e[33m"
RED="\e[31m"
NORMAL="\033[0;39m"

# terminate script
die() { echo -e $RED"$1"$NORMAL >&2; echo; exit 1; }

## Handle flags
while getopts "l:f:p:c" o; do
    case "${o}" in
        l)
            TIME_LOG=${OPTARG}
            ;;
        f)
            TMP_FAILED=${OPTARG}
            ;;
        p)
            MAX_PROC=${OPTARG}
            ;;
        c)
            CHECK_TOOLS=true
            ;;
        *)
            echo "[pipeline.sh] Incorrect args"
            exit 1
            ;;
    esac
done
shift $((OPTIND-1))

#====================User Start - Tool check====================

[ -z ${SLOW5TOOLS} ] && export SLOW5TOOLS=slow5tools
[ -z ${BLUECRAB} ] && export BLUECRAB=blue-crab
[ -z ${BUTTERY_EEL} ] && export BUTTERY_EEL=buttery-eel
[ -z ${MINIMAP2} ] && export MINIMAP2=minimap2
[ -z ${SAMTOOLS} ] && export SAMTOOLS=samtools
[ -z ${REALFREQ} ] && export REALFREQ=realfreq

$SLOW5TOOLS --version &> /dev/null || die "[pipeline.sh] slow5tools not found. Add to PATH or export SLOW5TOOLS=/path/to/slow5tools. Exiting."
$BLUECRAB --version &> /dev/null || die "[pipeline.sh] bluecrab not found. Add to PATH or export BLUECRAB=/path/to/bluecrab. Exiting."
command -v $BUTTERY_EEL &> /dev/null || die "[pipeline.sh] buttery-eel not found. Add to PATH or export BUTTERY_EEL=/path/to/buttery-eel. Exiting."
$MINIMAP2 --version &> /dev/null || die "[pipeline.sh] minimap2 not found. Add to PATH or export MINIMAP2=/path/to/minimap2. Exiting."
$SAMTOOLS --version &> /dev/null || die "[pipeline.sh] samtools not found. Add to PATH or export SAMTOOLS=/path/to/samtools. Exiting."

EEL="$BUTTERY_EEL -g $DORADO_BIN --port 5000 --use_tcp --device cuda:all"

[ -z ${DORADO_BIN} ] && die "[pipeline.sh] DORADO_BIN not set. export DORADO_BIN=/path/to/dorado. Exiting."
[ -z ${DORADO_MODEL} ] && die "[pipeline.sh] DORADO_MODEL not set. export DORADO_MODEL=/path/to/dorado_model. Exiting."
[ -z ${REF} ] && die "[pipeline.sh] REF not set. export REF=/path/to/ref.fa. Exiting."
[ -z ${REFIDX} ] && die "[pipeline.sh] REFIDX not set. export REFIDX=/path/to/ref.idx. Exiting."

#====================User End - Tool check====================

if [ "$CHECK_TOOLS" = true ]; then
    echo "[pipeline.sh] Tool check successful"
    exit 0
fi

echo "[pipeline.sh] Starting pipeline with $MAX_PROC max processes"
counter=0
while read PIPELINE_INPUT
do
(
    #=======================User Start - Pipeline=======================

    MODBAM_FILEPATH=$PIPELINE_INPUT # first argument
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

    t1=$(date)
    ${SAMTOOLS} index -@ 8 $MODBAM_FILEPATH >> $LOG_FILEPATH 2>&1 || die $RED"Indexing $MODBAM_FILEPATH failed. Please check log at $LOG_FILEPATH"$NORMAL
    t2=$(date)
    echo -e "$MODBAM_FILEPATH\tsam-index\t${t1}\t${t2}" >> ${TIME_LOG}

    PIPELINE_OUTPUT="$MODBAM_FILEPATH"
    #=======================User End - Pipeline=======================
    
    echo "realfreq-pipeline-output:$PIPELINE_OUTPUT"
)&
    ((counter++))
    if [ $counter -ge $MAX_PROC ]; then
        echo "[pipeline.sh] Waiting for $counter jobs to finish."
        wait
        counter=0
    fi
done
wait
