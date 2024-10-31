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

    P5_FILEPATH=$PIPELINE_INPUT # first argument
    P5_DIR=${P5_FILEPATH%/*} # strip filename from .pod5 filepath
    PARENT_DIR=${P5_DIR%/*} # get folder one heirarchy higher
    # name of the .pod5 file (strip the path and get only the name with extension)
    P5_FILENAME=$(basename $P5_FILEPATH)
    # name of the .pod5 file without the extension
    P5_PREFIX=${P5_FILENAME%.*}

    # deduce the directory for slow5 and sam files
    if [[ "$P5_DIR" =~ .*"pod5_pass".* ]]; then
        SLOW5_DIR=$(echo $P5_DIR | sed 's/pod5_pass/slow5_pass/g')
        SAM_DIR=$(echo $P5_DIR | sed 's/pod5_pass/sam_pass/g')
        LOG_DIR=$(echo $P5_DIR | sed 's/pod5_pass/pass_logs/g')
    elif [[ "$P5_DIR" =~ .*"pod5_fail".* ]]; then
        SLOW5_DIR=$(echo $P5_DIR | sed 's/pod5_fail/slow5_fail/g')
        SAM_DIR=$(echo $P5_DIR | sed 's/pod5_fail/sam_fail/g')
        LOG_DIR=$(echo $P5_DIR | sed 's/pod5_fail/fail_logs/g')
    elif [[ "$P5_DIR" =~ .*"pod5_skip".* ]]; then
        SLOW5_DIR=$(echo $P5_DIR | sed 's/pod5_skip/slow5_skip/g')
        SAM_DIR=$(echo $P5_DIR | sed 's/pod5_skip/sam_skip/g')
        LOG_DIR=$(echo $P5_DIR | sed 's/pod5_skip/skip_logs/g')
    else
        SLOW5_DIR=$PARENT_DIR/slow5/
        SAM_DIR=$PARENT_DIR/sam/
        LOG_DIR=$PARENT_DIR/logs/
    fi
    if [ -z "$SLOW5_DIR" ] || [ -z "$LOG_DIR" ] ; then
        SLOW5_DIR=$PARENT_DIR/slow5/
        SAM_DIR=$PARENT_DIR/sam/
        LOG_DIR=$PARENT_DIR/logs/
    fi

    test -d $SLOW5_DIR/ || { mkdir -p $SLOW5_DIR/; echo "[pipeline.sh] Created $SLOW5_DIR/. Converted SLOW5 files will be here."; }
    test -d $LOG_DIR/ || { mkdir -p $LOG_DIR/; echo "[pipeline.sh] Created $LOG_DIR/. individual logs for each conversion will be here."; }
    test -d $SAM_DIR/ || { mkdir -p $SAM_DIR/; echo "[pipeline.sh] Created $SAM_DIR/. SAM files will be here."; }

    SLOW5_FILEPATH=$SLOW5_DIR/$P5_PREFIX.blow5
    UNALN_SAM_FILEPATH=$SAM_DIR/$P5_PREFIX.remora.unaln.sam
    FASTQ_FILEPATH=$SAM_DIR/$P5_PREFIX.fastq
    UNSORTED_BAM_FILEPATH=$SAM_DIR/$P5_PREFIX.remora.unsorted.bam
    BAM_FILEPATH=$SAM_DIR/$P5_PREFIX.remora.bam
    LOG_FILEPATH=$LOG_DIR/$P5_PREFIX.log

    t1=$(date)
    ${BLUECRAB} p2s -p1 $P5_FILEPATH -o $SLOW5_FILEPATH >> $LOG_FILEPATH 2>&1 || die "blue-crab failed. check $LOG_FILEPATH"
    t2=$(date)
    echo -e "$P5_FILEPATH\tp2s\t${t1}\t${t2}" >> ${TIME_LOG}

    ${EEL} --log $PARENT_DIR/buttery_eel_logs --call_mods --config $DORADO_MODEL -i $SLOW5_FILEPATH -o $UNALN_SAM_FILEPATH >> $LOG_FILEPATH 2>&1 || die "buttery-eel failed. check $LOG_FILEPATH"
    t3=$(date)
    echo -e "$P5_FILEPATH\teel\t${t2}\t${t3}" >> ${TIME_LOG}

    ${SAMTOOLS} fastq -@ 8 -TMM,ML $UNALN_SAM_FILEPATH > $FASTQ_FILEPATH 2>> $LOG_FILEPATH || die "samtools fastq failed. check $LOG_FILEPATH"
    t4=$(date)
    echo -e "$P5_FILEPATH\tsam-fastq\t${t3}\t${t4}" >> ${TIME_LOG}

    ${MINIMAP2} -t 8 -ax map-ont --sam-hit-only -Y -y --secondary=no $REFIDX $FASTQ_FILEPATH > $UNSORTED_BAM_FILEPATH 2>> $LOG_FILEPATH || die "minimap2 failed. check $LOG_FILEPATH"
    t5=$(date)
    echo -e "$P5_FILEPATH\tminimap2\t${t4}\t${t5}" >> ${TIME_LOG}

    ${SAMTOOLS} sort -@ 8 -o $BAM_FILEPATH $UNSORTED_BAM_FILEPATH >> $LOG_FILEPATH 2>&1 || die "samtools sort failed. check $LOG_FILEPATH"
    t6=$(date)
    echo -e "$P5_FILEPATH\tsam-sort\t${t5}\t${t6}" >> ${TIME_LOG}

    ${SAMTOOLS} index -@ 8 $BAM_FILEPATH >> $LOG_FILEPATH 2>&1 || die "samtools index failed. check $LOG_FILEPATH"
    t7=$(date)
    echo -e "$P5_FILEPATH\tsam-index\t${t6}\t${t7}" >> ${TIME_LOG}

    ${F5C} index --slow5 $SLOW5_FILEPATH $FASTQ_FILEPATH >> $LOG_FILEPATH 2>&1 || die $RED"Running f5c index on $SLOW5_FILEPATH failed. Please check log at $LOG_FILEPATH"$NORMAL
    t8=$(date)
    echo -e "$P5_FILEPATH\tf5c-index\t${t7}\t${t8}" >> ${TIME_LOG}

    ${F5C} call-methylation --slow5 $SLOW5_FILEPATH -b $BAM_FILEPATH -g $REF -r $FASTQ_FILEPATH > $F5C_FILEPATH 2>> $LOG_FILEPATH || die $RED"Running f5c call-methylation on $SLOW5_FILEPATH failed. Please check log at $LOG_FILEPATH"$NORMAL
    t9=$(date)
    echo -e "$P5_FILEPATH\tf5c\t${t8}\t${t9}" >> ${TIME_LOG}

    PIPELINE_OUTPUT="$F5C_FILEPATH"
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
