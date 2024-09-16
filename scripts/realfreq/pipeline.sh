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

$SLOW5TOOLS --version &> /dev/null || die $RED"[pipeline.sh] slow5tools not found in path. Exiting."$NORMAL
$BLUECRAB --version &> /dev/null || die $RED"[pipeline.sh] bluecrab not found in path. Exiting."$NORMAL
command -v $BUTTERY_EEL &> /dev/null || die $RED"[pipeline.sh] buttery-eel not found in path. Exiting."$NORMAL
$MINIMAP2 --version &> /dev/null || die $RED"[pipeline.sh] minimap2 not found in path. Exiting."$NORMAL
$SAMTOOLS --version &> /dev/null || die $RED"[pipeline.sh] samtools not found in path. Exiting."$NORMAL

EEL="$BUTTERY_EEL -g $GUPPY_BIN --port 5000 --use_tcp --device cuda:all"

echo "[pipeline.sh] Starting pipeline with $MAX_PROC max processes"
#test -e ${LOG}  && rm ${LOG}
counter=0
while read FILE
do
(
    P5_FILEPATH=$FILE # first argument
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
        SLOW5_LOG_DIR=$(echo $P5_DIR | sed 's/pod5_pass/slow5_pass_logs/g')
        SAM_LOG_DIR=$(echo $P5_DIR | sed 's/pod5_pass/sam_pass_logs/g')
    elif [[ "$P5_DIR" =~ .*"pod5_fail".* ]]; then
        SLOW5_DIR=$(echo $P5_DIR | sed 's/pod5_fail/slow5_fail/g')
        SAM_DIR=$(echo $P5_DIR | sed 's/pod5_fail/sam_fail/g')
        SLOW5_LOG_DIR=$(echo $P5_DIR | sed 's/pod5_fail/slow5_fail_logs/g')
        SAM_LOG_DIR=$(echo $P5_DIR | sed 's/pod5_fail/sam_fail_logs/g')
    elif [[ "$P5_DIR" =~ .*"pod5_skip".* ]]; then
        SLOW5_DIR=$(echo $P5_DIR | sed 's/pod5_skip/slow5_skip/g')
        SAM_DIR=$(echo $P5_DIR | sed 's/pod5_skip/sam_skip/g')
        SLOW5_LOG_DIR=$(echo $P5_DIR | sed 's/pod5_skip/slow5_skip_logs/g')
        SAM_LOG_DIR=$(echo $P5_DIR | sed 's/pod5_skip/sam_skip_logs/g')
    else
        SLOW5_DIR=$PARENT_DIR/slow5/
        SAM_DIR=$PARENT_DIR/sam/
        SLOW5_LOG_DIR=$PARENT_DIR/slow5_logs/
        SAM_LOG_DIR=$PARENT_DIR/sam_logs/
    fi
    if [ -z "$SLOW5_DIR" ] || [ -z "$SLOW5_LOG_DIR" ] ; then
        SLOW5_DIR=$PARENT_DIR/slow5/
        SAM_DIR=$PARENT_DIR/sam/
        SLOW5_LOG_DIR=$PARENT_DIR/slow5_logs/
        SAM_LOG_DIR=$PARENT_DIR/sam_logs/
    fi

    test -d $SLOW5_DIR/ || { mkdir -p $SLOW5_DIR/; echo "[pipeline.sh] Created $SLOW5_DIR/. Converted SLOW5 files will be here."; }
    test -d $SLOW5_LOG_DIR/ || { mkdir -p $SLOW5_LOG_DIR/; echo "[pipeline.sh] Created $SLOW5_LOG_DIR/. SLOW5 individual logs for each conversion will be here."; }
    test -d $SAM_DIR/ || { mkdir -p $SAM_DIR/; echo "[pipeline.sh] Created $SAM_DIR/. SAM files will be here."; }
    test -d $SAM_LOG_DIR/ || { mkdir -p $SAM_LOG_DIR/; echo "[pipeline.sh] Created $SAM_LOG_DIR/. SAM individual logs for each conversion will be here."; }

    SLOW5_FILEPATH=$SLOW5_DIR/$P5_PREFIX.blow5
    UNALN_SAM_FILEPATH=$SAM_DIR/$P5_PREFIX.remora.unaln.sam
    FASTQ_FILEPATH=$P5_DIR/$P5_PREFIX.fastq
    UNSORTED_BAM_FILEPATH=$SAM_DIR/$P5_PREFIX.remora.unsorted.bam
    BAM_FILEPATH=$SAM_DIR/$P5_PREFIX.remora.bam
    LOG_FILEPATH=$SLOW5_LOG_DIR/$P5_PREFIX.log

    test -e $SLOW5_FILEPATH &&  { echo -e $RED"$SLOW5_FILEPATH already exists. Converting $FILE to $SLOW5_FILEPATH failed."$NORMAL; echo $FILE >> $TMP_FAILED; }

    START_TIME=$(date)
    echo "[pipeline.sh::${START_TIME}]  Converting $FILE to $SLOW5_FILEPATH"
    ${BLUECRAB} p2s -p1 $FILE -o $SLOW5_FILEPATH 2> $LOG_FILEPATH || die $RED"Converting $FILE to $SLOW5_FILEPATH failed. Please check log at $LOG_FILEPATH"$NORMAL
    t2=$(date)
    echo "[pipeline.sh::${t2}]  Finished converting $FILE to $SLOW5_FILEPATH"

    echo "[pipeline.sh::${t2}]  Running buttery-eel on $SLOW5_FILEPATH"
    ${EEL} --call_mods --config $MODEL -i $SLOW5_FILEPATH -o $UNALN_SAM_FILEPATH 2> $LOG_FILEPATH || die $RED"Running buttery-eel on $SLOW5_FILEPATH failed. Please check log at $LOG_FILEPATH"$NORMAL
    t3=$(date)
    echo "[pipeline.sh::${t3}]  Finished running buttery-eel on $SLOW5_FILEPATH"

    echo "[pipeline.sh::${t3}]  Converting $UNALN_SAM_FILEPATH to $FASTQ_FILEPATH"
    ${SAMTOOLS} fastq -@ 36 -TMM,ML $UNALN_SAM_FILEPATH > $FASTQ_FILEPATH 2> $LOG_FILEPATH || die $RED"Converting $UNALN_SAM_FILEPATH to $FASTQ_FILEPATH failed. Please check log at $LOG_FILEPATH"$NORMAL
    t4=$(date)
    echo "[pipeline.sh::${t4}]  Finished converting $UNALN_SAM_FILEPATH to $FASTQ_FILEPATH"

    echo "[pipeline.sh::${t4}]  Running minimap2 on $FASTQ_FILEPATH"
    ${MINIMAP2} -t 36 -ax map-ont --sam-hit-only -Y -y --secondary=no $REFIDX $FASTQ_FILEPATH > $UNSORTED_BAM_FILEPATH 2> $LOG_FILEPATH || die $RED"Running minimap2 on $FASTQ_FILEPATH failed. Please check log at $LOG_FILEPATH"$NORMAL
    t5=$(date)
    echo "[pipeline.sh::${t5}]  Finished running minimap2 on $FASTQ_FILEPATH"

    echo "[pipeline.sh::${t5}]  Sorting $UNSORTED_BAM_FILEPATH to $BAM_FILEPATH"
    ${SAMTOOLS} sort -@ 36 -o $BAM_FILEPATH $UNSORTED_BAM_FILEPATH 2> $LOG_FILEPATH || die $RED"Sorting $UNSORTED_BAM_FILEPATH to $BAM_FILEPATH failed. Please check log at $LOG_FILEPATH"$NORMAL
    t6=$(date)
    echo "[pipeline.sh::${t6}]  Finished sorting $UNSORTED_BAM_FILEPATH to $BAM_FILEPATH"

    echo "[pipeline.sh::${t6}]  Indexing $BAM_FILEPATH"
    ${SAMTOOLS} index -@ 36 $BAM_FILEPATH 2> $LOG_FILEPATH || die $RED"Indexing $BAM_FILEPATH failed. Please check log at $LOG_FILEPATH"$NORMAL
    t7=$(date)
    echo "[pipeline.sh::${t7}]  Finished indexing $BAM_FILEPATH"
    
    END_TIME=$(date)

    echo "[pipeline.sh::${END_TIME}]  Finished pipeline for $FILE modbam: $BAM_FILEPATH"
   
    # echo "[pipeline.sh] $P5_FILEPATH" >> $TMP_FILE
    # echo -e "${P5_FILEPATH}\t${SLOW5_FILEPATH}\t${START_TIME}\t${END_TIME}" >> ${LOG}
)&
    ((counter++))
    if [ $counter -ge $MAX_PROC ]; then
        echo "[pipeline.sh] Waiting for $counter jobs to finish."
        wait
        counter=0
    fi
done
wait
