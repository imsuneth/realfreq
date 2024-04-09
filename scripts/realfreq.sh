#!/bin/bash
# set -x

RED='\033[0;31m'
NC='\033[0m' # No Color

SCRIPT_NAME=$(basename $0)
info() {
    echo -e $GREEN"[$SCRIPT_NAME] $1"$NC
}

error() {
    echo -e $RED"[$SCRIPT_NAME] $1"$NC >&2
}

die() {
	error "$1"
	exit 1
}

usage() {
    echo "Usage: $0 [-h] [-y] [-b] -g <guppy_bin> -r <reference> -i <reference_index> -m <model> -d <monitor_dir>"
    echo "  -h  Show help message"
    echo "  -y  Say yes to all prompts"
    echo "  -b  Bedmethyl output"
    echo "  -g  Path to guppy binary"
    echo "  -r  Path to reference fasta"
    echo "  -i  Path to reference index"
    echo "  -m  Model name"
    echo "  -d  Directory to monitor for blow5 files"
    echo
}

# download realp2s if not present
if [ ! -d scripts/realtime-p2s ]; then
    mkdir -p scripts/realtime-p2s/monitor
    echo "Downloading realp2s"
    wget -O scripts/realtime-p2s/realp2s.sh https://raw.githubusercontent.com/Psy-Fer/blue-crab/main/scripts/realtime-p2s/realp2s.sh
    wget -O scripts/realtime-p2s/pipeline.sh https://raw.githubusercontent.com/Psy-Fer/blue-crab/main/scripts/realtime-p2s/pipeline.sh
    wget -O scripts/realtime-p2s/monitor/monitor.sh https://raw.githubusercontent.com/Psy-Fer/blue-crab/main/scripts/realtime-p2s/monitor/monitor.sh
    wget -O scripts/realtime-p2s/monitor/ensure.sh https://raw.githubusercontent.com/Psy-Fer/blue-crab/main/scripts/realtime-p2s/monitor/ensure.sh
    wget -O scripts/realtime-p2s/monitor/simulator.sh https://raw.githubusercontent.com/Psy-Fer/blue-crab/main/scripts/realtime-p2s/monitor/simulator.sh
    chmod +x scripts/realtime-p2s/realp2s.sh
    chmod +x scripts/realtime-p2s/pipeline.sh
    chmod +x scripts/realtime-p2s/monitor/monitor.sh
    chmod +x scripts/realtime-p2s/monitor/ensure.sh
    chmod +x scripts/realtime-p2s/monitor/simulator.sh
fi

say_yes=false
guppy_bin_set=false
reference_set=false
reference_index_set=false
model_set=false
monitor_dir_set=false
bedmethyl_out=false

while getopts "hybg:r:i:m:d:" opt; do
    case $opt in
        d) monitor_dir_set=true; MONITOR_DIR=$OPTARG;;
        g) guppy_bin_set=true; GUPPY_BIN=$OPTARG;;
        r) reference_set=true; REF=$OPTARG;;
        i) reference_index_set=true; REFIDX=$OPTARG;;
        m) model_set=true; MODEL=$OPTARG;;
        h) usage; exit 0;;
        y) say_yes=true;;
        b) bedmethyl_out=true;;
        \?) error "Invalid option: -$OPTARG" >&2
            exit 1;;
        :) error "Option -$OPTARG requires an argument." >&2
            exit 1;;
        *) usage; exit 1;;
    esac
done

if [ $guppy_bin_set = false ]; then
    error "guppy_bin is not set."; usage; exit 1
fi

if [ $reference_set = false ]; then
    error "reference is not set."; usage; exit 1
fi

if [ $reference_index_set = false ]; then
    error "reference_index is not set."; usage; exit 1
fi

if [ $model_set = false ]; then
    error "model is not set."; usage; exit 1
fi

yes_flag=""
if [ $say_yes=true ]; then
    yes_flag="-y"
fi

bed_flag=""
if [ $bedmethyl_out = true ]; then
    bed_flag="-b"
fi

SCRIPT_PATH="$( cd "$(dirname "$0")" ; pwd -P )" # Script's current path
if [ ! -f "$SCRIPT_PATH/pipeline.sh" ]; then
    info "pipeline.sh not found."
    exit 1
fi
PIPELINE=$SCRIPT_PATH/pipeline.sh

mkdir -p $MONITOR_DIR
mkdir -p $MONITOR_DIR/sam
chown -R $USER $MONITOR_DIR

SLOW5TOOLS=slow5tools
BLUECRAB=blue-crab
BUTTERY_EEL=buttery-eel
MINIMAP2=minimap2
SAMTOOLS=samtools
REALFREQ=realfreq
INOTIFYWAIT=inotifywait

command -v $INOTIFYWAIT &> /dev/null || die $RED"$INOTIFYWAIT command not found."$NC
command -v $SLOW5TOOLS &> /dev/null || die $RED"$SLOW5TOOLS command not found."$NC
command -v $BLUECRAB &> /dev/null || die $RED"$BLUECRAB command not found."$NC
command -v $BUTTERY_EEL &> /dev/null || die $RED"$BUTTERY_EEL command not found."$NC
command -v $MINIMAP2 &> /dev/null || die $RED"$MINIMAP2 command not found."$NC
command -v $SAMTOOLS &> /dev/null || die $RED"$SAMTOOLS command not found."$NC
command -v $REALFREQ &> /dev/null || die $RED"$REALFREQ command not found."$NC
command -v $BLUECRAB --version &> /dev/null || die $RED"$BLUECRAB command not found."$NC

EEL="$BUTTERY_EEL -g $GUPPY_BIN --port 5000 --use_tcp"

REALP2S=scripts/realtime-p2s/realp2s.sh

# required for bluecrab
export SLOW5TOOLS=$SLOW5TOOLS
export BLUECRAB=$BLUECRAB
export REALP2S_AUTO=0

pipeline() {
    while read blow5; do
        echo $(date) "Starting pipeline for $blow5" >> $SCRIPT_LOG
        START_TIME=$(date)
        echo -e "$START_TIME\t$blow5" >> $PIPELINE_LOG_ATTEMPTED
        ($PIPELINE -b $blow5 -g $GUPPY_BIN -r $REF -i $REFIDX -m $MODEL -o $MONITOR_DIR/sam 2>> $SCRIPT_LOG) || (echo -e "$(date)\t$blow5" >> $PIPELINE_LOG_FAILED && continue)
        END_TIME=$(date)
        echo -e "$END_TIME\t$blow5" >> $PIPELINE_LOG_DONE
        echo -e "$blow5\t$START_TIME\t$END_TIME" >> $PIPELINE_LOG_START_END
    done
}

catch_blow5() {
    while read line; do
        echo $(date) $line >> $SCRIPT_LOG
        if [[ $line == *"Finished converting"* ]]; then
            blow5=$(echo $line | grep -oP '(?<=to ).*')
            #check file extension
            if [[ $blow5 == *.blow5 ]]; then
                echo $blow5
            fi
            
        fi
    done
}

catch_bam() {
    while read line; do
        echo $(date) $line >> $SCRIPT_LOG
        if [[ $line == *"Finished"* ]]; then
            bam=$(echo $line | grep -oP '(?<=bam-file: ).*')
            #check file extension 
            if [[ $bam == *.bam ]]; then
                echo $bam
            fi
        fi
    done
}

REALFREQ_PROCESSED_LOG=$MONITOR_DIR/realfreq_processed.log
SCRIPT_LOG=$MONITOR_DIR/realfreq.log
PIPELINE_LOG_ATTEMPTED=$MONITOR_DIR/realfreq_pipeline_attempted.log
PIPELINE_LOG_DONE=$MONITOR_DIR/realfreq_pipeline_done.log
PIPELINE_LOG_START_END=$MONITOR_DIR/realfreq_pipeline_start_end_trace.log
PIPELINE_LOG_FAILED=$MONITOR_DIR/realfreq_pipeline_failed.log

clear_logs(){
    test -e $SCRIPT_LOG && rm $SCRIPT_LOG # Empty log file
    test -e $PIPELINE_LOG_ATTEMPTED && rm $PIPELINE_LOG_ATTEMPTED # Empty log file
    test -e $PIPELINE_LOG_DONE && rm $PIPELINE_LOG_DONE # Empty log file
    test -e $PIPELINE_LOG_START_END && rm $PIPELINE_LOG_START_END # Empty log file
    test -e $PIPELINE_LOG_FAILED && rm $PIPELINE_LOG_FAILED # Empty log file
}

#if say_yes is true, clear logs
if [ "$say_yes"=true ]; then
    clear_logs
else
    read -p "$0 Clear logs? [y/n] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        clear_logs
    fi
fi

$REALP2S -m $MONITOR_DIR $yes_flag | catch_blow5 | pipeline | catch_bam | $REALFREQ -r $REF -o $MONITOR_DIR/methfreq.tsv -l $REALFREQ_PROCESSED_LOG $yes_flag $bed_flag
