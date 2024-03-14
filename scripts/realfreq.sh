#!/bin/bash
# set -x
# check if env variables are set
if [ -z "$GUPPY_BIN" ]; then
    echo "GUPPY_BIN is not set."
    exit 1
fi

if [ -z "$REF" ]; then
    echo "REF is not set."
    exit 1
fi

if [ -z "$REFIDX" ]; then
    echo "REFIDX is not set."
    exit 1
fi

if [ -z "$MODEL" ]; then
    echo "MODEL is not set."
    exit 1
fi

if [ -z "$1" ]; then
	echo "Usage: $0 <experiment_dir>"
	exit 1
fi
EXP_DIR=$1


SCRIPT_PATH="$( cd "$(dirname "$0")" ; pwd -P )" # Script's current path
if [ ! -f "$SCRIPT_PATH/pipeline.sh" ]; then
    echo "pipeline.sh not found."
    exit 1
fi
PIPELINE=$SCRIPT_PATH/pipeline.sh

mkdir -p $EXP_DIR
mkdir -p $EXP_DIR/sam
chown -R $USER $EXP_DIR

RED='\033[0;31m'
NC='\033[0m' # No Color

die() {
	echo -e "$1" >&2
	echo
	exit 1
}

SLOW5TOOLS=slow5tools
BLUECRAB=blue-crab
BUTTERY_EEL=buttery-eel
MINIMAP2=minimap2
SAMTOOLS=samtools
F5C=f5c
REALFREQ=realfreq
INOTIFYWAIT=inotifywait

command -v $INOTIFYWAIT &> /dev/null || die $RED"$INOTIFYWAIT command not found."$NC
command -v $SLOW5TOOLS &> /dev/null || die $RED"$SLOW5TOOLS command not found."$NC
command -v $BLUECRAB &> /dev/null || die $RED"$BLUECRAB command not found."$NC
command -v $BUTTERY_EEL &> /dev/null || die $RED"$BUTTERY_EEL command not found."$NC
command -v $MINIMAP2 &> /dev/null || die $RED"$MINIMAP2 command not found."$NC
command -v $SAMTOOLS &> /dev/null || die $RED"$SAMTOOLS command not found."$NC
command -v $F5C &> /dev/null || die $RED"$F5C command not found."$NC
command -v $REALFREQ &> /dev/null || die $RED"$REALFREQ command not found."$NC
command -v $BLUECRAB --version &> /dev/null || die $RED"$BLUECRAB command not found."$NC

EEL="$BUTTERY_EEL -g $GUPPY_BIN --port 5000 --use_tcp"

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

REALP2S=scripts/realtime-p2s/realp2s.sh

# required for bluecrab
export SLOW5TOOLS=$SLOW5TOOLS
export BLUECRAB=$BLUECRAB
export REALP2S_AUTO=0

pipeline() {
    while read blow5; do
        $PIPELINE $blow5 $EXP_DIR/sam
        sleep 10
    done
}

catch_blow5() {
    while read line; do
        echo $(date) $line >> $EXP_DIR/realfreq.log
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
        echo $(date) $line >> $EXP_DIR/realfreq.log
        if [[ $line == *"Finished"* ]]; then
            bam=$(echo $line | grep -oP '(?<=bam-file: ).*')
            #check file extension 
            if [[ $bam == *.bam ]]; then
                echo $bam
            fi
        fi
    done
}

echo > $EXP_DIR/realfreq.log
$REALP2S -m $EXP_DIR | catch_blow5 | pipeline 2>&1 | catch_bam | $REALFREQ -o $EXP_DIR/freq.tsv 2>&1 >> $EXP_DIR/realfreq.log
