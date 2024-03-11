#!/bin/bash

# Configuration
GUPPY_BIN=/home/suneth/tools/ont-dorado-server/bin
MODEL_CONFIG="dna_r10.4.1_e8.2_400bps_5khz_modbases_5hmc_5mc_cg_fast.cfg"
REF=/media/suneth/Data/Workspace/validate/data/ref/hg38noAlt.fa
# End of Configuration

if [ -z "$1" ]; then
	echo "Usage: $0 <experiment_dir>"
	exit 1
fi

EXP_DIR=$1

RED='\033[0;31m'
NC='\033[0m' # No Color

die() {
	echo -e "$1" >&2
	echo
	exit 1
}

check_prev_cmd() {
    if [ $? -ne 0 ]; then
        die $RED"$1"$NC
    fi
}

check_cmd() {
    if ! command -v $1 &> /dev/null
    then
        die $RED"$1 command not found."$NC
    fi
}

SLOW5TOOLS=slow5tools
BLUECRAB=blue-crab
BUTTERY_EEL=buttery-eel
MINIMAP2=minimap2
SAMTOOLS=samtools
F5C=f5c
REALFREQ=realfreq
INOTIFYWAIT=inotifywait

check_cmd $INOTIFYWAIT
check_cmd $SLOW5TOOLS
check_cmd $BUTTERY_EEL
check_cmd $MINIMAP2
check_cmd $SAMTOOLS
check_cmd $F5C
check_cmd $REALFREQ

$BLUECRAB --version
check_prev_cmd "Error in running blue-crab. Make sure virtual environment is activated"

EEL="$BUTTERY_EEL -g $GUPPY_BIN --port 5000 --use_tcp"

# download realp2s if not present
REALP2S=scripts/realp2s.sh
if [ ! -f $REALP2S ]; then
    echo "Downloading realp2s"
    wget -O $REALP2S https://raw.githubusercontent.com/Psy-Fer/blue-crab/main/scripts/realtime-p2s/realp2s.sh
    chmod +x $REALP2S
fi

# required for bluecrab
export SLOW5TOOLS=$SLOW5TOOLS
export BLUECRAB=$BLUECRAB
export REALP2S_AUTO=0

pipeline() {
    while read blow5file; do
        readname=$(basename -s .blow5 $blow5file)
        
        echo "Basecalling $blow5file"
        unalignedsamfile="$EXP_DIR/sam/$readname.unaln.sam"
        $EEL--call_mods --config $MODEL_CONFIG -i $blow5file -o $unalignedsamfile --device cuda:all
        check_prev_cmd "Error in basecalling $blow5file"

        echo "Aligning $unalignedsamfile and get bam-file"
        bamfile="$EXP_DIR/sam/$readname.aln.mod.bam"
        samtools fastq -@ 36 -TMM,ML $unalignedsamfile | \
        minimap2 -t 36 -ax map-ont -y $REF - | \
        samtools view -@ 36 -Sb - | \
        samtools sort -@ 36 - > $bamfile
        check_prev_cmd "Error in aligning $unalignedsamfile"

        echo "Indexing $bamfile"
        samtools index -@ 36 $bamfile
        check_prev_cmd "Error in indexing $bamfile"

        echo "Finished. bam-file: $bamfile"
        
        sleep 10
    done
}

catch_blow5() {
    while read line; do
        echo $(date) $line >> $EXP_DIR/realp2s.log
        if [[ $line == *"Finished converting"* ]]; then
            blow5file=$(echo $line | grep -oP '(?<=to ).*')
            #check file extension
            if [[ $blow5file == *.blow5 ]]; then
                echo $blow5file
            fi
            
        fi
    done
}

catch_bam() {
    while read line; do
        echo $(date) $line >> $EXP_DIR/pipeline.log
        if [[ $line == *"Finished"* ]]; then
            bamfile=$(echo $line | grep -oP '(?<=bam-file: ).*')
            #check file extension 
            if [[ $bamfile == *.bam ]]; then
                echo $bamfile
            fi
        fi
    done
}

echo > $EXP_DIR/realp2s.log
echo > $EXP_DIR/pipeline.log
$REALP2S -m $EXP_DIR | catch_blow5 | pipeline | catch_bam | $REALFREQ > $EXP_DIR/freq.tsv
