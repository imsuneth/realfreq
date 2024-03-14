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

if [ -z "$1" ] || [ -z "$2" ]; then
	echo "Usage: $0 <blow5_file> <out_dir>"
	exit 1
fi

blow5=$1
OUT_DIR=$2

mkdir -p $OUT_DIR
chown -R $USER $OUT_DIR

RED='\033[0;31m'
NC='\033[0m' # No Color

die() {
	echo -e "$1" >&2
	echo
	exit 1
}

BUTTERY_EEL=buttery-eel
MINIMAP2=minimap2
SAMTOOLS=samtools

command -v $BUTTERY_EEL &> /dev/null || die $RED"$BUTTERY_EEL command not found."$NC
command -v $MINIMAP2 &> /dev/null || die $RED"$MINIMAP2 command not found."$NC
command -v $SAMTOOLS &> /dev/null || die $RED"$SAMTOOLS command not found."$NC

EEL="$BUTTERY_EEL -g $GUPPY_BIN --port 5000 --use_tcp"


read=$(basename -s .blow5 $blow5)
unalsam="$OUT_DIR/$read.remora.unaln.sam"
bam="$OUT_DIR/$read.remora.bam"

$EEL --call_mods --config $MODEL -i $blow5 -o $unalsam --device cuda:all || die "Error in basecalling $blow5"
samtools fastq -@ 36 -TMM,ML $unalsam | minimap2 -t 36 -x map-ont --sam-hit-only -Y -a -y --secondary=no $REFIDX - | samtools sort -@ 36 - > $bam || die "Error in aligning $unalsam"
samtools index -@ 36 $bam || die "Error in indexing $bam"

echo "Finished. bam-file: $bam"
exit 0