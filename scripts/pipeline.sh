#!/bin/bash
set -x
RED='\033[0;31m'
NC='\033[0m' # No Color

LOG=realfreq_pipline.log

info() {
    echo -e $GREEN"[$0] $1"$NC
}

error() {
    echo -e $RED"[$0] $1"$NC >&2
}

die() {
	error "$1"
	exit 1
}

log() {
    echo -e "read:$read, step: $1, status: $2, start: $3, end: $4" >> $LOG
}

usage() {
    (
    echo "Usage: $0 [-h] [-y] -b <blow5> -o <output_dir> -g <guppy_bin> -r <reference> -i <reference_index> -m <model> -d <monitor_dir>"
    echo "  -h  Show help message"
    echo "  -y  Say yes to all prompts"
    echo "  -b  Path to blow5 file"
    echo "  -o  Output directory"
    echo "  -g  Path to guppy binary"
    echo "  -r  Path to reference fasta"
    echo "  -i  Path to reference index"
    echo "  -m  Model name"
    echo
    ) >&2
}

blow5_file_set=false
out_dir_set=false
guppy_bin_set=false
reference_set=false
reference_index_set=false
model_set=false

while getopts "hyb:o:g:r:i:m:" opt; do
    case $opt in
        b) blow5_file_set=true; BLOW5=$OPTARG;;
        o) out_dir_set=true; OUT_DIR=$OPTARG;;
        g) guppy_bin_set=true; GUPPY_BIN=$OPTARG;;
        r) reference_set=true; REF=$OPTARG;;
        i) reference_index_set=true; REFIDX=$OPTARG;;
        m) model_set=true; MODEL=$OPTARG;;
        h) usage; exit 0;;
        \?) error "Invalid option: -$OPTARG" >&2
            exit 1;;
        :) error "Option -$OPTARG requires an argument." >&2
            exit 1;;
        *) usage; exit 1;;
    esac
done

if [ $blow5_file_set = false ]; then
    error "blow5 is not set."; usage; exit 1
fi

if [ $out_dir_set = false ]; then
    error "Output directory is not set."; usage; exit 1
fi

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

mkdir -p $OUT_DIR
chown -R $USER $OUT_DIR

BUTTERY_EEL=buttery-eel
MINIMAP2=minimap2
SAMTOOLS=samtools

command -v $BUTTERY_EEL &> /dev/null || die $RED"$BUTTERY_EEL command not found."$NC
command -v $MINIMAP2 &> /dev/null || die $RED"$MINIMAP2 command not found."$NC
command -v $SAMTOOLS &> /dev/null || die $RED"$SAMTOOLS command not found."$NC

EEL="$BUTTERY_EEL -g $GUPPY_BIN --port 5000 --use_tcp"

read=$(basename -s .blow5 $BLOW5)
unalsam="$OUT_DIR/$read.remora.unaln.sam"
fastq="$OUT_DIR/$read.remora.fastq"
unsortedbam="$OUT_DIR/$read.remora.unsorted.bam"
bam="$OUT_DIR/$read.remora.bam"

t0=$(date)
/usr/bin/time -v $EEL --call_mods --config $MODEL -i $BLOW5 -o $unalsam --device cuda:all || echo "read:$read, step: eel, status: failed, start: $t0, end: $t1"
t1=$(date)
echo "read:$read, step: eel, status: success, start: $t0, end: $t1"
/usr/bin/time -v samtools fastq -@ 36 -TMM,ML $unalsam > $fastq || echo "read:$read, step: samtools, status: failed, start: $t1, end: $t2"
t2=$(date)
echo "read:$read, step: samtools-fastq, status: success, start: $t1, end: $t2"
/usr/bin/time -v minimap2 -t 36 -x map-ont --sam-hit-only -Y -a -y --secondary=no $REFIDX $fastq > $unsortedbam || echo "read:$read, step: minimap2, status: failed, start: $t2, end: $t3"
t3=$(date)
echo "read:$read, step: minimap2, status: success, start: $t2, end: $t3"
/usr/bin/time -v samtools sort -@ 36 $unsortedbam > $bam || echo "read:$read, step: samtools, status: failed, start: $t3, end: $t4"
t4=$(date)
echo "read:$read, step: samtools-sort, status: success, start: $t3, end: $t4"
/usr/bin/time -v samtools index -@ 36 $bam || echo "read:$read, step: samtools, status: failed, start: $t4, end: $t5"
t5=$(date)
echo "read:$read, step: samtools-index, status: success, start: $t4, end: $t5"

echo "Finished. bam-file: $bam"
exit 0