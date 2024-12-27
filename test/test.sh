#!/bin/bash

BLUE='\033[0;34m'
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

# terminate script
die() {
	echo "$1" >&2
	echo
	exit 1
}

if [ "$1" = 'mem' ]; then
    mem=1
else
    mem=0
fi

ex() {
    if [ $mem -eq 1 ]; then
        valgrind --leak-check=full --error-exitcode=1 "$@"
    else
        "$@"
    fi
}

echo_bams() {
    echo "test/example-ont-part1.bam"
    sleep 1
    echo "test/example-ont-part2.bam"
}

echo_tsv() {
    echo "test/nanopolish.tsv"
}

mkdir -p test/tmp || die "Creating the test/tmp directory failed"

if [ ! -f test/tmp/genome_chr22.fa ]; then
    wget  -N -O test/tmp/genome_chr22.fa.gz "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz" || die "Downloading the genome chr22 failed"
    gzip -d test/tmp/genome_chr22.fa.gz || die "Unzipping the genome chr22 failed"
fi

echo "Test 1: tsv output"
echo_bams | ex ./realfreq -d test/tmp/test1.dump.tmp -o test/tmp/test1.tsv test/tmp/genome_chr22.fa > test/tmp/test1.log || die "Test 1 failed running realfreq"
sort -k1,1 -k2,2n test/expected/test1.tsv > test/tmp/test1.expected.sorted.tsv
sort -k1,1 -k2,2n test/tmp/test1.tsv > test/tmp/test1.sorted.tsv
diff -q test/tmp/test1.expected.sorted.tsv test/tmp/test1.sorted.tsv || die "Test 1: diff failed"
rm test/tmp/test1.tsv test/tmp/test1.log test/tmp/test1.expected.sorted.tsv test/tmp/test1.sorted.tsv test/tmp/test1.dump.tmp

echo "Test 2: bedmethyl output"
echo_bams | ex ./realfreq -d test/tmp/test2.dump.tmp -o test/tmp/test2.bedmethyl -b test/tmp/genome_chr22.fa > test/tmp/test2.log || die "Test 2 failed running realfreq"
sort -k1,1 -k2,2n test/expected/test2.bedmethyl > test/tmp/test2.expected.sorted.bedmethyl
sort -k1,1 -k2,2n test/tmp/test2.bedmethyl > test/tmp/test2.sorted.bedmethyl
diff -q test/tmp/test2.expected.sorted.bedmethyl test/tmp/test2.sorted.bedmethyl || die "Test 2: diff failed"
rm test/tmp/test2.bedmethyl test/tmp/test2.log test/tmp/test2.expected.sorted.bedmethyl test/tmp/test2.sorted.bedmethyl test/tmp/test2.dump.tmp

echo "Test 3a: dumping"
echo "test/example-ont.bam" | ex ./realfreq  -d test/tmp/test3.dump.tmp -o test/tmp/test3.out1.tsv test/tmp/genome_chr22.fa > test/tmp/test3.log1.log || die "Test 3a failed running realfreq - run"
echo "Test 3b: loading"
echo "" | ex ./realfreq -d test/tmp/test3.dump.tmp -o test/tmp/test3.out2.tsv -r test/tmp/genome_chr22.fa > test/tmp/test3.log2.log || die "Test 3b failed running realfreq - resume"
sort -k1,1 -k2,2n test/tmp/test3.out1.tsv > test/tmp/test3.out1.sorted.tsv
sort -k1,1 -k2,2n test/tmp/test3.out2.tsv > test/tmp/test3.out2.sorted.tsv
diff -q test/tmp/test3.out1.sorted.tsv test/tmp/test3.out2.sorted.tsv || die "Test 3: diff failed"
rm test/tmp/test3.dump.tmp test/tmp/test3.out1.tsv test/tmp/test3.out2.tsv test/tmp/test3.log1.log test/tmp/test3.log2.log test/tmp/test3.out1.sorted.tsv test/tmp/test3.out2.sorted.tsv

echo "Test 4: nanopore tsv input"
echo_tsv | ex ./realfreq --tsv -d test/tmp/test4.dump.tmp -o test/tmp/test4.tsv > test/tmp/test4.log || die "Test 4 failed running realfreq"
sort -k1,1 -k2,2n test/expected/test4.tsv > test/tmp/test4.expected.sorted.tsv
sort -k1,1 -k2,2n test/tmp/test4.tsv > test/tmp/test4.sorted.tsv
diff -q test/tmp/test4.expected.sorted.tsv test/tmp/test4.sorted.tsv || die "Test 4: diff failed"

echo -e "${GREEN}ALL TESTS PASSED !${NC}"