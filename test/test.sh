#!/bin/bash
set -x
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

mkdir -p test/tmp || die "Creating the test/tmp directory failed"

if [ ! -f test/tmp/genome_chr22.fa ]; then
    wget  -N -O test/tmp/genome_chr22.fa.gz "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz" || die "Downloading the genome chr22 failed"
    gzip -d test/tmp/genome_chr22.fa.gz || die "Unzipping the genome chr22 failed"
fi

echo "Test 1: tsv output"
echo_bams | ex ./realfreq -r test/tmp/genome_chr22.fa -o test/tmp/test1.tsv >> test/tmp/test1.log || die "Test 1 failed running realfreq"
sort -k1,1 -k2,2n test/expected/test1.tsv > test/tmp/test1.expected.sorted.tsv
sort -k1,1 -k2,2n test/tmp/test1.tsv > test/tmp/test1.sorted.tsv
diff test/tmp/test1.expected.sorted.tsv test/tmp/test1.sorted.tsv || die "Test 1: diff failed"

echo "Test 2: bedmethyl output"
echo_bams | ex ./realfreq -r test/tmp/genome_chr22.fa -o test/tmp/test2.bedmethyl -b >> test/tmp/test2.log || die "Test 2 failed running realfreq"
sort -k1,1 -k2,2n test/expected/test2.bedmethyl > test/tmp/test2.expected.sorted.bedmethyl
sort -k1,1 -k2,2n test/tmp/test2.bedmethyl > test/tmp/test2.sorted.bedmethyl
diff test/tmp/test2.expected.sorted.bedmethyl test/tmp/test2.sorted.bedmethyl || die "Test 2: diff failed"

echo "TESTS PASSED!"