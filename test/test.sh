#!/bin/bash

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

mkdir -p temp || die "Failed to create temp directory"

echo "Test 1"
echo_bams | ex ./realfreq -o temp/example-ont-freq.tsv || die "Test 1 failed running realfreq"
cat temp/example-ont-freq.tsv | tail -n +2 - | sort -n -k2,2 -k10,10 > temp/example-ont-freq-actual.tsv
cat test/example-ont-freq.tsv | tail -n +2 - | sort -n -k2,2 -k10,10 > temp/example-ont-freq-expected.tsv
diff -q temp/example-ont-freq-expected.tsv temp/example-ont-freq-actual.tsv || die "Test 1 failed"


echo "TESTS PASSED!"