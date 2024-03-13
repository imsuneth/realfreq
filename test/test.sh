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
echo_bams | ex ./realfreq -o temp/example-ont-real.tsv || die "Test 1 failed running realfreq"
cat temp/example-ont-real.tsv | sort -n -k2,2 -k10,10 > temp/example-ont-real-s.tsv
diff -q test/example-ont-freq-s.tsv temp/example-ont-real-s.tsv || die "Test 1 diff failed"


echo "TESTS PASSED!"