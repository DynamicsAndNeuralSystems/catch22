#!/bin/bash

# Run catch22 feature extraction over every test-data subdirectory.
# Each subdirectory of testData holds a set of .txt time-series files;
# outputs are written back alongside the inputs with an "_output" suffix.

testdir="$(dirname -- "$0")"
runner="${testdir}/../C/runAllTS.sh"

for dir in "${testdir}"/*/
do
    dir="${dir%/}"
    echo "=== Running tests in $(basename "$dir") ==="
    "$runner" -i "$dir" -o "$dir" -a "_output" -s 1
done
