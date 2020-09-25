#!/bin/bash

# https://stackoverflow.com/a/28905692

if [[ $# != 1 ]]; then
  echo "Usage: $0 [.vcf]"
  exit 1
fi

echo "Counting number of genotypes"

# Count number of called SNPs e.g. 0/0, 0/1, 1/1
grep -oE "[[:digit:]]+/[[:digit:]]+" $1 | sort | uniq -c | \
  awk '{sum += $1} END {print "Total called SNPs:", sum}'

# Count number of non-calls e.g. ./.
grep -oE "\.\/\." $1 | sort | uniq -c | \
  awk '{print "Total NO_CALL SNPs:", $1}'

