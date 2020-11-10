#!/bin/bash

BCFTOOLS=/home/meep/Desktop/Biocomputing/bcftools-1.10.2/bcftools

if [[ $# != 1 ]]
then
  echo "Usage: $0 [full path to *.vcf without extension]" 
  exit 1
fi

# Produce bcftools stats, per sample
echo "Running bcftools stats"
$BCFTOOLS stats -s - $1.vcf > $1.stats && \

# Subset per sample SNP counts
echo "Subsetting per-sample SNP counts"
awk -F '\t' '/^\# PSC\t/;/^PSC/' $1.stats | \
sed s/\#\ /''/ | sed -E s/\\[[0-9]+\\]//g \
> $1.psc
