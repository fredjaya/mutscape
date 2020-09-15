#!/bin/bash

BCFTOOLS=/home/meep/Desktop/Biocomputing/bcftools-1.10.2/bcftools

# Get count of variants
if [[ $# != 2 ]]
then
  echo "Usage: $0 [vcf] [out]" 
  exit 1
fi

$BCFTOOLS stats $1 > $2
