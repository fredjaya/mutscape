#!/bin/bash

INPUT_VCF=$1
OUT_STATS=$2

if [[ $# != 2 ]]
then
  echo "Usage: $0 [.vcf] [.stats]"
  exit 1
fi

# Get count of variants, generate .stats
echo "Running bcftools: generating .stats"
bin/countVariants.sh $INPUT_VCF $OUT_STATS && 

# Parse the output .stats file for plotting
echo "Parsing .stats"
bin/parseBcfStats.sh $OUT_STATS &&

# Plot .stats
echo "Plotting"
Rscript bin/plotBcfStats.R $OUT_STATS

