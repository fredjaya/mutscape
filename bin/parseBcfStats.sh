#!/bin/bash

INFO=('SN', 'QUAL', 'DP')
VCF=$1

if [[ $# != 1 ]]; then
  echo "Usage: $0 [.vcf]"
  exit 1
fi

# Output SN
grep "# SN" $1 | sed -n '2p' | tr -d '[:punct:]' | tr -d '[:digit:]' | awk 'BEGIN{OFS="\t"} {print $2, $3, $4}' > $1.SN 
grep ^SN $1 | awk -F'\t' 'BEGIN{OFS="\t"} {print $2,$3,$4}' | sed "s/number of //g" | sed "s/\://g" >> $1.SN

# Output QUAL
grep "# QUAL" $1 | sed -n '2p' | sed 's/1st/first/g' | tr -d '\[\]' | tr -d '[:digit:]' | awk -F'\t' 'BEGIN{OFS="\t"} {print $2, $3, $4, $5, $6, $7}' > $1.QUAL
grep '^QUAL' $1 | awk -F'\t' 'BEGIN{OFS="\t"} {print $2, $3, $4, $5, $6, $7}' >> $1.QUAL 

# Output DP
grep "# DP" $1 | sed -n '2p' | tr -d '\[\]\#' | tr -d '[:digit:]' | awk -F'\t' 'BEGIN{OFS="\t"} {print $2, $3, $4, $5, $6, $7}' > $1.DP
grep "^DP" $1 | awk -F'\t' 'BEGIN{OFS="\t"} {print $4, $3, $4, $5, $6, $7}' >> $1.DP

# NOTE: 501 == '>500'
sed -i 's/>500/501/g' $1.DP 
