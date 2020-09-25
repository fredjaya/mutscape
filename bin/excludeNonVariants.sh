#!/bin/bash

GATK=/home/meep/Desktop/Biocomputing/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
REF=$DIR/GCF_003254395.2_Amel_HAv3.1_genomic.fa

if [[ $# != 2 ]]; then
  echo "Usage: $0 [input .vcf] [output .vcf]"
  exit 1
fi

# Run GATK to remove sites with "NO_VARIATION" 
java -jar $GATK -T SelectVariants \
		-R $REF \
		-V $1 \
		--excludeNonVariants \
		-o $2
