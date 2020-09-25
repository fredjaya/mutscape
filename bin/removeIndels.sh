#!/bin/bash

GATK=/home/meep/Desktop/Biocomputing/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
REF=/home/meep/Desktop/People/fred/Dropbox/meep/bee/02_working/2005_vcf_test/GCF_003254395.2_Amel_HAv3.1_genomic.fa
INPUT=$1
OUTPUT=$2

if [[ $# != 2 ]]; then
  echo "Usage: $0 [input .vcf] [output .vcf]"
  exit 1
fi                                

# Run GATK to retain SNPs only
java -jar $GATK -T SelectVariants \
		-R $REF \
		-V $INPUT \
		-selectType SNP \
		-o $OUTPUT
