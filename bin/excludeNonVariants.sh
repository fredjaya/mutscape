#!/bin/bash

GATK=/home/meep/Desktop/Biocomputing/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
REF=/home/meep/Desktop/People/fred/Dropbox/meep/bee/02_working/2005_vcf_test/GCF_003254395.2_Amel_HAv3.1_genomic.fa
DATA_PATH=/home/meep/Desktop/People/fred/Dropbox/meep/bee/01_data
VCF=$DATA_PATH/combinedLarv.vcf

# Run GATK to merge samples
java -jar $GATK -T SelectVariants \
		-R $REF \
		-V $VCF \
		--exclude_sample_name Fdrone \
		--excludeNonVariants \
		-o $DATA_PATH/combined_nonVariants.vcf
