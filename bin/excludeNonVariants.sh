#!/bin/bash

DIR=/home/meep/Desktop/People/fred/Dropbox/meep/bee/02_working/2009_filter_vcf
GATK=/home/meep/Desktop/Biocomputing/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
REF=$DIR/GCF_003254395.2_Amel_HAv3.1_genomic.fa
VCF=$DIR/combined_exFiltered_exCommons.vcf

# Run GATK to merge samples
java -jar $GATK -T SelectVariants \
		-R $REF \
		-V $VCF \
		--excludeNonVariants \
		-o $DIR/combined_exFiltered_exCommons_exNonVar.vcf
