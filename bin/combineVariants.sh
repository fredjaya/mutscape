#!/bin/bash

GATK=/home/meep/Desktop/Biocomputing/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
REF=/home/meep/Desktop/People/fred/Dropbox/meep/bee/02_working/2005_vcf_test/GCF_003254395.2_Amel_HAv3.1_genomic.fa
DATA_PATH=/home/meep/Desktop/People/fred/Dropbox/meep/bee/01_data/

# Create list of .vcf files to merge
if [ ! -f "vcf.list" ]
then
    echo "Creating vcf.list"
    ls $DATA_PATH/*_raw_variants.vcf > vcf.list
fi

# Run GATK to merge samples
java -jar $GATK -T CombineVariants \
		-R $REF \
		-V vcf.list \
		-o combinedLarv.vcf
