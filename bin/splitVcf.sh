#!/bin/bash

BCFTOOLS=/home/meep/Desktop/Biocomputing/bcftools-1.10.2/bcftools 
GATK=/home/meep/Desktop/Biocomputing/gatk-4.1.8.1/gatk
DIR=/home/meep/Desktop/People/fred/Dropbox/meep/bee/02_working/2010_consolidate
REC=/home/meep/Desktop/People/fred/Dropbox/meep/bee/02_working/2011_recombination

for sample in $BCFTOOLS query -l $DIR/merged_variants.vcf; do
	echo $sample
	#$GATK SelectVariants -R $REF -V $DIR/merged_variants.vcf -O $REC/${sample}.vcf -sn $sample
done

