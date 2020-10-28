#!/bin/bash

DATA=/home/meep/Desktop/People/fred/Dropbox/meep/bee/01_data/2010_gvcf
REF=/home/meep/Desktop/People/fred/Dropbox/meep/bee/02_working/2009_filter_vcf/GCF_003254395.2_Amel_HAv3.1_genomic.fa

~/Desktop/Biocomputing/gatk-4.1.8.1/gatk --java-options "-Xmx32g -Xms32g" GenomicsDBImport \
	-R $REF \
	-V $DATA/Larv01.g.vcf \
	-V $DATA/Larv02.g.vcf \
	-V $DATA/Larv03.g.vcf \
	-V $DATA/Larv04.g.vcf \
	-V $DATA/Larv05.g.vcf \
	-V $DATA/Larv06.g.vcf \
	-V $DATA/Larv07.g.vcf \
	-V $DATA/Larv08.g.vcf \
	-V $DATA/Larv09.g.vcf \
	-V $DATA/Larv10.g.vcf \
	-V $DATA/Larv11.g.vcf \
	-V $DATA/Larv12.g.vcf \
	-V $DATA/Larv13.g.vcf \
	-V $DATA/Larv14.g.vcf \
	-V $DATA/Worker.g.vcf \
	--genomicsdb-workspace-path ./genomicsdb \
	-L intervals.list
