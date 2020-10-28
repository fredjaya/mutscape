#!/bin/bash

DATA=/home/meep/Desktop/People/fred/Dropbox/meep/bee/01_data/2010_gvcf
REF=/home/meep/Desktop/People/fred/Dropbox/meep/bee/02_working/2009_filter_vcf/GCF_003254395.2_Amel_HAv3.1_genomic.fa

~/Desktop/Biocomputing/gatk-4.1.8.1/gatk --java-options "-Xmx60g" GenotypeGVCFs \
	-R $REF \
	-V gendb://genomicsdb \
	-O joint_genotype.vcf.gz
