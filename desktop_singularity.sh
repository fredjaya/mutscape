#!/bin/bash

WORKDIR=/home/meep/Desktop/People/fred/mutscape

nextflow run main.nf \
	--ref    ${WORKDIR}/data/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
	--fai    ${WORKDIR}/data/GCF_003254395.2_Amel_HAv3.1_genomic.fna.fai \
	--dict   ${WORKDIR}/data/GCF_003254395.2_Amel_HAv3.1_genomic.dict \
	--bwaidx ${WORKDIR}/data/GCF_003254395.2_Amel_HAv3.1_genomic.fna.* \
	--reads  ${WORKDIR}/mutscape/data/*_{R1,R2}_paired.fastq.gz \
	--outDir ${WORKDIR}/results \
	-profile singularity \
	-resume 
