#!/bin/bash

#PBS -P RDS-FSC-Scape-RW
#PBS -l select=1:ncpus=4:mem=16GB
#PBS -l walltime=3:00:00

DIR=/scratch/Scape/fred/rtc_idx
READ_PATH=/project/Scape/Trimmomatic/Paired
SAMPLE=Larv09
READ_GROUP=@RG\\tID:${SAMPLE}\\tSM:${SAMPLE}\\tPL:ILLUMINA\\tLB:lib1\\tPU:unit1

module load bwa/0.7.17

cd $DIR

bwa mem \
    -R $READ_GROUP \
    -t 4 \
    GCF_003254395.2_Amel_HAv3.1_genomic.fna \
    $READ_PATH/${SAMPLE}_R1_paired.fastq.gz \
    $READ_PATH/${SAMPLE}_R2_paired.fastq.gz \
    > ${SAMPLE}_pe.sam
