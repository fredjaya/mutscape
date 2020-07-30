#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=1:mem=20GB
#PBS -l walltime=1:00:00

DIR='/scratch/Scape/fred/rtc_idx'

cd $DIR
module load samtools/1.9

samtools sort \
    -o Larv09_pe_sorted.bam \
    -O bam \
    Larv09_pe.sam
