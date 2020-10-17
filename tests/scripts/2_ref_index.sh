#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l walltime=00:01:00

DIR=/scratch/Scape/fred/2008_manual

module load samtools/1.9

cd $DIR

samtools faidx GCF_003254395.2_Amel_HAv3.1_genomic.fna
