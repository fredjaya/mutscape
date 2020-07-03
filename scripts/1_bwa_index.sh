#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l walltime=00:5:00

DIR='/scratch/Scape/fred/rtc_idx'

module load bwa/0.7.17

cd $DIR

bwa index GCF_003254395.2_Amel_HAv3.1_genomic.fna
