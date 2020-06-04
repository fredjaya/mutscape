#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=1:mem=16GB
#PBS -l walltime=12:00:00

module load bwa
cd /scratch/Scape/fred
bwa index GCA_003254395.2_Amel_HAv3.1_genomic.fna.gz
