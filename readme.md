# mutscape

Proof-of-concept pipeline to estimate _A. m. capensis_ mutation rates

## Aims
- To familiarise self with the processes and filter involved
- Produce a filtered `.vcf` file of Larv09 for comparison with Fdrone
- Produce full pipeline from `.fastq` to `.vcf` -> develop into multi-sample pipeline
	- I/O files
	- Resource usage

## Template
Following GATK 3.X best practices from Van der Auwera et al. (2013)

## Protocol

### Preparing sequence data (`.fastq` to `.bam`)
```
# Download the reference genome for _Apis mellifera_
cd /scratch/Scape/fred/
rsync --copy-links --times --verbose \
rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz \ 
/scratch/Scape/fred/ # Artemis HPC

# Decompress reference to avoid bgzip indexing issues
cd /scratch/Scape/fred/rtc_idx
bgzip -d GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz

# Generate BWA index
qsub scripts/1_bwa_index.sh

# Generate fasta index 
qsub scripts/2_ref_index.sh

# Create dictionary
qsub scripts/3_create_dict.sh

# Map reads
qsub scripts/4_bwa_mem.sh

# Sort and convert to bam
qsub scripts/5_sortsam.sh

# Mark duplicates
qsub scripts/6_markdups.sh

# Realign target creator
qsub scripts/7_rtc.sh

# Realign reads around indels
qsub scripts/8_realign_indels.sh
```

### Call SNPs on unrecalibrated data for BSQR
According to [documentation](https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-methods-and-algorithms/Base_Quality_Score_Recalibration_(BQSR).md).

ReduceReads skipped as deprecated in GATK3

```
# HaplotypeCaller - call all sites in discovery mode
qsub scripts/9_unrecal_hap_call.sh 
