# mutscape
Variant calling pipeline to estimate *A. m. capensis* mutation rates.

Took an iterative approach to pipeline development i.e. running nextflow processes individually to generate results, over developing a comprehensive pipeline with singularity (e.g. `nf-desktop` branch).

## Protocol
Following GATK 3.X best practices from Van der Auwera et al. (2013)

### Downloading and indexing the reference genome
These are run prior to nextflow

```
# Download reference genome, generate bwa and .fasta indexes
cd /scratch/Scape/fred/
rsync --copy-links --times --verbose \
rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz \ 
/scratch/Scape/fred/ # Artemis HPC

# Decompress reference to avoid bgzip indexing issues
cd /scratch/Scape/fred/rtc_idx
bgzip -d GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz

# Generate BWA index
qsub tests/scripts/1_bwa_index.sh

# Generate .fasta index 
qsub tests/scripts/2_ref_index.sh

# Create dictionary
qsub scripts/3_create_dict.sh
```

### Preparing sequence data (.fastq to .bam)
These are run through nextflow (`main.nf`). Each process is called by queueing nextflow to run individual processes. Processes need to be specified with `nextflow run --mode <process>`.

```
# 1. Map trimmed reads to reference genome
qsub run_scripts/1_bwaMapReads.sh   

# 2. Convert .sam into .sorted .bam
qsub run_scripts/2_samtoolsSort.sh

# 3. Mark duplicates
qsub run_scripts/3_markDuplicates.sh

# 4. Mark regions for realignment around indels 
qsub run_scripts/4_realignerTargetCreator.sh

# Added tests here for inputting bam/bai files.
# 5. Realign regions around indels
qsub run_scripts/5_indelRealigner
```

### Base Quality Score Recalibration
Bootstrapping known sites according to [documentation](https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-methods-and-algorithms/Base_Quality_Score_Recalibration_(BQSR).md). 

```
# Covariate tables and plots not generated here (TODO).
# ReduceReads skipped as deprecated in GATK3

# 6. Call confident sites for BSQR
qsub run_scripts/6_haplotypeCallerUnrecal

# 7. Generate covariate tables for BQSR
qsub run_scripts/7_bsqrTable

# 8. Recalibrate base quality scores
qsub run_scripts/8_recalibrateBQS.sh
```

### Variant calling
```
# 9. Call SNPs and indels from recalibrated reads 
qsub run_scripts/9_callVariants.sh
```

### Variant introspection and filtering
```
# Merge individual vcf files into a single multiple-sample vcf
bin/combineVariants.sh

# Remove sites where variants are all homozygous for the reference allele (0|0) 
bin/excludeNonVariants.sh

# Remove sites contain all common variants (e.g. all 0/0, 0/1, 1/1)
bin/filterCommonVariants.sh
```
