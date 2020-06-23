# mutscape  

Proof of concept pipeline to estimate _A. m. capensis_ mutation rates

## Aims
- To familiarise self with the processes and filtering involved
- Produce a filtered .vcf file of `Larv09` for comparison with `Fdrone` 
- Produce full pipeline from .fastq to .vcf --> develop into multi-sample pipeline
	- I/O files
	- Resource usage

## Template
Following GATK 3.X best practices from Van der Auwera et al. (2013)

## Protocol 

### Download the reference genome for _Apis mellifera_:
```
cd /scratch/Scape/fred/
rsync --copy-links --times --verbose \
rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz \ 
/scratch/Scape/fred/ # Artemis HPC
```

Output: `GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz` 

### Indexing

#### 1. Index reference
```
qsub scripts/1_ref_idx.sh
mv /scratch/Scape/fred/GCF_* 1_ref_idx
```

Input: `GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz`

Outputs:
```
      859 Jun 22 16:04 GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz.amb
    27911 Jun 22 16:04 GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz.ann
225250988 Jun 22 16:04 GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz.bwt
 56312723 Jun 22 16:04 GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz.pac
112625496 Jun 22 16:05 GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz.sa
```

#### 2. Index fasta

Decompress `GCA_003254395.2_Amel_HAv3.1_genomic.fna.gz`:
```
cp GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz 2_fa_idx
bgzip -d 2_fa_idx/GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz
```

And index it
```
module load samtools/1.9
samtools faidx 2_fa_idx/GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz 
```

Outputs: `GCF_003254395.2_Amel_HAv3.1_genomic.fna.fai`

### 3. Create dictionary
```
qsub scripts/3_create_dict.sh
```

Input: `GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz`

Output: `GCF_003254395.2_Amel_HAv3.1_genomic.dict`

Note: no read group info in *.dict
- may not need to use downstream?
- bwa mem step below will be different to source (read group not specified)

### 4. BWA-MEM

Suggested use for Illumina paired-end reads longer than ~70:
```
qsub scripts/4_bwa_mem.sh
```

Note: `-M` option [not used](https://gatkforums.broadinstitute.org/gatk/discussion/21351/bwa-mem-m-option)

Inputs:
```
GCA_3254396sn.2_Amel_HAv3.1_genomic.fna.gz.bwt
GCA_003254395.2_Amel_HAv3.1_genomic.fna.gz.sa
GCA_003254395.2_Amel_HAv3.1_genomic.fna.gz.pac
GCA_003254395.2_Amel_HAv3.1_genomic.fna.gz.ann
GCA_003254395.2_Amel_HAv3.1_genomic.fna.gz.amb
```

Output: `Larv09_pe.sam`

### Sorting, converting to BAM, and marking dupes

#### 5. SortSam
```
qsub scripts/5_sortsam.sh
```

Input: `Larv09_pe_sorted.sam`

Output: `Larv09_pe_sorted.bam`

#### 6. MarkDuplicates
```
qsub scripts/6_markdups.sh
``` 

Input: `Larv09_pe_sorted.bam`

Outputs:
```
Larv09_pe_marked_dups.bam
marked_dup_metrics.txt
```

