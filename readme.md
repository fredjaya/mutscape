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

#### 3. Create dictionary
```
qsub scripts/3_create_dict.sh
```

Input: `GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz`

Output: `GCF_003254395.2_Amel_HAv3.1_genomic.dict`

Note: no read group info in *.dict
- read groups are added in step 7, after duplicates are marked

#### 4. BWA-MEM

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

Input: `Larv09_pe.sam`

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

### Local realignment around indels

#### 7. Add read groups 
```
qsub scripts/7_add_rg.sh
```

Input: `Larv09_pe_marked_dups.bam`

Output: `Larv09_marked_dups_rg.bam`

#### 8. RealignerTargetCreator

`-L 20` flag omitted

```
qsub scripts/8_rtc.sh
```

```
##### ERROR MESSAGE: Invalid command line: Cannot process the provided BAM/CRAM fi
le(s) because they were not indexed.  The GATK does offer limited processing of un
indexed BAM/CRAMs in --unsafe mode, but this feature is unsupported -- use it at y
our own risk!                                                                     
```
- symlinks dont work?
- are real files needed to be in the same dir, with no option to specify path of .idx?
- try:

```
# Copy reference, index, and dictionary into the same dir
cd /scratch/Scape/fred
mkdir 8_rtc_cp
cp 2_fai_dx/GCF_003254395.2_Amel_HAv3.1_genomic.fna 8_rtc_cp
cp 2_fai_dx/GCF_003254395.2_Amel_HAv3.1_genomic.fna.fai 8_rtc_cp
cp 3_create_dict/GCF_003254395.2_Amel_HAv3.1_genomic.dict 8_rtc_cp  
cp 7_add_rg/Larv09_marked_dups_rg.bam 8_rtc_cp

# Run RealignTargetCreator
cd ~/mutscape
qsub scripts/8_rtc_cp.sh
```
- nope, same error ):
- try rerunning all the index steps in one dir?

```
# Copy ref into subdir, as output files automatically write to same dir
cd /scratch/Scape/fred
cp GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz 8_rtc_idx

# Index reference
cd ~/mutscape
qsub scripts/8_ref_idx.sh

# Decompress ref.gz, as samtools cannot read this compression format
bgzip -d < GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz > GCF_003254395.2_Amel_HAv3.1_genomic.fna

# And index fasta
module load samtools/1.9
samtools faidx GCF_003254395.2_Amel_HAv3.1_genomic.fna 

# Create dictionary
qsub scripts/8_create_dict.sh

# 
qsub scripts/8_rtc_singledir.sh
```
Inputs: 
```
GCF_003254395.2_Amel_HAv3.1_genomic.fna
GCF_003254395.2_Amel_HAv3.1_genomic.dict
Larv09_marked_dups_rg.bam
# Possibly more of the output files when indexing?
```

According to [this thread](https://gatkforums.broadinstitute.org/gatk/discussion/2617/bam-is-not-indexed-or-bad-input) it may be an issue of needing to `create_index=true` during markdupes

Going back and adding `CREATE_INDEX=true` when marking duplicates
