# mutscape

Variant calling pipeline for estimating *A. m. capensis* mutation rates.

## Branch notes

This branch (gvcf) revises the variant calling step in the previous pipeline, follow [GATK4X best practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-)

Variant calling will be conducted using HaplotypeCaller in GVCF mode with GATK4.1.8.0 

BQSR and first pass variant calls for known sites will remain unchanged from the one using GATK3.8

## Downloading and indexing the reference genome

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

## Preparing sequence data (.fastq to .bam)

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

## Base Quality Score Recalibration
Bootstrapping known sites according to [documentation](https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-methods-and-algorithms/Base_Quality_Score_Recalibration_(BQSR).md). 

```
# Covariate tables and plots not generated here (TODO).

# 6. Call confident sites for BSQR
qsub run_scripts/6_haplotypeCallerUnrecal

# 7. Generate covariate tables for BQSR
qsub run_scripts/7_bsqrTable

# 8. Recalibrate base quality scores
qsub run_scripts/8_recalibrateBQS.sh
```

## Variant calling

```
# 9. Run HaplotypeCaller in GVCF mode (single-sample) 
qsub run_scripts/9_gvcf.sh

# 10. Consolidate GVCFs
runscripts/10_consolidateGVCFs.sh

# 11. Joint genotyping
runscripts/11_jointGenotyping.sh
```

## Variant filtering
```
export DIR=/home/meep/Desktop/People/fred/Dropbox/meep/bee/02_working/2010_consolidate
export GATK=/home/meep/Desktop/Biocomputing/gatk-4.1.8.1/gatk
export REF=/home/meep/Desktop/People/fred/Dropbox/meep/bee/02_working/2009_filter_vcf/GCF_003254395.2_Amel_HAv3.1_genomic.fa

# Get initial variant counts
bin/countVariants.sh $DIR/joint_genotype.vcf.gz $DIR/joint_genotype.stats

# SNP sites: 3156296
# Indel sites: 1520510

# Output SNPs 
$GATK SelectVariants \
	-R $REF \
	-V $DIR/joint_genotype.vcf.gz \
	--select-type-to-include SNP \
	-O $DIR/joint_genotype_SNPs.vcf.gz

# Count sites
bin/countVariants.sh $DIR/joint_genotype_SNPs.vcf.gz $DIR/joint_genotype.stats	

# Sites: 3114386

# Variants to table for diagnostic plots
$GATK VariantsToTable \
	-R $REF \
	-V $DIR/joint_genotype_SNPs.vcf.gz \
	-F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR -GF GQ  \	
	-O $DIR/joint_genotype_SNPs.table

# Count SNPs per sample
bin/countVariants.sh $DIR/joint_genotype_SNPs.vcf.gz $DIR/joint_genotype.stats
plot-vcfstats -p $DIR $DIR/joint_genotype.stats
bin/get_sample_SNPs.sh $DIR/joint_genotype.stats $DIR/joint_genotype.psc

# Ran diagnosis.R here

# Site filters according to [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants)

$GATK VariantFiltration \
	-R $REF \
	-V $DIR/joint_genotype_SNPs.vcf.gz \
	-O $DIR/sitesFiltered_SNPs.vcf.gz \
	--filter-name "site_filter_QD" \
	--filter-expression "QD < 2.0" \
	--filter-name "site_filter_FS" \
	--filter-expression "FS > 60.0" \
	--filter-name "site_filter_SOR" \
	--filter-expression "SOR > 5.0" \
	--filter-name "site_filter_MQ" \
	--filter-expression "MQ < 40.0" \
	--filter-name "site_filter_MQRankSum" \
	--filter-expression "MQRankSum < -12.5" \
	--filter-name "site_filter_ReadPosRankSum" \
	--filter-expression "ReadPosRankSum > 10.0"

# Separate PASS and FAIL sites
gzip -d $DIR/sitesFiltered_SNPs.vcf.gz
grep -E '^#|PASS' $DIR/sitesFiltered_SNPs.vcf > $DIR/sitesFiltered_PASS.vcf

grep -E '^#' $DIR/sitesFiltered_SNPs.vcf > $DIR/sitesFiltered_FAIL.vcf && \
grep -Ev 'PASS' $DIR/sitesFiltered_SNPs.vcf >> $DIR/sitesFiltered_FAIL.vcf

bin/countVariants.sh $DIR/sitesFiltered_PASS

# Sites: 2984762
# Removed 1129624

## Genotype filtering DP < 5 and GQ < 30 (Yang, 2017; Yagound, 2020)

# Create table
$GATK VariantsToTable \
	-R $REF \
	-V $DIR/sitesFiltered_PASS.vcf \
	-F CHROM -F POS -GF GT -GF DP -GF GQ \
	-O $DIR/sitesFiltered_PASS.table

# Generate diagnostic plots for DP and GQ

# Genotype filtering
# DP > 80 removed, likely mapping error based on highest average SNP depth form joint_genotype.vcf
$GATK VariantFiltration \
	-R $REF \
	-V $DIR/sitesFiltered_PASS.vcf \
	-O $DIR/sitesFiltered_PASS_gtFiltered.vcf \
	--genotype-filter-expression "GQ < 30" \
	--genotype-filter-name "GQ_30" \
	--genotype-filter-expression "DP < 5 || DP > 80" \
	--genotype-filter-name "DP_5-80"

# Mask SNPs that failed genotype filters
$GATK SelectVariants \
	-R $REF \ 
	-V $DIR/sitesFiltered_PASS_gtFiltered.vcf \
	-O $DIR/sitesPASS_gtPASS.vcf \
	--set-filtered-gt-to-nocall

# Get counts
bin/countVariants.sh $DIR/sitesPASS_gtPASS 

# Retain sites with one unique genotype
# This removes sites with identical genotypes across samples
# And sites where the same mutation is found in >1 sample
python3 bin/genofreq.py genofreq $DIR/sitesPASS_gtPASS.vcf $DIR/genofreq.vcf
bin/countVariants.sh $DIR/genofreq

# Sites: 220902
# Removed 2763860

# Remove sites without a SNP called in Worker
python3 bin/genofreq.py orphan $DIR/genofreq.vcf $DIR/genofreq_orphan.vcf 
bin/countVariants.sh $DIR/genofreq_orphan

# Sites: 160830
# Removed 60072

# Filter 
```

