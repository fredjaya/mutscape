#!/usr/bin/env nextflow

// Define global parameters 
params.scpath = "/scratch/Scape/fred/2008_manual" 
params.mode   = false // As processes are run individually, important to specify correct mode
params.ref    = "${params.scpath}/GCF_003254395.2_Amel_HAv3.1_genomic.fna" 
params.outdir = "${params.scpath}"

sampleId_ch = Channel
    .of('Larv01', 'Larv02','Larv03','Larv04','Larv05','Larv06','Larv07',
        'Larv08', 'Larv09','Larv10','Larv11','Larv12','Larv13','Larv14',
        'Larv01', 'Larv01', 'Fdrone', 'Worker')

if (params.mode == 'bwaMapReads') {
// Map trimmed reads (paired) to reference genome

    params.bwaidx = "${params.scpath}/GCF_003254395.2_Amel_HAv3.1_genomic.fna.*"
    params.reads  = "/project/Scape/Trimmomatic/Paired/*_{R1,R2}_paired.fastq.gz"
    
    log.info """\
    
    ====================
    ref    : ${params.ref}
    bwaidx : ${params.bwaidx}
    reads  : ${params.reads}
    outdir : ${params.outdir}
    ====================
    """
    
    reads_ch  = Channel.fromFilePairs(params.reads)
    bwaidx_ch = Channel.fromPath(params.bwaidx)
     
    process bwaMapReads {
        
        cpus = 4
        memory = 16.GB
        time = '3h'
        publishDir "${params.outdir}" 
        tag "$sampleId"
    
        input:
            path ref from params.ref
            path bwaidx from bwaidx_ch.collect()
            tuple val(sampleId), path(reads) from reads_ch
    
        output:
            tuple val(sampleId), path("${sampleId}.sam") into sam_ch
    
        script:
        readGroup = "@RG\\tID:${sampleId}\\tSM:${sampleId}\\tPL:ILLUMINA\\tLB:lib1\\tPU:u"
        """
        module load bwa/0.7.17
    
        bwa mem -R \"${readGroup}\" -t 24 ${ref} ${reads} > ${sampleId}.sam
        """
    }
}

if (params.mode == 'samtoolsSort') {
// Convert .sam into sorted .bam 

    params.sam = "${params.scpath}/*.sam"
    
    log.info """\
    
    ====================
    sam    : ${params.sam}
    outdir : ${params.outdir}
    ====================
    """
    
    sam_ch  = Channel.fromPath(params.sam)
     
    process samtoolsSort {
        
        cpus = 1
        memory = 20.GB
        time = '1h'
        publishDir "${params.outdir}" 
        tag "sam"
    
        input:
            path sam from sam_ch
            // Should've used dynamic input sampleId names, as samples are known
        output:
            path "*_sorted.bam"
    
        script:
        """
        module load samtools/1.9
    
        samtools sort \
            -o ${sam}_sorted.bam \
            -O bam \
            ${sam} 
        """
    }
}

if (params.mode == 'markDuplicates') {
// Mark duplicates 

    params.sortedbam = "${params.scpath}/*.sam_sorted.bam"
    
    log.info """\
    
    ====================
    sortedbam : ${params.sortedbam}
    outdir    : ${params.outdir}
    ====================
    """
    
    sortedbam_ch = Channel.fromPath(params.sortedbam)
     
    process markDuplicates {
        
        cpus = 1
        memory = 32.GB
        time = '1h'
        publishDir "${params.outdir}" 
        tag '${sampleId}'
    
        input:
            val sampleId from sampleId_ch
            path "${sampleId}.sam_sorted.bam" from sortedbam_ch
 
        output:
            path "${sampleId}_marked_dups.bam"
            path "${sampleId}_marked_dup_metrics.txt" 
            path "${sampleId}_marked_dups.bai"
 
        script:
        """
        module load picard/2.7.1 
    
        picard MarkDuplicates \
            I=${sampleId}.sam_sorted.bam \
            O=${sampleId}_marked_dups.bam \
            M=${sampleId}_marked_dup_metrics.txt \
            CREATE_INDEX=true
        """
    }
}

