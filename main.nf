#!/usr/bin/env nextflow

// Define default parameters 
//// Left on false to avoid conflicts between processes
params.scpath = "/scratch/Scape/fred/2008_manual" 
params.mode   = false
params.ref    = "${params.scpath}/GCF_003254395.2_Amel_HAv3.1_genomic.fna" 
params.reads  = "/project/Scape/Trimmomatic/Paired/*_{R1,R2}_paired.fastq.gz"
params.outdir = "${params.scpath}"

// Map trimmed reads (paired) to reference genome
if (params.mode == 'bwaMapReads') {

log.info """\
====================
ref    : ${params.ref}
reads  : ${params.reads}
outdir : ${params.outdir}
====================
"""

reads_ch = Channel.fromFilePairs(params.reads)

process bwaMapReads {
    
    cpus = 4
    memory = 20.GB
    time = '3h'
    publishDir '${params.outdir}' 
    tag '$sampleId'

    input:
        path ref from params.ref
        tuple val(sampleId), path(reads) from reads_ch

    output:
        tuple val(sampleId), path("${sampleId}_pe_sorted.bam")

    script:
    readGroup = "@RG\\tID:${sampleId}\\tSM:${sampleId}\\tPL:ILLUMINA\\tLB:lib1\\tPU:u"
    """
    module load bwa/0.7.17

    bwa mem -R \"${readGroup}\" -t 24 ${ref} ${reads}
    """
}
