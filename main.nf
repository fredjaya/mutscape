#!/usr/bin/env nextflow

// Define global parameters 
params.scpath       = "/scratch/Scape/fred/2008_manual" 
params.mode         = false // As processes are run individually, important to specify correct mode
params.ref          = "${params.scpath}/GCF_003254395.2_Amel_HAv3.1_genomic.fna" 
params.fai          = "${params.ref}.fai"
params.dict         = "${params.scpath}/GCF_003254395.2_Amel_HAv3.1_genomic.dict"
params.markedbam    = "${params.scpath}/*.sam_sorted.bam_marked_dups.{bam,bai}"
params.realignedbam = "${params.scpath}/*_realigned.{bam,bai}"
params.outdir       = "${params.scpath}"
    
markedbam_ch    = Channel.fromFilePairs(params.markedbam)
realignedbam_ch = Channel.fromFilePairs(params.realignedbam)
 
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
    sam_ch  = Channel.fromPath(params.sam)
    
    log.info """\
    
    ====================
    sam    : ${params.sam}
    outdir : ${params.outdir}
    ====================
    """
     
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
    sortedbam_ch = Channel.fromPath(params.sortedbam)
    
    log.info """\
    
    ====================
    sortedbam : ${params.sortedbam}
    outdir    : ${params.outdir}
    ====================
    """
    
    process markDuplicates {
        
        cpus = 1
        memory = 32.GB
        time = '2h'
        publishDir "${params.outdir}" 
        tag "${bam}"
    
        input:
            path bam from sortedbam_ch

        output:
            path "*_marked_dups.bam"
            path "*_marked_dups_metrics.txt" 
            path "*_marked_dups.bai"
 
        script:
        """
        module load picard/2.7.1 
    
        picard MarkDuplicates \
            I=${bam} \
            O=${bam}_marked_dups.bam \
            M=${bam}_marked_dups_metrics.txt \
            CREATE_INDEX=true
        """
    }
}

if (params.mode == 'realignerTargetCreator') {
// Mark targets to be realigned around indels
    
    log.info """\
    
    ====================
    ref       : ${params.ref}
    fai       : ${params.fai}
    dict      : ${params.dict}
    markedbam : ${params.markedbam}
    outdir    : ${params.outdir}
    ====================
    """
    
    process realignerTargetCreator {
        
        cpus = 1
        memory = 8.GB
        time = '1h'
        publishDir "${params.outdir}" 
        tag "$sampleId"
    
        input:
            path ref  from params.ref
            path fai  from params.fai 
            path dict from params.dict
            tuple val(sampleId), path(bamfiles) from markedbam_ch

        output:
            tuple val(sampleId), path("${sampleId}_target_intervals.list")
 
        script:
        // https://github.com/nextflow-io/nextflow/issues/1246#issuecomment-515919438 
        def bam = bamfiles.findAll{ it.toString() =~ /.bam$/ }.join('')        
        """
        module load gatk/3.8.1

        gatk -T RealignerTargetCreator \
             -R ${ref} \
             -I ${bam} \
             -o ${sampleId}_target_intervals.list        
        """
    }
}

if (params.mode == 'indelRealigner') {
// Realign regions around indels
    
    params.intervals = "${params.scpath}/*_target_intervals.list"                   
    intervals_ch = Channel                                                          
                    .fromPath(params.intervals)                                     
                    .map { file ->                                                  
                            def sampleId = file.name.toString().tokenize('_').get(0)
                            return tuple(sampleId, file)                            
                    }                                                               
                                                                                    
    realigntargets_ch = markedbam_ch.join(intervals_ch)                             
     
    log.info """\
    
    ====================
    ref       : ${params.ref}
    fai       : ${params.fai}
    dict      : ${params.dict}
    markedbam : ${params.markedbam}
    intervals : ${params.intervals}
    outdir    : ${params.outdir}
    ====================
    """
                                                                                    
    process indelRealigner {

        cpus = 2
        memory = 14.GB
        time = '5h'
        publishDir "${params.outdir}"
        tag "$sampleId"

        input:                                                                      
            path ref  from params.ref                                               
            path fai  from params.fai                                               
            path dict from params.dict                                              
            tuple val(sampleId), path(bamfiles), path(intervals) from realigntargets_ch
        
        output:
            tuple val(sampleId), path("${sampleId}_realigned.bam") 
            tuple val(sampleId), path("${sampleId}_realigned.bai") 

        script:                                                                     
        def bam = bamfiles.findAll{ it.toString() =~ /.bam$/ }.join('')             
        """                                                                         
        module load gatk/3.8.1                                                      
                                                                                    
        gatk -T IndelRealigner \
             -R ${ref} \
             -I ${bam} \
             -targetIntervals ${intervals} \
             -o ${sampleId}_realigned.bam 
        """                                                                         
    }                                                                               
 
}

if (params.mode == 'haplotypeCallerUnrecal') {
// Call confident sites for recalibration 
     
    log.info """\
    
    ====================
    ref          : ${params.ref}
    fai          : ${params.fai}
    dict         : ${params.dict}
    realignedbam : ${params.realignedbam}
    outdir       : ${params.outdir}
    ====================
    """
                                                                                    
    process haplotypeCallerUnrecal {                                            
    
        cpus = 4                                                           
        memory = 16.GB                                                     
        time = '14h'                                                       
        publishDir "$params.outdir"                                        
        tag "$sampleId"                                                    
                                                                           
        input:                                                             
            path ref  from params.ref                                      
            path fai  from params.fai                                      
            path dict from params.dict                                     
            tuple val(sampleId), path(bamfiles) from realignedbam_ch       
                                                                           
        output:                                                            
            tuple val(sampleId), path("${sampleId}_unrecal_conf_sites.vcf")
                                                                           
        script:                                                            
        def bam = bamfiles.findAll{ it.toString() =~ /.bam$/ }.join('')    
        """                                                                
        module load gatk/3.8.1                                             
                                                                           
        gatk -T HaplotypeCaller \
             -R ${ref} \
             -I ${bam} \
             --genotyping_mode DISCOVERY \
             --output_mode EMIT_ALL_CONFIDENT_SITES \
             -stand_call_conf 30 \
             -o ${sampleId}_unrecal_conf_sites.vcf \
             -nt 1 -nct 24                                                 
        """                                                                
         
    }                                                                               
 
}

if (params.mode == 'bqsrTable') {
// Generate table of covariance for bqsr

    params.knownsites = "${params.scpath}/*_unrecal_conf_sites.vcf"    
    knownsites_ch = Channel
                    .fromPath(params.knownsites)
                    .map { file ->
                            def sampleId = file.name.toString().tokenize('_').get(0)
                            return tuple(sampleId, file)
                    }
    bsqr_ch = realignedbam_ch.join(knownsites_ch)                             

    log.info """\
    
    ====================
    ref          : ${params.ref}
    fai          : ${params.fai}
    dict         : ${params.dict}
    realignedbam : ${params.realignedbam}
    knownsites   : ${params.knownsites}
    outdir       : ${params.outdir}
    ====================
    """
                                                                                    
    process bqsrTable {                                            
    
        cpus = 2
        memory = 12.GB                                                     
        time = '4h'                                                       
        publishDir "$params.outdir"                                        
        tag "$sampleId"                                                    
                                                                           
        input:                                                             
            path ref  from params.ref                                      
            path fai  from params.fai                                      
            path dict from params.dict                                     
            tuple val(sampleId), path(bamfiles), path(knownsites) from bsqr_ch      
                                                                        
        output:                                                            
            tuple val(sampleId), path("${sampleId}_unrecal_conf_sites.vcf.idx")
            tuple val(sampleId), path("${sampleId}_pre_recal.table")
                                                                           
        script:                                                            
        def bam = bamfiles.findAll{ it.toString() =~ /.bam$/ }.join('')    
        """                                                                
        module load gatk/3.8.1                                             
                                                                           
        gatk -T BaseRecalibrator \
             -R ${ref} \
             -I ${bam} \
             -knownSites ${knownsites} \
             -o ${sampleId}_pre_recal.table \
             -nt 1 -nct 24
        """                                                                
         
    }

}

if (params.mode == 'recalibrateBQS') {
// Recalibrate base quality scores 

    params.table = "${params.scpath}/*_pre_recal.table"
    table_ch = Channel
                .fromPath(params.table)
                .map { file ->
                        def sampleId = file.name.toString().tokenize('_').get(0)
                        return tuple(sampleId, file)
                    }
    
    recal_ch = realignedbam_ch.join(table_ch)

    log.info """\
    
    ====================
    ref          : ${params.ref}
    fai          : ${params.fai}
    dict         : ${params.dict}
    realignedbam : ${params.realignedbam}
    table        : ${params.table}
    outdir       : ${params.outdir}
    ====================
    """
                                                                                    
    process recalibrateBQS {                                            
    
        cpus = 3
        memory = 18.GB                                                     
        time = '3h'                                                       
        publishDir "$params.outdir"                                        
        tag "$sampleId"                                                    
                                                                           
        input:                                                             
            path ref  from params.ref                                      
            path fai  from params.fai                                      
            path dict from params.dict                                     
            tuple val(sampleId), path(bamfiles), path(table) from recal_ch      
                                                                        
        output:                                                            
            tuple val(sampleId), path("${sampleId}_recalibrated_reads.bam")
            tuple val(sampleId), path("${sampleId}_recalibrated_reads.bai")
                                                                           
        script:                                                            
        def bam = bamfiles.findAll{ it.toString() =~ /.bam$/ }.join('')    
        """                                                                
        module load gatk/3.8.1                                             
                                                                           
        gatk -T PrintReads \
             -R ${ref} \
             -I ${bam} \
             -BQSR ${table} \
             -o ${sampleId}_recalibrated_reads.bam \
             -nt 1 -nct 3
        """                                                                
         
    }

}

