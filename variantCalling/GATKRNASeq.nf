#!/usr/bin/env nextflow

/*
* MIT License
*
* Copyright (c) 2017 Tobias Neumann
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/

log.info " GATK basic variant calling pipeline on RNA-Seq data "
log.info "====================================================="
log.info "bams              : ${params.bams}"
log.info "genome            : ${params.genome}"
log.info "threads           : ${params.threads}"

Channel
    .fromFilePairs( params.bams, size: 1 )
    .ifEmpty { error "Cannot find any reads matching: ${params.bams}" }
    .set { bams }
    
if (params.genome == "human") {
	genome = ""
} else {
	genome		= "/groups/zuber/zubarchive/USERS/tobias/mareike/rnaseq/GATK/common/Mus_musculus.GRCm38.dna.primary_assembly.chr.fa"
	knownIndels	= "/groups/zuber/zubarchive/USERS/tobias/mareike/rnaseq/GATK/common/C57BL_6NJ.mgp.v5.indels.dbSNP142.normed.vcf"
	dbsnp		= "/groups/zuber/zubarchive/USERS/tobias/mareike/rnaseq/GATK/common/GRCm38_dbsnp142_C57BL_6NJ.mgp.v5.vcf"
}

process preprocess {

	tag { name }
	
	module params.samtools
	module params.picard
     
    input:
    set val(name), file(bam) from bams
     
    output:
    set val(name), "${name}.ro.dedup.bam" into preprocessChannel
    file("${name}.ro.dedup.bam.bai") into preprocessBaiChannel
    
    shell:
    ''' 
    shopt -s expand_aliases
	
	java -jar /biosw/debian7-x86_64/picard-tools/2.6.0/picard.jar AddOrReplaceReadGroups INPUT=!{bam} OUTPUT=!{name}.rg.bam RGID=!{name} RGLB=!{name} RGPL="Illumina" RGPU=!{name} RGSM=!{name}

	java -jar /biosw/debian7-x86_64/picard-tools/2.6.0/picard.jar ReorderSam INPUT=!{name}.rg.bam OUTPUT=!{name}.ro.bam REFERENCE=!{genome} ALLOW_CONTIG_LENGTH_DISCORDANCE=false

	samtools sort -@ 10 -o !{name}.ro.sorted.bam !{name}.ro.bam

	java -jar /biosw/debian7-x86_64/picard-tools/2.6.0/picard.jar MarkDuplicates INPUT=!{name}.ro.sorted.bam OUTPUT=!{name}.ro.dedup.bam METRICS_FILE=!{name}.ro.dedup.metrics.txt

	java -jar /biosw/debian7-x86_64/picard-tools/2.6.0/picard.jar ValidateSamFile I=!{name}.ro.dedup.bam O=!{name}.ro.dedup.val.txt VALIDATION_STRINGENCY=STRICT IGNORE="MISSING_TAG_NM"

	samtools index !{name}.ro.dedup.bam
	
    '''
}

process preformat {

	tag { name }
	
	module params.gatk
     
    input:
    set val(name), file(bam) from preprocessChannel
    file(bai) from preprocessBaiChannel
     
    output:
    set val(name), "${name}.recal.bam" into preformatChannel
    
    shell:
    ''' 
    shopt -s expand_aliases
	
	java -jar /biosw/generic-x86_64/gatk/3.11/GenomeAnalysisTK.jar \
    -T SplitNCigarReads \
    -o !{name}.snt.bam \
    -I !{bam} \
    -R !{genome} \
    -rf ReassignOneMappingQuality \
    -RMQF 255 \
    -RMQT 60 \
    -U ALLOW_N_CIGAR_READS

    java -jar /biosw/generic-x86_64/gatk/3.11/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R !{genome} \
    -I !{name}.snt.bam \
    -known !{knownIndels} \
    -o !{name}.target_intervals.list \
    -nt !{params.threads}

    java -Xmx4g -jar /biosw/generic-x86_64/gatk/3.11/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R !{genome} \
    -I !{name}.snt.bam \
    -targetIntervals !{name}.target_intervals.list \
    -known !{knownIndels} \
    -maxReads 1000000 \
    -maxInMemory 1000000 \
    -LOD 0.4 \
    -o !{name}.realign.bam

    java -jar /biosw/generic-x86_64/gatk/3.11/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R !{genome} \
    -I !{name}.realign.bam \
    -knownSites !{dbsnp} \
    -nct !{params.threads} \
    -o !{name}.recal_data.table
    
    java -jar /biosw/generic-x86_64/gatk/3.11/GenomeAnalysisTK.jar \
    -T PrintReads \
    -R !{genome} \
    -I !{name}.realign.bam \
    -BQSR !{name}.recal_data.table \
    -nct !{params.threads} \
    -o !{name}.recal.bam
	
    '''
}


 
 
workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
