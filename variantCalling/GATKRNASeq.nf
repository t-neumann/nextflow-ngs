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
    
    
snpEff	= "/groups/zuber/zubarchive/USERS/tobias/bin/snpEff/snpEff.jar"
snpSift	= "/biosw/generic-x86_64/snpeff/SnpSift.jar"
    
if (params.genome == "human") {
	snpEffConf  = "/groups/zuber/zubarchive/USERS/tobias/bin/snpEff/snpEff.config"
	genome 		= "/groups/zuber/zubarchive/USERS/tobias/hg38/GATK/Homo_sapiens_assembly38.fasta"
	genomeBuild = "hg38"
	knownIndels = ["/groups/zuber/zubarchive/USERS/tobias/hg38/GATK/Mills_and_1000G_gold_standard.indels.hg38.vcf",
					"/groups/zuber/zubarchive/USERS/tobias/hg38/GATK/Homo_sapiens_assembly38.known_indels.vcf"
					]
	snpEffIndels = "/groups/zuber/zubarchive/USERS/tobias/hg38/GATK/Homo_sampiesn_assembly38.known_indels_Mills_and_100G_gold_standard.merged.hg38.vcf"
	snpEffSnps 	= "/groups/zuber/zubarchive/USERS/tobias/hg38/GATK/1000G_hapmap_Axiom_Exome.merged.hg38.vcf"
	knownSnps 	= ["/groups/zuber/zubarchive/USERS/tobias/hg38/GATK/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
				"/groups/zuber/zubarchive/USERS/tobias/hg38/GATK/1000G_omni2.5.hg38.vcf.gz",
				"/groups/zuber/zubarchive/USERS/tobias/hg38/GATK/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz",
				"/groups/zuber/zubarchive/USERS/tobias/hg38/GATK/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf",
				"/groups/zuber/zubarchive/USERS/tobias/hg38/GATK/hapmap_3.3.hg38.vcf.gz"
				]
	dbSnp 		= "/groups/zuber/zubarchive/USERS/tobias/hg38/GATK/Homo_sapiens_assembly38.dbsnp138.vcf"
} else {
	snpEffConf  = "/groups/zuber/zubarchive/USERS/tobias/bin/snpEff/snpEff.config"
	genome		= "/groups/zuber/zubarchive/USERS/tobias/mareike/rnaseq/GATK/common/Mus_musculus.GRCm38.dna.primary_assembly.chr.fa"
	genomeBuild = "GRCm38.75"
	knownIndels	= ["/groups/zuber/zubarchive/USERS/tobias/mareike/rnaseq/GATK/common/C57BL_6NJ.mgp.v5.indels.dbSNP142.normed.vcf"]
	knownSnps	= ["/groups/zuber/zubarchive/USERS/tobias/mareike/rnaseq/GATK/common/GRCm38_dbsnp142_C57BL_6NJ.mgp.v5.vcf"]
	dbSnp		= "/groups/zuber/zubarchive/USERS/tobias/mareike/rnaseq/GATK/common/GRCm38_dbsnp142.vcf"
	snpEffIndels = "/groups/zuber/zubarchive/USERS/tobias/mareike/rnaseq/GATK/common/C57BL_6NJ.mgp.v5.indels.dbSNP142.normed.vcf"
	snpEffSnps 	= "/groups/zuber/zubarchive/USERS/tobias/mareike/rnaseq/GATK/common/GRCm38_dbsnp142_C57BL_6NJ.mgp.v5.vcf"
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

	#java -jar /biosw/debian7-x86_64/picard-tools/2.6.0/picard.jar ValidateSamFile I=!{name}.ro.dedup.bam O=!{name}.ro.dedup.val.txt VALIDATION_STRINGENCY=STRICT IGNORE="MISSING_TAG_NM"

	samtools index !{name}.ro.dedup.bam
	
    '''
}

process preformat {

	tag { name }
	
	module params.gatk
	module params.samtools
     
    input:
    set val(name), file(bam) from preprocessChannel
    file(bai) from preprocessBaiChannel
     
    output:
    set val(name), "${name}.recal.bam" into preformatChannel
    file("${name}.recal.bam.bai") into preformatBaiChannel
    
    shell:
    
    known = knownIndels.collect{"-known $it"}.join(' ')
    snps = knownSnps.collect{"-knownSites $it"}.join(' ')
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
    !{known} \
    -o !{name}.target_intervals.list \
    -nt !{params.threads}

    java -jar /biosw/generic-x86_64/gatk/3.11/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R !{genome} \
    -I !{name}.snt.bam \
    -targetIntervals !{name}.target_intervals.list \
    !{known} \
    -maxReads 1000000 \
    -maxInMemory 1000000 \
    -LOD 0.4 \
    -o !{name}.realign.bam

    java -jar /biosw/generic-x86_64/gatk/3.11/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R !{genome} \
    -I !{name}.realign.bam \
    !{snps} \
    -nct !{params.threads} \
    -o !{name}.recal_data.table
    
    java -jar /biosw/generic-x86_64/gatk/3.11/GenomeAnalysisTK.jar \
    -T PrintReads \
    -R !{genome} \
    -I !{name}.realign.bam \
    -BQSR !{name}.recal_data.table \
    -nct !{params.threads} \
    -o !{name}.recal.bam
    
    samtools index !{name}.recal.bam
	
    '''
}

process haplotype {

	tag { name }
	
	module params.gatk
     
    input:
    set val(name), file(bam) from preformatChannel
    file(bai) from preformatBaiChannel
     
    output:
    set val(name), "${name}.flt.vcf" into haplotypeChannel
    
    shell:
    ''' 
    shopt -s expand_aliases
	
	java -jar /biosw/generic-x86_64/gatk/3.11/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -R !{genome} \
    -I !{bam} \
    -recoverDanglingHeads \
    -dontUseSoftClippedBases \
    -stand_call_conf 20.0 \
    -stand_emit_conf 20.0 \
    -o !{name}.vcf

    java -jar /biosw/generic-x86_64/gatk/3.11/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R !{genome} \
    -V !{name}.vcf \
    -window 35 \
    -cluster 3 \
    -filterName FS \
    -filter "FS > 30.0" \
    -filterName QD \
    -filter "QD < 2.0" \
    -o !{name}.flt.vcf
	
    '''
}


process annotate {

	tag { name }
	
	module params.gatk
     
    input:
    set val(name), file(vcf) from haplotypeChannel
     
    output:
    file("${name}.spf.txt") into spfChannel
    file("${name}.final.vcf") into vcfChannel
    file("${name}_snpEff_genes.txt") into genesChannel
    file("${name}_snpEff_summary.html") into reportChannel
    
    shell:
    
    ''' 
    shopt -s expand_aliases
	
	grep "^#" !{vcf} > !{name}.flt_passed.vcf
	grep -v "^#" !{vcf} | awk \'BEGIN {FS = "\\t"} ($7 == "PASS")\' | grep -v "DP=0;" >> !{name}.flt_passed.vcf
	
	java -jar !{snpSift} annotate -id !{dbSnp} !{name}.flt_passed.vcf > !{name}.annotated.dbsnp.vcf
	
	java -jar !{snpSift} annotateMem -info INDEL -name "indels.mgp." \
	!{snpEffIndels} !{name}.annotated.dbsnp.vcf > !{name}.annotated.mpg_indels.vcf
	
	sed -i \'s/ID=INDEL,/ID=indels.mgp.INDEL,/\' !{name}.annotated.mpg_indels.vcf
	
	java -jar !{snpSift} annotateMem -info DP -name "snps.mgp." !{snpEffSnps} \
	!{name}.annotated.mpg_indels.vcf > !{name}.annotated.mpg.vcf

	grep "^##" !{name}.annotated.mpg.vcf > !{name}.annotated.mpg.tmp.vcf
	echo "##INFO=<ID=snps.mgp.DP,Number=1,Type=Integer,Description=\\"Raw read depth for SNPs in mouse genome project\\">" >> !{name}.annotated.mpg.tmp.vcf
	grep -v "^##" !{name}.annotated.mpg.vcf >> !{name}.annotated.mpg.tmp.vcf
	
	mv !{name}.annotated.mpg.tmp.vcf !{name}.annotated.mpg.vcf
	
	java -jar !{snpEff} -c !{snpEffConf} -v -o gatk !{genomeBuild} !{name}.annotated.mpg.vcf > !{name}.snpeff.vcf
	
	mv snpEff_genes.txt !{name}_snpEff_genes.txt
	mv snpEff_summary.html !{name}_snpEff_summary.html
	
	java -jar /biosw/generic-x86_64/gatk/3.11/GenomeAnalysisTK.jar \
    -T VariantAnnotator \
    -R !{genome} \
    -A SnpEff \
    --variant !{name}.annotated.mpg.vcf \
    --snpEffFile !{name}.snpeff.vcf \
    -L !{name}.annotated.mpg.vcf \
    -o !{name}.final.vcf

	java -jar /biosw/generic-x86_64/gatk/3.11/GenomeAnalysisTK.jar \
    -T VariantsToTable \
    -R !{genome} \
    -V !{name}.final.vcf \
    -SMA -AMD -o !{name}.spf.txt \
    -F CHROM -F POS -F ID -F REF -F DP -F ALT -F QUAL -F AC -F AF -F AN \
    -F SNPEFF_AMINO_ACID_CHANGE -F SNPEFF_CODON_CHANGE \
    -F SNPEFF_EFFECT -F SNPEFF_EXON_ID -F SNPEFF_FUNCTIONAL_CLASS -F SNPEFF_GENE_BIOTYPE \
    -F SNPEFF_GENE_NAME -F SNPEFF_IMPACT -F SNPEFF_TRANSCRIPT_ID \
    -F indels.mgp.INDEL -F snps.mgp.DP \
    -F MQRankSum -F QD -F ReadPosRankSum \
    -GF GT -GF AD -GF DP -GF GQ -GF PL
	
    '''
}

 
 
workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
