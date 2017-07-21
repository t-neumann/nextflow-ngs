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

log.info " VBCF jnomics star align "
log.info "========================="
log.info "Automatically detects    "
log.info "SE vs PE                 "
log.info "PE: _1, _2               "
log.info "SE: No \"1\" or \"2\"    "
log.info "SE: as last character!   "
log.info "                         "
log.info "readDir              : ${params.readDir}"

pairedEndRegex = params.readDir + "/*_{1,2}.fastq.gz"
SERegex = params.readDir + "/*[!12].fastq.gz"

pairFiles = Channel.fromFilePairs(pairedEndRegex)
singleFiles = Channel.fromFilePairs(SERegex, size: 1){ file -> file.baseName.replaceAll(/.fastq/,"") }

readsChannel  = singleFiles.mix(pairFiles)

process align {

	cpus = 1
	
	stageInMode = 'link'
	
	tag { name }
	 
    input:
    set val(name), file(reads) from readsChannel

    script:
    
    date = new java.util.Date().format( 'yyyyMMdd' )
    
    def single = reads instanceof Path
    
	if (!single)
    
        """
        module load cutadapt
        
        /groups/vbcf-ngs/bin/funcGen/jnomicss.sh alignStarAndDeploy \
                 --inputFile1 ${reads[0]}  \
                 --inputFile2 ${reads[1]} \
                 --indexFile /groups/vbcf-ngs/misc/genomes/human/ncbi38_hg20/star \
                 --cleanIndex /groups/vbcf-ngs/misc/genomes/contaminants/human/bowtie2/rrna \
                 --statGTF /groups/vbcf-ngs/misc/genomes/human/ncbi38_hg20/annotation/Homo_sapiens.GRCh38.78.gtf\
                 --geneIndex groups/vbcf-ngs/misc/genomes/human/ncbi38_hg20/annotation/Homo_sapiens.GRCh38.cds.all.fa \
                 --mappingPath "${name}/" \
                 --fastaRef /groups/vbcf-ngs/misc/genomes/human/ncbi38_hg20/fasta/hg20.fa \
                 --ensemblTranscriptsFasta /groups/vbcf-ngs/misc/genomes/human/ncbi38_hg20/annotation/Homo_sapiens.GRCh38.genes.fa \
                 --ensemblGTF /groups/vbcf-ngs/misc/genomes/human/ncbi38_hg20/annotation/Homo_sapiens.GRCh38.78.gtf \
                 --cutWith standard \
                 --htSeqStranded reverse \
                 --libraryType fr-firststrand \
                 --tags trans \
                 --skipregister \
                 --date ${date} \
                 --sampleIdName \
                 --sampleId ${task.index}
        """
    
    else 
    
        """
        module load cutadapt
        
        /groups/vbcf-ngs/bin/funcGen/jnomicss.sh alignStarAndDeploy \
                     --inputFile1 $reads  \
                     --indexFile /groups/vbcf-ngs/misc/genomes/human/ncbi38_hg20/star \
                     --cleanIndex /groups/vbcf-ngs/misc/genomes/contaminants/human/bowtie2/rrna \
                     --statGTF /groups/vbcf-ngs/misc/genomes/human/ncbi38_hg20/annotation/Homo_sapiens.GRCh38.78.gtf \
                     --geneIndex /groups/vbcf-ngs/misc/genomes/human/ncbi38_hg20/annotation/Homo_sapiens.GRCh38.cds.all.fa \
                     --mappingPath "${name}/" \
                     --fastaRef  /groups/vbcf-ngs/misc/genomes/human/ncbi38_hg20/fasta/hg20.fa \
                     --ensemblTranscriptsFasta /groups/vbcf-ngs/misc/genomes/human/ncbi38_hg20/annotation/Homo_sapiens.GRCh38.genes.fa \
                     --ensemblGTF  /groups/vbcf-ngs/misc/genomes/human/ncbi38_hg20/annotation/Homo_sapiens.GRCh38.78.gtf \
                     --cutWith standard \
                     --htSeqStranded reverse \
                     --libraryType fr-firststrand \
                     --tags trans \
                     --skipregister \
                     --date ${date} \
                     --sampleIdName \
                 	 --sampleId ${task.index}
        """
}

 
 
workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
