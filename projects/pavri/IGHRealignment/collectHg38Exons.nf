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

log.info " Collect VDJ exons from hg38 Gencode annotation "
log.info "=================================="
log.info "VDJTab              : ${params.VDJTab}"
log.info "outDir               : ${params.outDir}"

vdjFile = file(params.VDJTab)

if( !vdjFile.exists() ) exit 1, "Missing VDJ file file: ${vdjFile}"

process collectVDJs {

	publishDir params.outDir

     
    input:
    file vdjFile
     
    output:
    file "hg38_vdj_exons.bed" into outChannel
 
    """
    grep IGHV $vdjFile > hg38_vdj_exons.tmp
    grep IGHD $vdjFile | grep -v \$'IGHD\t' >> hg38_vdj_exons.tmp
    grep IGHJ $vdjFile >> hg38_vdj_exons.tmp
    awk 'BEGIN{FS="\t"; OFS="\t"} {print \$2, \$4, \$5, \$11, 0, \$3}' hg38_vdj_exons.tmp > hg38_vdj_exons.bed
    rm hg38_vdj_exons.tmp
    """
}
 
 
workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
