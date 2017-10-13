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

log.info " SWEMBL/SICER initiation site and zone peak calling pipeline "
log.info "============================================================="

masterChannel = Channel
         		.from(params.files)
         		.map { condition, list -> 
            		def files = list.collect{ file(it)} 
            		return tuple(condition, files) 
        		}
        		
masterChannel.into{sicerInputChannel; swemblInputChannel}

        		        		
SICERPath = "/groups/vbcf-ngs/bin/peaks/SICER_V1.1/SICER/SICER.sh"

process sicer {

	echo true

	tag { name }
	
	module params.bedtools
	module params.python
     
    input:
    set val(name), file(reads) from sicerInputChannel
     
    output:
    set val(name),  "${name}_SICER.bed" into sicerChannel
 
    shell:
    '''
    bamToBed -i !{reads[0]} > chip.bed
    bamToBed -i !{reads[1]} > control.bed
    
    !{SICERPath} $PWD chip.bed control.bed . mm9 1 200 600 0.74 600 0.01
    
    echo $\'#chrom\\tstart\\tend\\tname\\tfold_change\\tstrand\' > !{name}_SICER.bed
    
    sort -n -k 8 chip-W200-G600-islands-summary-FDR0.01 | awk \'BEGIN{FS="\\t"; OFS="\\t"} {print $1,$2,$3,NR,$7,"+"}\' >> !{name}_SICER.bed

    '''
}

process swembl {

	echo true

	tag { name }
	
	module params.samtools
     
    input:
    set val(name), file(reads) from swemblInputChannel
     
    output:
    set val(name), "${name}_SWEMBL.bed" into swemblChannel
 
    shell:
    '''
    samtools index !{reads[0]}
    
    samtools view -hb !{reads[0]} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
    chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
    | samtools sort -@ !{params.threads} -o !{name}_noContig.bam
    
    /groups/zuber/zubarchive/USERS/tobias/bin/SWEMBL/SWEMBL -F -i !{name}_noContig.bam -r !{reads[1]} -R 0.002 -o !{name}_SWEMBL_raw.bed
    
    grep -v "#" !{name}_SWEMBL_raw.bed | grep -v "Ref." | awk \'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "swembl_peak_"NR, $7, "+"}\' \\
    | sort -k1,1 -k2,2n > !{name}_SWEMBL.bed

    '''
}

process zone_site_calling {

echo true

	tag { name }
	
	module params.bedtools
     
    input:
    set val(name), file(peaks) from sicerChannel.phase(swemblChannel).map { left, right -> tuple(left[0], [left[1], right[1]]) }
     
    output:
    file "${name}_IS.bed" into outISChannel
    file "${name}_IZ.bed" into outIZChannel
 
    shell:
    '''
    bedtools intersect -a !{peaks[0]} -b !{peaks[1]} -u -wa > !{name}_IZ.bed
    bedtools intersect -b !{peaks[0]} -a !{peaks[1]} -u -wa > !{name}_IS.bed

    '''
}
 
workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
