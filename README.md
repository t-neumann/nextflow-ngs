Nextflow-NGS
============
[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.24.0-brightgreen.svg)](http://nextflow.io)

Generic pipelines for next-generation sequencing data analysis using [Nextflow](https://www.nextflow.io/). These pipelines comprise robust and reasonable tool selection for various
analysis steps.

File conversion
---------------
* **SRA conversion** - File conversion of \*.sra files from the [Read Sequencing Archive](https://www.ncbi.nlm.nih.gov/sra) to regular fastq-files using [sra-tools](https://github.com/ncbi/sra-tools).

Preprocessing
-------------
* **Adapter trimming** - Standard adapter trimming and QC using [trim_galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

Genome tools
------------
* **Mappability** - Genome-wide mappability assessment using the [GEM-library](http://algorithms.cnag.cat/wiki/The_GEM_library)

Mapping
-------
* **star_align** - [VBCF](http://www.vbcf.ac.at/facilities/next-generation-sequencing/) RNA-seq alignment and QA pipeline using [STAR](https://github.com/alexdobin/STAR)
* **ChIPSeq_align** - [VBCF](http://www.vbcf.ac.at/facilities/next-generation-sequencing/) ChIP-seq alignment and QA pipeline using [bowtie](http://bowtie-bio.sourceforge.net/index.shtml)

Variant calling
---------------
* **GATK RNA-Seq** - Variant calling on RNA-Seq data using Broad's Genome Analysis Toolkit following their [best practices](https://software.broadinstitute.org/gatk/documentation/article.php?id=3891).

Peak calling
------------
* **SICER** - Standard broad peak calling using [SICER](http://home.gwu.edu/~wpeng/Software.htm).
* **MACS2** - Standard general peak calling using [MACS2](https://github.com/taoliu/MACS).
* **SWEMBL** - Standard tight peak calling using [SWEMBL](http://www.ebi.ac.uk/~swilder/SWEMBL/).
* **SNSPeak** - Calling of initiation sites (IS) and initation zones (IZ) according to *[Cayrou et al, Genome Research 2015](http://genome.cshlp.org/content/25/12/1873)*.
* **SNSPeakMACS** - Calling of initiation sites (IS) and initation zones (IZ) using MACS2.


