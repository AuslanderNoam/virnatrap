## 1. Pre-processing

A Requirements

(a) bowtie2: https://github.com/BenLangmead/bowtie2

(b) bam2fastq: https://github.com/jts/bam2fastq

(c) gdc-client: https://github.com/topics/gdc-client

B Instructions

B.1 To obtain the TCGA BAM files given dbGAP access, please follow the GDC instruction using gdc-client:

>gdc-client download -m manifest_file.txt -t path-to-token

B.2 Use the align_extract_reads.sh script to extract unaligned reads from the BAM files to fastq files

B.2 Use the cat_paired.py scripts to copy paired-end files into one fastq file of all reads


## 2. Post-processing

A Requirements

(a) Blast: https://www.ncbi.nlm.nih.gov/books/NBK569861/

(b) STAR: https://github.com/alexdobin/STAR

divergent virus database (C-RVDBv16.0.txt): https://herv.img.cas.cz

Use run_blast_os_par.py script for post-processing and running Blast for the resulting contigs

