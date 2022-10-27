## 1. Pre-processing

A Requirements

(a) bowtie2: https://github.com/BenLangmead/bowtie2

(b) bam2fastq: https://github.com/jts/bam2fastq

(c) gdc-client: https://github.com/topics/gdc-client

B Instructions
We provide an example manifest for download with 10 LIHC samples (manifest_LIHC_example.txt), where the token file path is defined as 'path-to-token'

B.1 To obtain the TCGA BAM files given dbGAP access, please follow the GDC instruction using gdc-client:

>gdc-client download -m manifest_LIHC_example.txt -t path-to-token

B.2 Use bowtie to align against hg19 and the phix phage and get the unmapped reads. For that use the align_extract_reads.sh script . Set the paths to  bowtie2 and bam2fastq, and change the following 3 lines in the align_extract_reads.sh script to define: (a) workdir - the working directory, where bam files are located. (b) bowtie2Indexhg19 - path to the hg19 index, and (c) bowtie2Indexphage - path to the phix index. 

Then, running the align_extract_reads.sh script, it will create in workdir a directory called Results, with the paired-end files

B.2 Use the cat_paired.py scripts to copy paired-end files into one fastq file of all reads


## 2. Post-processing

A Requirements

(a) Blast: https://www.ncbi.nlm.nih.gov/books/NBK569861/

(b) STAR: https://github.com/alexdobin/STAR

divergent virus database (C-RVDBv16.0.txt): https://herv.img.cas.cz

Use run_blast_os_par.py script for post-processing and running Blast for the resulting contigs

