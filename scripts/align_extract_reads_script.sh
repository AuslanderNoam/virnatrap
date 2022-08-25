#!/bin/bash -l
#SBATCH -o out.txt
#SBATCH --job-name=alext
#SBATCH -e error.txt
#SBATCH  -n 15
         


#workdir line below has to be set to your working directory
#bowtie2Indexhg19 line below has to be set to where the hg19 index is
#bowtie2Indexphage line below has to be set to where the phix index is
######Change following as needed######
workdir=**path-do-work-directory**
bowtie2Indexhg19=**path-do-hg19-index**
bowtie2Indexphage=**path-do-phix-index**
alignHg19=Yes


CWD=$(pwd)
rst=$CWD/Results
mkdir $rst
output=$CWD/tmp_manifest.txt
firstout=$CWD/unmapped_reads#.fastq
firstout1=$CWD/unmapped_reads_1.fastq
firstout2=$CWD/unmapped_reads_2.fastq
firstout3=$CWD/unmapped_reads_M.fastq
bowtieout1=$CWD/bowtie1out.fastq
bowtieout1_1=$CWD/bowtie1out.1.fastq
bowtieout1_2=$CWD/bowtie1out.2.fastq

mapped=$CWD/Mapped.sam

for infile in $workdir/*/*.bam
	do name=$(basename $infile .bam) 
	if [ "$alignHg19" == "yes" ]; then
	{
   
		bam2fastq -o $firstout -f --no-aligned --unaligned --no-filtered $infile

		if [ -s $firstout1 ]; then
			bowtie2 -p 32 -N 1  --local -X 1500 --no-unal --un-conc $bowtieout1 -x $bowtie2Indexhg19 -1 $firstout1 -2 $firstout2 -S $mapped
			bowtie2 -p 32 -N 1  --local -X 1500 --no-unal --un-conc $rst/${name}_unmapped.fastq -x $bowtie2Indexphage -1 $bowtieout1_1 -2 $bowtieout1_2 -S $mapped
		else 
			bowtie2 -p 32 -N 1  --local --no-unal --un $bowtieout1 -x $bowtie2Indexhg19 $firstout3 -S $mapped
			bowtie2 -p 32 -N 1  --local --no-unal --un $rst/${name}_unmapped.fastq -x $bowtie2Indexphage $bowtieout1 -S $mapped
		fi
  
	} 
	else

    	{

     bam2fastq -o $firstout -f --no-aligned --unaligned --no-filtered $infile

     if [ -s $firstout1 ]; then
       
        bowtie2 -p 32 -N 1  --local -X 1500 --no-unal --un-conc $rst/${name}_unmapped.fastq -x $bowtie2Indexphage -1 $firstout1 -2 $firstout2 -S $mapped
     else 
        
        bowtie2 -p 32 -N 1  --local --no-unal --un $rst/${name}_unmapped.fastq -x $bowtie2Indexphage $firstout3 -S $mapped
     fi
  
 
   } 

fi
done
