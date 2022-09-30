#!/bin/bash
#SBATCH --job-name=RNAseqBigWig
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --mail-type=end
#SBATCH --mail-user=megannj@email.unc.edu
#SBATCH --time=120:00:00
#SBATCH -e RNAseqBigWig.%job.err

module load samtools
module load deeptools/2.4.1

##Here are the variables you need to change before running the script. PINE_DIR should be the path to your pine scratch space. FASTQ_FILE should be the unzipped fastq file. FILE_NAME should be an identifier used for the data. The output name from HTSF can be long, so pick something simpler, such as the name of the antibody. 
FASTQ_DIR=/proj/dowenlab/users/megan/WizKO-RNAseq/*.bam

for f in ${FASTQ_DIR}
	do

##Now we need to do a bowtie alignment of the fastq file. -a designates all. -v is the alignment mode. -p is number of threads -S indicates input is in sam. -M indicates number of duplicate alignments reported. --best and --strata ensure you're reporting the best alignment. 
	#bowtie -a -v 2 -p 24 -S -M 1 --best --strata /proj/seq/data/MM10_UCSC/Sequence/BowtieIndex/genome ${FASTQ_FILE}.fastq > ${FASTQ_FILE}.sam

##The sam file must be converted to a sorted bam file. 
	#samtools view -S -b -@ 24 -o ${FASTQ_FILE}.bam ${FASTQ_FILE}.sam
	#samtools sort -@ 24 -o ${FASTQ_FILE}_sorted.bam ${FASTQ_FILE}.bam

##Next, remove PCR duplicates.			
	#samtools sort -n -@ 24 ${f}_sorted.bam -o ${f}_sortedbyname.bam
	#samtools fixmate -m -@ 24 ${FA}_sortedbyname.bam ${FASTQ_FILE}_fixmate.bam
	#samtools sort -@ 24 ${FASTQ_FILE}_fixmate.bam -o ${FASTQ_FILE}_fixmate_sorted.bam
	#samtools markdup -r -s -@ 24 ${FASTQ_FILE}_fixmate_sorted.bam ${FASTQ_FILE}_markdup.bam

##Finally, index the bam file and use deeptools to create a bigwig file for uploading into UCSC.
#module load deeptools
	#samtools index ${f} ${f}.bai	
	bamCoverage -p 24 -bs 10 --normalizeUsingRPKM -ignore chrX -b ${f} -o ${f}_coverage.bw
	done
