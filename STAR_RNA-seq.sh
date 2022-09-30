#!/bin/bash
#SBATCH --job-name=STARseq-pipeline
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --mail-type=end
#SBATCH --mail-user=narruda@email.unc.edu
#SBATCH --time=120:00:00
#SBATCH -e STARSeqPipeline.%j.err

#align RNA-seq data to produce genecount files to take into DESeq - need to make a combined excel file of all genecount files (WT_rep1, WT_rep2, WT_rep3 etc)

module load star
module load samtools

star --runThreadN 24 --genomeDir /proj/seq/data/STAR_genomes/mm10/ --readFilesIn /path/to/fastq/JMDc3_R1.fastq.gz /path/to/fastq/JMDc3_R2.fastq.gz --outFileNamePrefix /path/to/NEW/output/directory/JMDc3 --readFilesCommand gunzip -c --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile /proj/seq/data/MM10_UCSC/Annotation/Genes/genes.gtf --outWigType wiggle 

