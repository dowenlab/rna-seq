######### descriptions of files in this directory: ##################


######### RNA-seq analysis ################

RNA-SeqProcessing_2022.R #R script for RNA-seq processing. Counts_bbmap input should resemble RNAseq_ReadsPerGene_Master.txt in this folder where gene counts for every sample are listed (generated from STAR_RNA-seq.sh). Sample_Info input should also resemble Sample_info.csv in this folder. Contains basic DESeq running and plot making

RNAseq_ReadsPerGene_Master.txt #example input file for RNAseq data to be analyzed with RNA-SeqProcessing_2022.R script. gene counts from STAR alignment for all replicates of RNAseq listed in tab delimited format. make sure when opening in excel to File>Open>Columns as General (not text) to prevent genes turning into dates

Sample_info.csv #example sample input file for RNAseq data to be analyzed with RNA-SeqProcessing_2022.R script. this tells DESeq what to group together based on replicates

STAR_RNA-seq.sh #Bash script to align RNAseq fastq files and get gene count files to process in R. creates a directory for each sample where genecount file can be found. 


######### Reference files ################

AllGenesBioMart_Strand_NoNA_FIX.bed #Gene coordinates with strand info - origin BioMart/Ensembl

Genes-in-PDs_2014_Dowen.bed #coordinates of polycomb domains from Dowen et al 2014

Genes-in-SDs_2014_Dowen.bed #coordinates of super-enhancer domains from Dowen et al 2014

genes.gtf #gtf (Gene transfer format) file for mm10 genes

UCSC-TSS-mm10_PROMOTER-COORDS.bed #coordinates for gene promoters in mm10 from UCSC. Use for getting coordinates in python script from a list of genes



########## Downstream analysis #############

MakeBed-DESeqLFC-Coords.py #python script to get coordinates for a list of genes. change your file name (file1) within script before using. uses UCSC-TSS-mm10_PROMOTER-COORDS.txt as reference coordinates. to use: module load python| python3 MakeBed-DESeqLFC-Coords.py > output.bed

RNAseqBigWig.sh #sbatch script for making a bigwig file from RNA-seq bam file. this realigns using bowtie, but could also use STAR aligned bam
