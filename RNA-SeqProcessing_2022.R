##Load DESeq2
source("https://bioconductor.org/biocLite.R")
library(DESeq2)
library(ggplot2)

#Set your working directory on your desktop
setwd("~/Desktop/RNAseq002")

DESeq2::DESeqDataSetFromMatrix()

##Read in your data of gene counts for all replicates (output from STAR) and sample info table
Counts_bbmap <- read.table("RNAseq002_ReadsPerGene_Master.txt", header = TRUE, row.names=1)
sampleInfo <- read.csv("Sample_Info.csv")
sampleInfo <- data.frame(sampleInfo)
dds <- DESeqDataSetFromMatrix(countData = Counts_bbmap, colData = sampleInfo, design = ~Genotype)

dds <- DESeq(dds)
resultsNames(dds)
vsd <- vst(dds, blind=FALSE)

#Normalized counts, standard from DESeq2
test <- counts(dds, normalized=T)
write.table(test, "WT_Normalized_counts.txt")

##PCA Plot from DESeq2
print(colData(vsd))
plotPCA(vsd, intgroup=c("Sample", "Genotype"))

##Customized PCA using ggplot
pcaData <- plotPCA(vsd, intgroup=c("Sample", "Genotype"), returnData=TRUE)
ggplot(pcaData, aes(PC1, PC2, color=Genotype))+
  geom_point(size=4) + scale_color_manual(values=c("#ba6000", "#e07400", "#ff8503", "#fca447", "#550091", "#7102bf", "#9000f5", "#c076f5", "#018c2d", "#00ba3b", "#00fc50", "#84f5a6")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

##MA Plot from DESeq2 - change genotypes you want to compare, listed first = top of MA
res <- results(dds, contrast=c("Genotype", "P5A-siSA2", "WT-siGLO"))
resLFC <- lfcShrink(dds, contrast=c("Genotype", "P5A-siSA2", "WT-siGLO"))
summary(resLFC)
plotMA(resLFC, ylim=c(-3, 3), colSig="blue")

##Print differential gene values to a file based on padj < 0.1
resSign <- subset(resLFC, resLFC$padj < 0.1)                  
write.csv(as.data.frame(resSign), file="P5A-siSA2_WT-siGLO-DE.txt")

#Print values of all genes regardless of signifiance
resSign <- resLFC
write.csv(as.data.frame(resSign), file="P5A-siSA2_WT-siGLO-ALL-GENES.txt")         

##Heatmap of combined log2 fold changes, tables made using dplyr - create matrix first, then use pheatmap or heatmap
LFC <- read.table("P5A_and_P5B_All_DE_LFC.txt", header=TRUE, sep="")

row.names(LFC) <- LFC$Gene
LFC=subset(LFC,select =-Gene)
LFC_matrix <- data.matrix(LFC)

library("RColorBrewer")
heatmap(LFC_matrix, Colv = NA, na.rm=TRUE, col = brewer.pal(n=8, name = 'Purples'))

library("pheatmap")

#make a 'smooth' heatmap #two colors
#n is always sum(lengths)-1
my_palette<- colorRampPalette(c("gray", "mediumslateblue"))(n=1000)
col_breaks = c(seq(0,2000, length = 500), seq(2001,4000, length = 501))
pheatmap(LFC_matrix, cellwidth = 50, breaks= col_breaks, cluster_rows = TRUE,  cluster_cols = FALSE,show_rownames = FALSE, color = my_palette )


#make a 'smooth' heatmap #three colors examples
#n is always sum(lengths)-1
my_palette<- colorRampPalette(c("dodgerblue","black", "yellow"))(n=49)
col_breaks = c(seq(-3,-0.200000000001, length=10), seq(-0.20,0.20, length=30), seq(0.200000001,3, length=10))
pheatmap(LFC_matrix, cellwidth = 50, breaks = col_breaks, cluster_rows = FALSE,  cluster_cols = FALSE,show_rownames = FALSE, color=my_palette )



##Get Help
browseVignettes("DESeq2")
