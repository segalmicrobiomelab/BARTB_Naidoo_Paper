### Scripts used for transcriptome analysis

#############################################################################################
#############################################################################################
#############################################################################################

# Set working directory
setwd("C:/Users/Charissa/Dropbox/bartb/")

# Load Packages
library(DESeq2)
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(pathfindR)
library(scales)
library(data.table)
library(fBasics)
library(forcats)
library(vegan)
library(dplyr)
library(MetaboSignal)

# Load counts
count_data = read.table(file = "rna.rawcounts.csv", header = T, sep = ",", row.names=1)
head(count_data)

# Load metadata
sampleinfo <- read.table("rna.map.txt", header=T, sep="\t", row.names=1)

# Sampleinfo <- as.matrix((sampleinfo))
colnames(sampleinfo)

# Round off counts
counts <- round(count_data)
head(counts)

#########################################################
#########################################################
#########################################################

### DIFFERENTIALLY ABUNDANT GENES

# Make CountData and MetaData into DESEq Object; Choose the comparison Variable as design
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sampleinfo,
                              design= ~ tb_status_biome)

# Estimate Factors of DESeq Object
dds <- estimateSizeFactors(dds)

# Pruning before DESeq:
#This would filter out genes where there are less than 3 samples with normalized counts greater than or equal to 100.
idx100 <- rowSums( counts(dds, normalized=TRUE) >= 100 ) >= 3
dds <- dds[idx100,]
dds

# Differential Expression Analysis
# Drop Unencessary Levels
dds$tb_status_biome <- droplevels(dds$tb_status_biome)

# Set the baseline level --> positive is upregulated in cases; negative is downregulated
dds$tb_status_biome <- relevel(dds$tb_status_biome, ref ="No_TB")

# Run DESeq
dds<- DESeq(dds, test="Wald", fitType="local")

res <- results(dds, cooksCutoff=FALSE)

res = res[order(res$padj, na.last = NA), ]

# Annotation
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db) #To get a list of all available key types

convertIDs <- function( ids, fromKey, toKey, db, ifMultiple=c( "putNA", "useFirst" ) ) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select( 
    db, keys=ids, keytype=fromKey, columns=c(fromKey,toKey) ) )
  if( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]   
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ] }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}

# ENTREZID:
res$entrezgene = convertIDs(row.names(res), "ENSEMBL", "ENTREZID", org.Hs.eg.db)

# Gene symbols:
res$hgnc_symbol = convertIDs(row.names(res), "ENSEMBL", "SYMBOL", org.Hs.eg.db)

# Gene name:
res$gene_name = convertIDs(row.names(res), "ENSEMBL", "GENENAME", org.Hs.eg.db)

# Remove NAs in all newly annotated columns:
res = res[order(res$hgnc_symbol, na.last = NA), ]
res = res[order(res$entrezgene, na.last = NA), ]
res = res[order(res$gene_name, na.last = NA), ]

# Save image
save.image(file="C:/Users/Charissa/Documents/16s/seq_batch_2/bar_tb/final/rna/rna.RData")

# Convert resuts table into a data.frame
res <- as.data.frame(res)

# Select genes to save
gene.to.save <-as.character(rownames(res))
gene.to.save_entrezgene <-as.character(res$entrezgene)
gene.to.save_hgnc_symbol <-as.character(res$hgnc_symbol)

# From main table we should get the mean across the row of the table
ko.table.df <- data.frame(assay(dds))
ko.table.df.meanRA <- rowMeans(ko.table.df)

# Subset and reorder just the IDs that we have
ko.table.df.meanRA.save <- ko.table.df.meanRA[gene.to.save]

# Add the abundance data for the res dataframe
res$abundance <- ko.table.df.meanRA.save

# Subset table for groups being analysed to get their separate abundances
# Extract and subset count data
countdata = assay(dds)
coldata = colData(dds)
cases.pruned.table = countdata[, coldata$tb_status_biome %in% c("TB")]
sc.pruned.table = countdata[, coldata$tb_status_biome %in% c("No_TB")]

# From pruned table we should get the mean across the row of the table
cases.pruned.table.df <- data.frame(cases.pruned.table)
cases.pruned.table.df.meanRA <- rowMeans(cases.pruned.table.df)

sc.pruned.table.df <- data.frame(sc.pruned.table)
sc.pruned.table.df.meanRA <- rowMeans(sc.pruned.table.df)

# Subset AND reorder just the genes that we have 
cases.pruned.table.df.meanRA.save <- cases.pruned.table.df.meanRA[gene.to.save]
sc.pruned.table.df.meanRA.save <- sc.pruned.table.df.meanRA[gene.to.save]

# Add the abundance data for the res dataframe
res$abundance.cases <- cases.pruned.table.df.meanRA.save
res$abundance.sc <- sc.pruned.table.df.meanRA.save

# Set Names of Results Table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","entrez_ID", "hgnc_symbol", "gene_name", "abundance", "abundance.cases", "abundance.sc"))
write.table(res,file="rna.abundance.txt", sep="\t", col.names = NA, row.names = TRUE)

# Keep only the variables you need for pathway analysis
res1 <- res[,c("entrez_ID","logFC","pvalue","adj.P.Val")]

# Write Tables to TXT file
write.table(res1,file="rna.IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Set Theme for Figures
theme <-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")

# Set alpha
alpha = 0.01

# Compute FDR in a log scales
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plat
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "red"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "darkgreen"

# Create a Variable for the size of the dots in the Volcano Plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Select genes with a defined p-value (DESeq2 assigns NA to some genes)
genes.to.plot <- !is.na(res$pvalue)

# Check the range of the LogFC
range(res[genes.to.plot, "logFC"])

# Volcano plot
pdf(file="rna.volcano.fdr.0.01.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=hgnc_symbol)) +
  geom_point(color=cols, alpha=0.7) + #Chose Colors and size for dots
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$logFC>2, as.character(res$hgnc_symbol),'')),size=2,force=25,segment.colour="darkgrey",segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$logFC<=-2, as.character(res$hgnc_symbol),'')),size=2,force=25,segment.colour="darkgrey",segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  theme(legend.position = "none") +
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") +
  xlab("Effect size: log2(fold-change)") +
  ylab("-log10(adjusted p-value)") +
  #ylim(0,20)+
  theme
dev.off()

#########################################################
#########################################################
#########################################################

### BRAY-CURTIS ANALYSIS

# Create Distance Matrix
vegdist   = vegdist(t(assay(dds)), method = "bray")

# Calculate Adonis for p-value for TB status
adonis(vegdist ~ tb_status_biome, data=data.frame(colData(dds)))

# Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

# Calculate sample variance for each PC
vars <- apply(CmdScale, 2, var)

# Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

# Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(dds), by = "row.names", all.x = TRUE)

# Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the Centroid Value - for TB
centroids <- aggregate(cbind(PC1,PC2)~ tb_status_biome,data= newResults, mean)

# Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="tb_status_biome",suffixes=c("",".centroid"))

# Plot
pdf("rna.bray.curtis.pdf", height = 10, width = 10)
ggplot(newResults, aes(PC1, PC2, color= tb_status_biome)) + 
  geom_point(size=5) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values=c("darkgreen", "red", "grey", "blue", "purple")) + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= tb_status_biome), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= tb_status_biome))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("SCs", "Cases")), size=10) + 
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour= "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()