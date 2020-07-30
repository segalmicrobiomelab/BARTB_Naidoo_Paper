### Scripts used for PICRUSt analysis

#############################################################################################
#############################################################################################
#############################################################################################

### In the UCT server

# Activate the  QIIME environment
source activate_qiime.sh

# Normalize
/opt/exp_soft/python-2.7.3/bin/normalize_by_copy_number.py -i /home/cnaidoo/bar_tb/picrust/otu_table.stool.biom -o /home/cnaidoo/bar_tb/picrust/otu_table.stool.norm.biom

# Predict
/opt/exp_soft/python-2.7.3/bin/predict_metagenomes.py -i /home/cnaidoo/bar_tb/picrust/otu_table.stool.norm.biom -c /home/cnaidoo/ko_13_5_precalculated.tab -o /home/cnaidoo/bar_tb/picrust/otu_table.stool.norm.predicted.metagenome.biom 

# Categorize
/opt/exp_soft/python-2.7.3/bin/categorize_by_function.py -i /home/cnaidoo/bar_tb/picrust/otu_table.stool.norm.predicted.metagenome.biom -c KEGG_Pathways -l 3 -o /home/cnaidoo/bar_tb/picrust/otu_table.stool.norm.predicted.metagenome.level_3.biom
/opt/exp_soft/python-2.7.3/bin/categorize_by_function.py -i /home/cnaidoo/bar_tb/picrust/otu_table.stool.norm.predicted.metagenome.biom -c KEGG_Pathways -l 2 -o /home/cnaidoo/bar_tb/picrust/stool.norm.predicted.metagenome.level_2.biom
/opt/exp_soft/python-2.7.3/bin/categorize_by_function.py -i /home/cnaidoo/bar_tb/picrust/otu_table.stool.norm.predicted.metagenome.biom -c KEGG_Pathways -l 1 -o /home/cnaidoo/bar_tb/picrust/stool.norm.predicted.metagenome.level_1.biom

# Convert BIOM files to TXT 
biom convert -i /home/cnaidoo/bar_tb/picrust/otu_table.stool.norm.predicted.metagenome.level_1.biom -o /home/cnaidoo/bar_tb/picrust/otu_table.stool.norm.predicted.metagenome.level_1.txt --to-tsv --header-key KEGG_Pathways 
biom convert -i /home/cnaidoo/bar_tb/picrust/otu_table.stool.norm.predicted.metagenome.level_2.biom -o /home/cnaidoo/bar_tb/picrust/otu_table.stool.norm.predicted.metagenome.level_2.txt --to-tsv --header-key KEGG_Pathways 
biom convert -i /home/cnaidoo/bar_tb/picrust/otu_table.stool.norm.predicted.metagenome.level_3.biom -o /home/cnaidoo/bar_tb/picrust/otu_table.stool.norm.predicted.metagenome.level_3.txt --to-tsv --header-key KEGG_Pathways 

### In Excel, remove all pathways except those involved in metabolism and subset in 3 files (cases and SCs, cases and CCCs, SCs and CCSCs)

### In R

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
library(phyloseq)

### CASES VS. SYMPTOMATIC CONTROLS (SCs)

# Load file
count_data = read.table(file = "otu_table.stool.norm.predicted.metagenome.level_3.patients.csv", header = T, sep = ",", row.names=1)
head(count_data)

# Load metadata
sampleinfo <- read.table("stool.map.patients.txt", header=T, sep="\t", row.names=1)
colnames(sampleinfo)

# Round off counts
counts <- round(count_data)
head(counts)

# Make CountData and MetaData into DESEq Object; Choose the comparison Variable as design
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sampleinfo,
                              design= ~ group)

# Estimate Factors of DESeq Object
dds <- estimateSizeFactors(dds)

# Drop Unencessary Levels
dds$group <- droplevels(dds$group)

# Set the baseline level
dds$group <- relevel(dds$group, ref ="SCs")

# Run DESeq
dds<- DESeq(dds, test="Wald", fitType="local")

res <- results(dds, cooksCutoff=FALSE)

res = res[order(res$padj, na.last = NA), ]

# Convert resuts table into a data.frame
res <- as.data.frame(res)

gene.to.save <-as.character(rownames(res))

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
cases.pruned.table = countdata[, coldata$group %in% c("Cases")]
sc.pruned.table = countdata[, coldata$group %in% c("SCs")]

# From relative table we should get the mean across the row of the OTU table
cases.pruned.table.df <- data.frame(cases.pruned.table)
cases.pruned.table.df.meanRA <- rowMeans(cases.pruned.table.df)

sc.pruned.table.df <- data.frame(sc.pruned.table)
sc.pruned.table.df.meanRA <- rowMeans(sc.pruned.table.df)

# Subset AND reorder just the otus that we have 
cases.pruned.table.df.meanRA.save <- cases.pruned.table.df.meanRA[gene.to.save]
sc.pruned.table.df.meanRA.save <- sc.pruned.table.df.meanRA[gene.to.save]

# Add the abundance data for the res dataframe
res$abundance.cases <- cases.pruned.table.df.meanRA.save
res$abundance.sc <- sc.pruned.table.df.meanRA.save

# Set Names of Results Table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","abundance", "abundance.cases", "abundance.sc"))
write.table(res,file="pathways.abundance.patients.table.txt", sep="\t", col.names = NA, row.names = TRUE)

# Set Theme for Figures
theme <-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
alpha = 0.20

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

#Check the range of the LogFC
range(res[genes.to.plot, "logFC"])

# Volcano plot
pdf(file="stool.picrust.patients.volcano.FDR.0.2.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
  geom_point(color=cols, size = ifelse(res$logFC>=0.1 & res$adj.P.Val < alpha, 0.0001 * res$abundance.cases, ifelse(res$logFC<=-0.1 & res$adj.P.Val < alpha, 0.0001 * res$abundance.sc,2)), alpha=0.7) + #Chose Colors and size for dots
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$logFC>2, as.character(res$Gene.symbol),'')),size=2,force=25,segment.colour="darkgrey",segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$logFC<=-2, as.character(res$Gene.symbol),'')),size=2,force=25,segment.colour="darkgrey",segment.alpha=0.5) + 
  theme(legend.position = "none") +
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)") + #label Y Axis
  theme #Set Theme
dev.off() 

### CASES VS. CLOSE CONTACTS (CCCs)

# Load file
count_data = read.table(file = "otu_table.stool.norm.predicted.metagenome.level_3.cases.vs.ccc.csv", header = T, sep = ",", row.names=1)
head(count_data)

# Load metadata
sampleinfo <- read.table("stool.map.cases.vs.ccc.txt", header=T, sep="\t", row.names=1)
colnames(sampleinfo)

# Round off counts
counts <- round(count_data)
head(counts)

# Make CountData and MetaData into DESEq Object; Choose the comparison Variable as design
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sampleinfo,
                              design= ~ group)

# Estimate Factors of DESeq Object
dds <- estimateSizeFactors(dds)

# Drop Unencessary Levels
dds$group <- droplevels(dds$group)

# Set the baseline level
dds$group <- relevel(dds$group, ref ="CCCs")

# Run DESeq
dds<- DESeq(dds, test="Wald", fitType="local")

res <- results(dds, cooksCutoff=FALSE)

res = res[order(res$padj, na.last = NA), ]

# Convert resuts table into a data.frame
res <- as.data.frame(res)

gene.to.save <-as.character(rownames(res))

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
cases.pruned.table = countdata[, coldata$group %in% c("Cases")]
ccc.pruned.table = countdata[, coldata$group %in% c("CCCs")]

# From relative table we should get the mean across the row of the OTU table
cases.pruned.table.df <- data.frame(cases.pruned.table)
cases.pruned.table.df.meanRA <- rowMeans(cases.pruned.table.df)

ccc.pruned.table.df <- data.frame(ccc.pruned.table)
ccc.pruned.table.df.meanRA <- rowMeans(ccc.pruned.table.df)

# Subset AND reorder just the otus that we have 
cases.pruned.table.df.meanRA.save <- cases.pruned.table.df.meanRA[gene.to.save]
ccc.pruned.table.df.meanRA.save <- ccc.pruned.table.df.meanRA[gene.to.save]

# Add the abundance data for the res dataframe
res$abundance.cases <- cases.pruned.table.df.meanRA.save
res$abundance.ccc <- ccc.pruned.table.df.meanRA.save

# Set Names of Results Table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","abundance", "abundance.cases", "abundance.ccc"))
write.table(res,file="pathways.abundance.cases.vs.CCCs.table.txt", sep="\t", col.names = NA, row.names = TRUE)

# Set Theme for Figures
theme <-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
alpha = 0.20

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

#Check the range of the LogFC
range(res[genes.to.plot, "logFC"])

# Volcano plot
pdf(file="stool.picrust.cases.vs.CCCs.volcano.FDR.0.2.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
  geom_point(color=cols, size = ifelse(res$logFC>=0.1 & res$adj.P.Val < alpha, 0.0001 * res$abundance.cases, ifelse(res$logFC<=-0.1 & res$adj.P.Val < alpha, 0.0001 * res$abundance.ccc,2)), alpha=0.7) + #Chose Colors and size for dots
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$logFC>2, as.character(res$Gene.symbol),'')),size=2,force=25,segment.colour="darkgrey",segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$logFC<=-2, as.character(res$Gene.symbol),'')),size=2,force=25,segment.colour="darkgrey",segment.alpha=0.5) + 
  theme(legend.position = "none") +
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)") + #label Y Axis
  theme #Set Theme
dev.off() 


### SYMPTOMATIC CONTROLS (SCs) VS. CLOSE CONTACTS (CCSCs)

# Load file
count_data = read.table(file = "otu_table.stool.norm.predicted.metagenome.level_3.sc.vs.ccsc.csv", header = T, sep = ",", row.names=1)
head(count_data)

# Load metadata
sampleinfo <- read.table("stool.map.sc.vs.ccsc.txt", header=T, sep="\t", row.names=1)
colnames(sampleinfo)

# Round off counts
counts <- round(count_data)
head(counts)

# Make CountData and MetaData into DESEq Object; Choose the comparison Variable as design
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sampleinfo,
                              design= ~ group)

# Estimate Factors of DESeq Object
dds <- estimateSizeFactors(dds)

# Drop Unencessary Levels
dds$group <- droplevels(dds$group)

# Set the baseline level
dds$group <- relevel(dds$group, ref ="CCSCs")

# Run DESeq
dds<- DESeq(dds, test="Wald", fitType="local")

res <- results(dds, cooksCutoff=FALSE)

res = res[order(res$padj, na.last = NA), ]

# Convert resuts table into a data.frame
res <- as.data.frame(res)

gene.to.save <-as.character(rownames(res))

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
sc.pruned.table = countdata[, coldata$group %in% c("SCs")]
ccsc.pruned.table = countdata[, coldata$group %in% c("CCSCs")]

# From relative table we should get the mean across the row of the OTU table
sc.pruned.table.df <- data.frame(sc.pruned.table)
sc.pruned.table.df.meanRA <- rowMeans(sc.pruned.table.df)

ccsc.pruned.table.df <- data.frame(ccsc.pruned.table)
ccsc.pruned.table.df.meanRA <- rowMeans(ccsc.pruned.table.df)

# Subset AND reorder just the otus that we have 
sc.pruned.table.df.meanRA.save <- sc.pruned.table.df.meanRA[gene.to.save]
ccsc.pruned.table.df.meanRA.save <- ccsc.pruned.table.df.meanRA[gene.to.save]

# Add the abundance data for the res dataframe
res$abundance.sc <- sc.pruned.table.df.meanRA.save
res$abundance.ccsc <- ccsc.pruned.table.df.meanRA.save

# Set Names of Results Table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","abundance", "abundance.sc", "abundance.ccsc"))
write.table(res,file="pathways.abundance.SCs.vs.CCSCs.table.txt", sep="\t", col.names = NA, row.names = TRUE)

# Set Theme for Figures
theme <-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
alpha = 0.20

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

#Check the range of the LogFC
range(res[genes.to.plot, "logFC"])

# Volcano plot
pdf(file="stool.picrust.SCs.vs.CCSCs.volcano.FDR.0.2.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
  geom_point(color=cols, size = ifelse(res$logFC>=0.1 & res$adj.P.Val < alpha, 0.0001 * res$abundance.sc, ifelse(res$logFC<=-0.1 & res$adj.P.Val < alpha, 0.0001 * res$abundance.ccsc,2)), alpha=0.7) + #Chose Colors and size for dots
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$logFC>2, as.character(res$Gene.symbol),'')),size=2,force=25,segment.colour="darkgrey",segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$logFC<=-2, as.character(res$Gene.symbol),'')),size=2,force=25,segment.colour="darkgrey",segment.alpha=0.5) + 
  theme(legend.position = "none") +
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)") + #label Y Axis
  theme #Set Theme
dev.off() 