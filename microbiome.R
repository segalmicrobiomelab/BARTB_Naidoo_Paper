### Scripts used for microbiome analysis

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
library(phyloseq)

# Load files needed
file = "otu_table.biom"
map = "master.map.txt"

# Load the abundance table and mapping table 
abundance.table = import_biom(file, taxaPrefix=F) 
mapping.table=sample_data(read.table(map, header=T, sep="\t", row.names=1))   

# Create phyloseq object
physeq=phyloseq(otu_table(abundance.table),tax_table(abundance.table), mapping.table)

# Give a colnames to separate different taxonomic levels
colnames(tax_table(physeq))=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "OTU")

# Load the tree file 
treefile = "97_otus_unannotated.tree"
tree.obj = import_qiime(treefilename = treefile)

# Now merge the three separate phyloseq objects into a single object
otu.table = merge_phyloseq(physeq, mapping.table, tree.obj)
rownames(sample_data(otu.table))
colnames(sample_data(otu.table))

# Remove taxa with 0 abundance
otu.table = subset_taxa(otu.table, rowSums(otu_table(otu.table)) != 0)

# To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}

# Create OTU tables
otu.relative.table = transformSampleCounts(otu.table, normalizeSample)
Phylum.rel.table = tax_glom(otu.relative.table, taxrank = "Phylum")
Class.rel.table = tax_glom(otu.relative.table, taxrank = "Class")
Order.rel.table = tax_glom(otu.relative.table, taxrank = "Order")
Family.rel.table = tax_glom(otu.relative.table, taxrank = "Family")
Genus.rel.table = tax_glom(otu.relative.table, taxrank = "Genus")
OTU.rel.table = tax_glom(otu.relative.table, taxrank = "OTU")

# Save 
save.image(file="all.samples.RData")

### ALPHA DIVERSITY

alpha <- estimate_richness(physeq)
write.csv(alpha, "alpha.diversity.csv")

# Plots done in GraphPad Prism 

#############################################################################################
#############################################################################################
#############################################################################################

### ORAL WASHES

# Subset oral washes
oral.otu.relative.table = subset_samples(otu.relative.table, sample_c %in% c(1))

## BETA DIVERSITY

# Create Distance Matrix
vegdist = distance(oral.otu.relative.table, "wunifrac")

# Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

# Calculate Sample variance for each PC
vars <- apply(CmdScale, 2, var)

# Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

# Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(oral.otu.relative.table), by = "row.names", all.x = TRUE)

# Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ group,data= newResults, mean)

# Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="group",suffixes=c("",".centroid"))

pdf("oral.wash.wuniF.pdf", height = 10, width = 15)
ggplot(newResults, aes(PC1, PC2, color= group)) +
  geom_point(size=3,alpha=0.7) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values=c("red", "purple", "cyan", "forestgreen", "grey")) + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= group), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= group, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=group, size=10)) +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Subset cases and SCs
patients.oral.otu.relative.table = subset_samples(oral.otu.relative.table, group_c %in% c(0,1))
# Create distance object
patients = distance(patients.oral.otu.relative.table, "wunifrac")
# Statistics
adonis(patients ~ sample_data(patients.oral.otu.relative.table)$group)

# Subset cases and CCCs
cases.ccc.oral.otu.relative.table = subset_samples(oral.otu.relative.table, group_c %in% c(1,2))
# Create distance object
cases = distance(cases.ccc.oral.otu.relative.table, "wunifrac")
# Statistics
adonis(cases ~ sample_data(cases.ccc.oral.otu.relative.table)$group)

# Subset SCs and CCSCs
sc.ccs.oral.otu.relative.table = subset_samples(oral.otu.relative.table, group_c %in% c(0,3))
# Create distance object
sc = distance(sc.ccs.oral.otu.relative.table, "wunifrac")
# Statistics
adonis(sc ~ sample_data(sc.ccs.oral.otu.relative.table)$group)

#########################################################
#########################################################
#########################################################

### BETA DIVERSITY IN EACH GROUP BY HIV STATUS

## Cases
# Subset cases
cases.oral.otu.relative.table = subset_samples(oral.otu.relative.table, group_c %in% c(1))

# Create Distance Matrix
vegdist = distance(cases.oral.otu.relative.table, "wunifrac")

# Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

# Calculate Sample variance for each PC
vars <- apply(CmdScale, 2, var)

# Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

# Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(cases.oral.otu.relative.table), by = "row.names", all.x = TRUE)

# Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ hiv,data= newResults, mean)

# Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="hiv",suffixes=c("",".centroid"))

pdf("oral.wash.wunif.hiv.cases.pdf", height = 10, width = 10)
ggplot(newResults, aes(PC1, PC2, color= hiv)) + # Graph PC1 and PC2
  geom_point(size=3,alpha=0.7) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  scale_color_manual(values=c("orange","red", "grey")) + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= hiv), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= hiv, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=hiv, size=10)) +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Statistics
adonis(vegdist ~ sample_data(cases.oral.otu.relative.table)$hiv)


## Symptomatic controls (SCs)
# Subset SCs
sc.oral.otu.relative.table = subset_samples(oral.otu.relative.table, group_c %in% c(0))

# Create Distance Matrix
vegdist = distance(sc.oral.otu.relative.table, "wunifrac")

# Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

# Calculate Sample variance for each PC
vars <- apply(CmdScale, 2, var)

# Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

# Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(sc.oral.otu.relative.table), by = "row.names", all.x = TRUE)

# Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ hiv,data= newResults, mean)

# Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="hiv",suffixes=c("",".centroid"))

# Plot
pdf("oral.wash.wunif.hiv.SCs.pdf", height = 10, width = 10)
ggplot(newResults, aes(PC1, PC2, color= hiv)) + # Graph PC1 and PC2
  geom_point(size=3,alpha=0.7) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  scale_color_manual(values=c("orange","red", "grey")) + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= hiv), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= hiv, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=hiv, size=10)) +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Statistics
adonis(vegdist ~ sample_data(sc.oral.otu.relative.table)$hiv)

## Close contacts of cases (CCCs)
# Subset CCCs
ccc.oral.otu.relative.table = subset_samples(oral.otu.relative.table, group_c %in% c(2))

#subset CCCs with only known hiv status
ccc.oral.otu.relative.table = subset_samples(ccc.oral.otu.relative.table, known_hiv_status %in% c(1))

# Create Distance Matrix
vegdist = distance(ccc.oral.otu.relative.table, "wunifrac")

# Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

# Calculate Sample variance for each PC
vars <- apply(CmdScale, 2, var)

# Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

# Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(ccc.oral.otu.relative.table), by = "row.names", all.x = TRUE)

# Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ hiv,data= newResults, mean)

# Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="hiv",suffixes=c("",".centroid"))

# Plot
pdf("oral.wash.wunif.hiv.CCCs.pdf", height = 10, width = 10)
ggplot(newResults, aes(PC1, PC2, color= hiv)) + # Graph PC1 and PC2
  geom_point(size=3,alpha=0.7) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  scale_color_manual(values=c("orange","red", "grey")) + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= hiv), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= hiv, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=hiv, size=10)) +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Statistics
adonis(vegdist ~ sample_data(ccc.oral.otu.relative.table)$hiv)

## Close contacts of symptomatic controls (CCSCs)
# Subset CCSCs
ccsc.oral.otu.relative.table = subset_samples(oral.otu.relative.table, group_c %in% c(3))

#subset CCSCs with only known hiv status
ccsc.oral.otu.relative.table = subset_samples(ccsc.oral.otu.relative.table, known_hiv_status %in% c(1))

# Create Distance Matrix
vegdist = distance(ccsc.oral.otu.relative.table, "wunifrac")

# Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

# Calculate Sample variance for each PC
vars <- apply(CmdScale, 2, var)

# Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

# Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(ccsc.oral.otu.relative.table), by = "row.names", all.x = TRUE)

# Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ hiv,data= newResults, mean)

# Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="hiv",suffixes=c("",".centroid"))

# Plot
pdf("oral.wash.wunif.hiv.CCSCs.pdf", height = 10, width = 10)
ggplot(newResults, aes(PC1, PC2, color= hiv)) + # Graph PC1 and PC2
  geom_point(size=3,alpha=0.7) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  scale_color_manual(values=c("orange","red", "grey")) + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= hiv), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= hiv, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=hiv, size=10)) +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

#Statistics not done because HIV-positive group only had one patient

#########################################################
#########################################################
#########################################################

### BETA DIVERSITY IN CLOSE CONTACTS OF CASES (CCCs) VS. CLOSE CONTACTS OF SYMPTOMATIC CONTROLS (CCSCs)

# Subset CCCs and CCSCs
contacts.oral.otu.relative.table = subset_samples(oral.otu.relative.table, group_c %in% c(2,3))
rownames(sample_data(contacts.oral.otu.relative.table))

# Create Distance Matrix
vegdist = distance(contacts.oral.otu.relative.table, "wunifrac")

# Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

# Calculate Sample variance for each PC
vars <- apply(CmdScale, 2, var)

# Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#  PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(contacts.oral.otu.relative.table), by = "row.names", all.x = TRUE)

# Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ group,data= newResults, mean)

# Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="group",suffixes=c("",".centroid"))

# Plot
pdf("oral.wash.wuniF.contacts.pdf", height = 10, width = 15)
ggplot(newResults, aes(PC1, PC2, color= group)) + # Graph PC1 and PC2
  geom_point(size=3,alpha=0.7) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  scale_color_manual(values=c("purple", "cyan")) + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= group), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= group, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=group, size=10)) + 
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Statistics
adonis(vegdist ~ sample_data(contacts.oral.otu.relative.table)$group)

#########################################################
#########################################################
#########################################################

### DIFFERENTIAL ABUNDANCE ANALYSIS WITH DESEQ2

# Set Theme For Figures
theme <-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),
              axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), 
              legend.position="none")

# Load file
#load(file="all.samples.RData")

# Set alpha
alpha <- 0.20


## CASES VS SYMPTOMATIC CONTROLS (SCs)

# Subset OTU table to oral washes
subset.otu.table = subset_samples(otu.table, sample_c %in% c(1))

# Subset OTU table to cases and SCs
subset.otu.table = subset_samples(subset.otu.table, group_c %in% c(0,1))

# Subset OTU table to genus level
subset.genus.table = tax_glom(subset.otu.table, taxrank = "Genus")

# Prune data to less than 100 genera remaining
filtered.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.15 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(filtered.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

# Create relative abundance table
# To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# genus table: pruned.genus.table
# genus relative table: pruned.genus.rel.table
# variable: group

# Convert Phyloseq Object to DESEq object, correcting for any potential confounders
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ group)

# Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$group <- droplevels(diagdds$group)

# Choose the 'control' in your comparison variable --> positive is upregulated in cases, negative is down-regulated
diagdds$group <- relevel(diagdds$group, ref ="SCs")

# Run the differential Analysis
diagdds<- DESeq(diagdds)

# Output the results from DESeq into a table
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant genes
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

# Tables

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$OTU,rownames(res))
# Replace spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

# Convert resuts table into a data.frame
res <- as.data.frame(res)

#Set names of results table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))

# Make the full trail the first column
res$names <- res$Gene.symbol

res$Gene.symbol <- res$row2

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
cases.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("1"))
sc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("0"))

# From relative table we should get the mean across the row of the otu table
cases.pruned.genus.rel.table.df <- data.frame(otu_table(cases.pruned.genus.rel.table))
cases.pruned.genus.rel.table.df.meanRA <- rowMeans(cases.pruned.genus.rel.table.df)

sc.pruned.genus.rel.table.df <- data.frame(otu_table(sc.pruned.genus.rel.table))
sc.pruned.genus.rel.table.df.meanRA <- rowMeans(sc.pruned.genus.rel.table.df)

# Subset AND reorder just the OTUs that we have 
cases.pruned.genus.rel.table.df.meanRA.save <- cases.pruned.genus.rel.table.df.meanRA[otu.to.save]
sc.pruned.genus.rel.table.df.meanRA.save <- sc.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res dataframe
res$abundance.cases <- cases.pruned.genus.rel.table.df.meanRA.save
res$abundance.sc <- sc.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need for pathway analysis
res.IPA <- res[,c("Gene.symbol", "logFC", "pvalue", "adj.P.Val")]
res.1 <- res[,c("Gene.symbol", "abundance.cases", "abundance.sc", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Write tables to TXT file
write.table(res.IPA,file="taxa.oral.wash.patients.IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
write.table(res.1,file="taxa.oral.wash.patients.abundance.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Compute FDR in a log scales
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "red"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "darkgreen"

# Create a Variable for the size of the dots in the Volcano Plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="oral.wash.patients.volcano.fdr.0.2.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 2000 * res$abundance.cases, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 2000 * res$abundance.sc,2)),alpha=0.7) + #Chose colors and size for dots
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=3,force=25, segment.colour="grey",segment.alpha=0.2) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)") + #label Y Axis
  theme #Set Theme
dev.off()


## CASES VS THEIR CLOSE CONTACTS (CCCs)

# Subset OTU table to oral washes
subset.otu.table = subset_samples(otu.table, sample_c %in% c(1))

# Subset OTU table to cases and CCCs
subset.otu.table = subset_samples(subset.otu.table, group_c %in% c(1,2))

# Subset OTU table to genus level
subset.genus.table = tax_glom(subset.otu.table, taxrank = "Genus")

# Prune data to less than 100 genera remaining
filtered.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.15 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(filtered.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

# Create relative abundance table
# To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# genus table: pruned.genus.table
# genus relative table: pruned.genus.rel.table
# variable: group

# Convert Phyloseq Object to DESEq object, correcting for any potential confounders
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ pairs + group)

# Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$group <- droplevels(diagdds$group)

# Choose the 'control' in your comparison variable 
diagdds$group <- relevel(diagdds$group, ref ="CCCs")

# Run the differential Analysis
diagdds<- DESeq(diagdds)

# Output the results from DESeq into a table
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant genes
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

## Tables

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$OTU,rownames(res))
# Replace spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

# Convert resuts table into a data.frame
res <- as.data.frame(res)

#Set names of results table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))

# Make the full trail the first column
res$names <- res$Gene.symbol

res$Gene.symbol <- res$row2

# Get abundance data - use oral.otu.relative.table - THIS CODE ADDS REL ABUNDANCE AS THE DOT SIZE

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
cases.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("1"))
ccc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("2"))

# From relative table we should get the mean across the row of the otu table
cases.pruned.genus.rel.table.df <- data.frame(otu_table(cases.pruned.genus.rel.table))
cases.pruned.genus.rel.table.df.meanRA <- rowMeans(cases.pruned.genus.rel.table.df)

ccc.pruned.genus.rel.table.df <- data.frame(otu_table(ccc.pruned.genus.rel.table))
ccc.pruned.genus.rel.table.df.meanRA <- rowMeans(ccc.pruned.genus.rel.table.df)

# Subset AND reorder just the OTUs that we have 
cases.pruned.genus.rel.table.df.meanRA.save <- cases.pruned.genus.rel.table.df.meanRA[otu.to.save]
ccc.pruned.genus.rel.table.df.meanRA.save <- ccc.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res dataframe
res$abundance.cases <- cases.pruned.genus.rel.table.df.meanRA.save
res$abundance.ccc <- ccc.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need for pathway analysis
res.IPA <- res[,c("Gene.symbol", "logFC", "pvalue", "adj.P.Val")]
res.1 <- res[,c("Gene.symbol", "abundance.cases", "abundance.ccc", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Write tables to TXT file
write.table(res.IPA,file="taxa.oral.wash.cases.vs.CCCs.IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
write.table(res.1,file="taxa.oral.wash.cases.vs.CCCs.abundance.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Volcano plot

# Compute FDR in a log scales
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "red"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "darkgreen"

# Create a Variable for the size of the dots in the Volcano Plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="oral.wash.cases.vs.CCCs.volcano.fdr.0.2.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 2000 * res$abundance.cases, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 2000 * res$abundance.ccc,2)),alpha=0.7) + #Chose colors and size for dots
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=3,force=25, segment.colour="grey",segment.alpha=0.2) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)") + #label Y Axis
  theme #Set Theme
dev.off()



## SYMPTOMATIC CONTROLS (SCs) VS THEIR CLOSE CONTACTS (CCCs)

# Subset OTU table to oral washes
subset.otu.table = subset_samples(otu.table, sample_c %in% c(1))

# Subset OTU table to SCs and CCSCs
subset.otu.table = subset_samples(subset.otu.table, group_c %in% c(0,3))

# Subset OTU table to genus level
subset.genus.table = tax_glom(subset.otu.table, taxrank = "Genus")

# Prune data to less than 100 genera remaining
filtered.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.20 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(filtered.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

# Create relative abundance table
# To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# genus table: pruned.genus.table
# genus relative table: pruned.genus.rel.table
# variable: group

# Convert Phyloseq Object to DESEq object, correcting for any potential confounders
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ pairs + group)

# Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$group <- droplevels(diagdds$group)

# Choose the 'control' in your comparison variable 
diagdds$group <- relevel(diagdds$group, ref ="CCSCs")

# Run the differential Analysis
diagdds<- DESeq(diagdds)

# Output the results from DESeq into a table
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant genes
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

## TABLES

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$OTU,rownames(res))
# Replace spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

# Convert resuts table into a data.frame
res <- as.data.frame(res)

#Set names of results table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))

# Make the full trail the first column
res$names <- res$Gene.symbol

res$Gene.symbol <- res$row2

# Get abundance data - use oral.otu.relative.table - THIS CODE ADDS REL ABUNDANCE AS THE DOT SIZE

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
sc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("0"))
ccsc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("3"))

# From relative table we should get the mean across the row of the otu table
sc.pruned.genus.rel.table.df <- data.frame(otu_table(sc.pruned.genus.rel.table))
sc.pruned.genus.rel.table.df.meanRA <- rowMeans(sc.pruned.genus.rel.table.df)

ccsc.pruned.genus.rel.table.df <- data.frame(otu_table(ccsc.pruned.genus.rel.table))
ccsc.pruned.genus.rel.table.df.meanRA <- rowMeans(ccsc.pruned.genus.rel.table.df)

# Subset AND reorder just the OTUs that we have 
sc.pruned.genus.rel.table.df.meanRA.save <- sc.pruned.genus.rel.table.df.meanRA[otu.to.save]
ccsc.pruned.genus.rel.table.df.meanRA.save <- ccsc.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res dataframe
res$abundance.sc <- sc.pruned.genus.rel.table.df.meanRA.save
res$abundance.ccsc <- ccsc.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need for pathway analysis
res.IPA <- res[,c("Gene.symbol", "logFC", "pvalue", "adj.P.Val")]
res.1 <- res[,c("Gene.symbol", "abundance.sc", "abundance.ccsc", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Write tables to TXT file
write.table(res.IPA,file="taxa.oral.wash.SCs.vs.CCSCs.IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
write.table(res.1,file="taxa.oral.wash.SCs.vs.CCSCs.abundance.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

## Volcano plot

# Compute FDR in a log scales
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "red"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "darkgreen"

# Create a Variable for the size of the dots in the Volcano Plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="oral.wash.SCs.vs.CCSCs.volcano.fdr.0.2.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 500 * res$abundance.sc, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 500 * res$abundance.ccsc,2)),alpha=0.7) + #Chose colors and size for dots
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=3,force=25, segment.colour="grey",segment.alpha=0.2) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)") + #label Y Axis
  theme #Set Theme
dev.off()


## ADJUSTING FOR HIV: CASES VS SYMPTOMATIC CONTROLS

# Subset OTU table to oral washes
subset.otu.table = subset_samples(otu.table, sample_c %in% c(1))

# Subset OTU table to cases and SCs
subset.otu.table = subset_samples(subset.otu.table, group_c %in% c(0,1))

# Subset OTU table to genus level
subset.genus.table = tax_glom(subset.otu.table, taxrank = "Genus")

# Prune data to less than 100 genera remaining
filtered.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.15 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(filtered.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

# Create relative abundance table
# To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# genus table: pruned.genus.table
# genus relative table: pruned.genus.rel.table
# variable: group

# Convert Phyloseq Object to DESEq object, correcting for any potential confounders
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ hiv + group)

# Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$group <- droplevels(diagdds$group)

# Choose the 'control' in your comparison variable --> positive is upregulated in cases, negative is down-regulated
diagdds$group <- relevel(diagdds$group, ref ="SCs")

# Run the differential Analysis
diagdds<- DESeq(diagdds)

# Output the results from DESeq into a table
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant genes
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

## Tables

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$OTU,rownames(res))
# Replace spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

# Convert resuts table into a data.frame
res <- as.data.frame(res)

#Set names of results table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))

# Make the full trail the first column
res$names <- res$Gene.symbol

res$Gene.symbol <- res$row2

# Get abundance data - use oral.otu.relative.table - THIS CODE ADDS REL ABUNDANCE AS THE DOT SIZE

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
cases.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("1"))
sc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("0"))

# From relative table we should get the mean across the row of the otu table
cases.pruned.genus.rel.table.df <- data.frame(otu_table(cases.pruned.genus.rel.table))
cases.pruned.genus.rel.table.df.meanRA <- rowMeans(cases.pruned.genus.rel.table.df)

sc.pruned.genus.rel.table.df <- data.frame(otu_table(sc.pruned.genus.rel.table))
sc.pruned.genus.rel.table.df.meanRA <- rowMeans(sc.pruned.genus.rel.table.df)

# Subset AND reorder just the OTUs that we have 
cases.pruned.genus.rel.table.df.meanRA.save <- cases.pruned.genus.rel.table.df.meanRA[otu.to.save]
sc.pruned.genus.rel.table.df.meanRA.save <- sc.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res dataframe
res$abundance.cases <- cases.pruned.genus.rel.table.df.meanRA.save
res$abundance.sc <- sc.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need for pathway analysis
res.IPA <- res[,c("Gene.symbol", "logFC", "pvalue", "adj.P.Val")]
res.1 <- res[,c("Gene.symbol", "abundance.cases", "abundance.sc", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Write tables to TXT file
write.table(res.IPA,file="taxa.oral.wash.patients.hiv.IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
write.table(res.1,file="taxa.oral.wash.patients.hiv.abundance.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Volcano plot

# Compute FDR in a log scales
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "red"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "darkgreen"

# Create a Variable for the size of the dots in the Volcano Plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="oral.wash.patients.hiv.volcano.fdr.0.2.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 2000 * res$abundance.cases, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 2000 * res$abundance.sc,2)),alpha=0.7) + #Chose colors and size for dots
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=3,force=25, segment.colour="grey",segment.alpha=0.2) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)") + #label Y Axis
  theme #Set Theme
dev.off()


## ADJUSTING FOR HIV: CASES VS THEIR CLOSE CONTACTS (CCCs)

# Subset OTU table to oral washes
subset.otu.table = subset_samples(otu.table, sample_c %in% c(1))

# Subset OTU table to cases and CCCs
subset.otu.table = subset_samples(subset.otu.table, group_c %in% c(1,2))

# Subset OTU table to genus level
subset.genus.table = tax_glom(subset.otu.table, taxrank = "Genus")

# Prune data to less than 100 genera remaining
filtered.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.15 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(filtered.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

# Create relative abundance table
# To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# genus table: pruned.genus.table
# genus relative table: pruned.genus.rel.table
# variable: group

# Convert Phyloseq Object to DESEq object, correcting for any potential confounders
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ hiv + pairs + group)

# Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$group <- droplevels(diagdds$group)

# Choose the 'control' in your comparison variable 
diagdds$group <- relevel(diagdds$group, ref ="CCCs")

# Run the differential Analysis
diagdds<- DESeq(diagdds)

# Output the results from DESeq into a table
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant genes
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

## Tables

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$OTU,rownames(res))
# Replace spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

# Convert resuts table into a data.frame
res <- as.data.frame(res)

#Set names of results table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))

# Make the full trail the first column
res$names <- res$Gene.symbol

res$Gene.symbol <- res$row2

# Get abundance data - use oral.otu.relative.table - THIS CODE ADDS REL ABUNDANCE AS THE DOT SIZE

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
cases.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("1"))
ccc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("2"))

# From relative table we should get the mean across the row of the otu table
cases.pruned.genus.rel.table.df <- data.frame(otu_table(cases.pruned.genus.rel.table))
cases.pruned.genus.rel.table.df.meanRA <- rowMeans(cases.pruned.genus.rel.table.df)

ccc.pruned.genus.rel.table.df <- data.frame(otu_table(ccc.pruned.genus.rel.table))
ccc.pruned.genus.rel.table.df.meanRA <- rowMeans(ccc.pruned.genus.rel.table.df)

# Subset AND reorder just the OTUs that we have 
cases.pruned.genus.rel.table.df.meanRA.save <- cases.pruned.genus.rel.table.df.meanRA[otu.to.save]
ccc.pruned.genus.rel.table.df.meanRA.save <- ccc.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res dataframe
res$abundance.cases <- cases.pruned.genus.rel.table.df.meanRA.save
res$abundance.ccc <- ccc.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need for pathway analysis
res.IPA <- res[,c("Gene.symbol", "logFC", "pvalue", "adj.P.Val")]
res.1 <- res[,c("Gene.symbol", "abundance.cases", "abundance.ccc", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Write tables to TXT file
write.table(res.IPA,file="taxa.oral.wash.cases.vs.CCCs.hiv.IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
write.table(res.1,file="taxa.oral.wash.cases.vs.CCCs.hiv.abundance.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Volcano plot

# Compute FDR in a log scales
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "red"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "darkgreen"

# Create a Variable for the size of the dots in the Volcano Plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="oral.wash.cases.vs.CCCs.hiv.volcano.fdr.0.2.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 2000 * res$abundance.cases, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 2000 * res$abundance.ccc,2)),alpha=0.7) + #Chose colors and size for dots
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=3,force=25, segment.colour="grey",segment.alpha=0.2) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)") + #label Y Axis
  theme #Set Theme
dev.off()



## ADJUSTING FOR HIV: SYMPTOMATIC CONTROLS (SCs) VS THEIR CLOSE CONTACTS (CCCs)

# Subset OTU table to oral washes
subset.otu.table = subset_samples(otu.table, sample_c %in% c(1))

# Subset OTU table to SCs and CCSCs
subset.otu.table = subset_samples(subset.otu.table, group_c %in% c(0,3))

# Subset OTU table to genus level
subset.genus.table = tax_glom(subset.otu.table, taxrank = "Genus")

# Prune data to less than 100 genera remaining
filtered.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.20 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(filtered.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

# Create relative abundance table
# To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# genus table: pruned.genus.table
# genus relative table: pruned.genus.rel.table
# variable: group

# Convert Phyloseq Object to DESEq object, correcting for any potential confounders
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ hiv + pairs + group)

# Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$group <- droplevels(diagdds$group)

# Choose the 'control' in your comparison variable 
diagdds$group <- relevel(diagdds$group, ref ="CCSCs")

# Run the differential Analysis
diagdds<- DESeq(diagdds)

# Output the results from DESeq into a table
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant genes
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

## TABLES

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$OTU,rownames(res))
# Replace spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

# Convert resuts table into a data.frame
res <- as.data.frame(res)

#Set names of results table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))

# Make the full trail the first column
res$names <- res$Gene.symbol

res$Gene.symbol <- res$row2

# Get abundance data - use oral.otu.relative.table - THIS CODE ADDS REL ABUNDANCE AS THE DOT SIZE

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
sc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("0"))
ccsc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("3"))

# From relative table we should get the mean across the row of the otu table
sc.pruned.genus.rel.table.df <- data.frame(otu_table(sc.pruned.genus.rel.table))
sc.pruned.genus.rel.table.df.meanRA <- rowMeans(sc.pruned.genus.rel.table.df)

ccsc.pruned.genus.rel.table.df <- data.frame(otu_table(ccsc.pruned.genus.rel.table))
ccsc.pruned.genus.rel.table.df.meanRA <- rowMeans(ccsc.pruned.genus.rel.table.df)

# Subset AND reorder just the OTUs that we have 
sc.pruned.genus.rel.table.df.meanRA.save <- sc.pruned.genus.rel.table.df.meanRA[otu.to.save]
ccsc.pruned.genus.rel.table.df.meanRA.save <- ccsc.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res dataframe
res$abundance.sc <- sc.pruned.genus.rel.table.df.meanRA.save
res$abundance.ccsc <- ccsc.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need for pathway analysis
res.IPA <- res[,c("Gene.symbol", "logFC", "pvalue", "adj.P.Val")]
res.1 <- res[,c("Gene.symbol", "abundance.sc", "abundance.ccsc", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Write tables to TXT file
write.table(res.IPA,file="taxa.oral.wash.SCs.vs.CCSCs.hiv.IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
write.table(res.1,file="taxa.oral.wash.SCs.vs.CCSCs.hiv.abundance.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Volcano plot

# Compute FDR in a log scales
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "red"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "darkgreen"

# Create a Variable for the size of the dots in the Volcano Plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="oral.wash.SCs.vs.CCSCs.hiv.volcano.fdr.0.2.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 1000 * res$abundance.sc, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 1000 * res$abundance.ccsc,2)),alpha=0.7) + #Chose colors and size for dots
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=3,force=25, segment.colour="grey",segment.alpha=0.2) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)") + #label Y Axis
  theme #Set Theme
dev.off()


## CLOSE CONTACTS OF CASES (CCCs) VS. CLOSE CONTACTS OF SYMPTOMATIC CONTROLS (CCSCs)

# Subset OTU table to oral washes
subset.otu.table = subset_samples(otu.table, sample_c %in% c(1))

# Subset OTU table to CCCs and CCSCs
subset.otu.table = subset_samples(subset.otu.table, group_c %in% c(2,3))

# Subset OTU table to genus level
subset.genus.table = tax_glom(subset.otu.table, taxrank = "Genus")

# Prune data to less than 100 genera remaining
filtered.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.15 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(filtered.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

# Create relative abundance table
# To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# genus table: pruned.genus.table
# genus relative table: pruned.genus.rel.table
# variable: group

# Convert Phyloseq Object to DESEq object, correcting for any potential confounders
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ group)

# Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$group <- droplevels(diagdds$group)

# Choose the 'control' in your comparison variable 
diagdds$group <- relevel(diagdds$group, ref ="CCSCs")

# Run the differential Analysis
diagdds<- DESeq(diagdds)

# Output the results from DESeq into a table
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant genes
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

## TABLES

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$OTU,rownames(res))
# Replace spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

# Convert resuts table into a data.frame
res <- as.data.frame(res)

#Set names of results table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))

# Make the full trail the first column
res$names <- res$Gene.symbol

res$Gene.symbol <- res$row2

# Get abundance data - use oral.otu.relative.table - THIS CODE ADDS REL ABUNDANCE AS THE DOT SIZE

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
ccc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("2"))
ccsc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("3"))

# From relative table we should get the mean across the row of the otu table
ccc.pruned.genus.rel.table.df <- data.frame(otu_table(ccc.pruned.genus.rel.table))
ccc.pruned.genus.rel.table.df.meanRA <- rowMeans(ccc.pruned.genus.rel.table.df)

ccsc.pruned.genus.rel.table.df <- data.frame(otu_table(ccsc.pruned.genus.rel.table))
ccsc.pruned.genus.rel.table.df.meanRA <- rowMeans(ccsc.pruned.genus.rel.table.df)

# Subset AND reorder just the OTUs that we have 
ccc.pruned.genus.rel.table.df.meanRA.save <- ccc.pruned.genus.rel.table.df.meanRA[otu.to.save]
ccsc.pruned.genus.rel.table.df.meanRA.save <- ccsc.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res dataframe
res$abundance.ccc <- ccc.pruned.genus.rel.table.df.meanRA.save
res$abundance.ccsc <- ccsc.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need for pathway analysis
res.IPA <- res[,c("Gene.symbol", "logFC", "pvalue", "adj.P.Val")]
res.1 <- res[,c("Gene.symbol", "abundance.ccc", "abundance.ccsc", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Write tables to TXT file
write.table(res.IPA,file="taxa.oral.wash.CCCs.vs.CCSCs.IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
write.table(res.1,file="taxa.oral.wash.CCCs.vs.CCSCs.abundance.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Volcano plot

# Compute FDR in a log scales
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "red"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "darkgreen"

# Create a Variable for the size of the dots in the Volcano Plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="oral.wash.CCCs.vs.CCSCs.volcano.fdr.0.2.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 2000 * res$abundance.ccc, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 2000 * res$abundance.ccsc,2)),alpha=0.7) + #Chose colors and size for dots
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=3,force=25, segment.colour="grey",segment.alpha=0.2) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)") + #label Y Axis
  theme #Set Theme
dev.off()

#############################################################################################
#############################################################################################
#############################################################################################

### BACKGROUND CONTROLS

# Subset oral washes with paired background controls
oral.bkg.otu.relative.table = subset_samples(otu.relative.table, subset_paired %in% c(1))

## BETA DIVERSITY

# Create Distance Matrix
vegdist = distance(oral.bkg.otu.relative.table, "wunifrac")

# Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

# Calculate Sample variance for each PC
vars <- apply(CmdScale, 2, var)

# Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

# Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(oral.bkg.otu.relative.table), by = "row.names", all.x = TRUE)

# Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ sample,data= newResults, mean)

# Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="sample",suffixes=c("",".centroid"))

# Plot
pdf("oral.vs.water.wuniF.pdf", height = 10, width = 15)
ggplot(newResults, aes(PC1, PC2, color= sample)) + # Graph PC1 and PC2
  geom_point(size=3,alpha=0.7) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  scale_color_manual(values=c("mediumseagreen", "darkviolet")) + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= sample), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= sample, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=sample, size=10)) + 
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Statistics
adonis(vegdist ~ sample_data(oral.bkg.otu.relative.table)$sample_c)

# Subset induced sputum with paired background controls
sputum.bkg.otu.relative.table = subset_samples(otu.relative.table, subset_paired %in% c(2))

# Create Distance Matrix
vegdist = distance(sputum.bkg.otu.relative.table, "wunifrac")

# Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

# Calculate Sample variance for each PC
vars <- apply(CmdScale, 2, var)

# Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

# Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(sputum.bkg.otu.relative.table), by = "row.names", all.x = TRUE)

# Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ sample,data= newResults, mean)

# Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="sample",suffixes=c("",".centroid"))

# Plot
pdf("sputum.vs.saline.wuniF.pdf", height = 10, width = 15)
ggplot(newResults, aes(PC1, PC2, color= sample)) + # Graph PC1 and PC2
  geom_point(size=3,alpha=0.7) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  scale_color_manual(values=c("gold1", "firebrick1")) + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= sample), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= sample, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=sample, size=10)) + 
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Statistics
adonis(vegdist ~ sample_data(sputum.bkg.otu.relative.table)$sample_c)


#########################################################
#########################################################
#########################################################

### HEATMAPS 

# Load packages
library("ggplot2")
theme_set(theme_bw())
library("RColorBrewer")
library("gplots")
library('ape')
library("plyr")
library("d3heatmap")
library("vegan")
library("Heatplus")
library("igraph")
library('Hmisc')
library("reshape2")
theme_set(theme_bw())

# Subset oral wash and water
oral_Genus.rel.table = subset_samples(Genus.rel.table, subset_paired %in% c(1))
colnames(sample_data(oral_Genus.rel.table))
rownames(sample_data(oral_Genus.rel.table))

# Subset sputum and saline
sputum_Genus.rel.table = subset_samples(Genus.rel.table, subset_paired %in% c(2))
colnames(sample_data(sputum_Genus.rel.table))
rownames(sample_data(sputum_Genus.rel.table))

## Oral wash and water

# Select for genera present in >5% relative abundance in 0.5% of the samples (this approach brings in > 75% of the data in almost all samples) WE LIKE THIS APPROACH
Genus.Rel.wh1 = genefilter_sample(oral_Genus.rel.table, filterfun_sample(function(x) x > 0.05), A = 0.005 * nsamples(oral_Genus.rel.table))
Genus.Rel.table1B = prune_taxa(Genus.Rel.wh1, oral_Genus.rel.table)
colnames(sample_data(Genus.Rel.table1B))
plot_bar(Genus.Rel.table1B, fill="Genus")

# Set data tables  
GenusData <-otu_table(Genus.Rel.table1B) #pruned to selected Genera based on abundance

# Create vector to label by sample type 
sample_data(Genus.Rel.table1B)$sample
SampleVector = sample_data(Genus.Rel.table1B)$sample_c
sample_data(Genus.Rel.table1B)$sample_c

# Duplicate to create a color vector and replace value with color 
# Colorvector can only replace numbers! 
Colorvector <-SampleVector
Colorvector <- replace(Colorvector, which (Colorvector == "1"), "mediumseagreen")
Colorvector <- replace(Colorvector, which (Colorvector == "5"), "darkviolet")

# Cluster Bray Heatmap
# Cluster Genuses(row)
GenusData.Bray.dist <-vegdist(GenusData, method = "bray")
Genus.Bray.clus <-hclust(GenusData.Bray.dist, "aver")

# Cluster samples(Col) 
Bray.dist = distance(GenusData, method="bray")  
cluster.Bray = hclust(Bray.dist, "aver")

# Change the names for genuses that are labelled as "g__" --> Come back to this
tax_table(Genus.Rel.table1B)
Genus.Rel.table1B.New.Names = prune_taxa(tail(names(sort(taxa_sums(Genus.Rel.table1B))), ntaxa(Genus.Rel.table1B)), Genus.Rel.table1B)

# Add a new rank, Strain, with the Genus ids
tax_table(Genus.Rel.table1B.New.Names) <- cbind(tax_table(Genus.Rel.table1B.New.Names), Strain=taxa_names(Genus.Rel.table1B.New.Names))

# Define the ranks you want to include
myranks = c("Family", "Genus")
mylabels = apply(tax_table(Genus.Rel.table1B.New.Names)[, myranks], 1, paste, sep="", collapse="_")

# Add concatenated labels as a new rank after strain
tax_table(Genus.Rel.table1B.New.Names) <- cbind(tax_table(Genus.Rel.table1B.New.Names), catglab=mylabels)
tax_table(Genus.Rel.table1B.New.Names)   

# Plot heatmap with dendograms
mypalette <- colorRampPalette(c('#ffffff','#4169E1','#0000CD'))

# Save to PDF
pdf("heatmap.oral.vs.water.pdf", height = 15, width = 30)
heatmap.2(GenusData, margins = c(10,20), 
          density.info = "none",
          trace = "none",
          keysize = 0.75,
          key.title = "Relative abundance",
          offsetRow = 1, offsetCol = 1,
          dendrogram = "both",
          Rowv = as.dendrogram(Genus.Bray.clus),
          Colv = as.dendrogram(cluster.Bray),
          labRow=tax_table(Genus.Rel.table1B.New.Names)[,"catglab"],
          cexRow = 1.5,
          labCol = sample_data(Genus.Rel.table1B)$heatmap,
          cexCol = 1.5,
          col = mypalette(17),
          symm=F,symkey=F,symbreaks=T, scale="none",
          breaks =c(seq(0,.1,length=10),seq(.11,0.3,length=4),seq(0.31,.7,length=4)),
          ColSideColors=Colorvector,
          main = "Oral wash vs Water",
)
dev.off()

## Induced sputum and saline

# Select for genera present in >5% relative abundance in 0.5% of the samples (this approach brings in > 75% of the data in almost all samples) WE LIKE THIS APPROACH
Genus.Rel.wh1 = genefilter_sample(sputum_Genus.rel.table, filterfun_sample(function(x) x > 0.05), A = 0.005 * nsamples(sputum_Genus.rel.table))
Genus.Rel.table1B = prune_taxa(Genus.Rel.wh1, sputum_Genus.rel.table)
colnames(sample_data(Genus.Rel.table1B))
plot_bar(Genus.Rel.table1B, fill="Genus")

# Set data tables  
GenusData <-otu_table(Genus.Rel.table1B) 

# Create vector to label by sample type 
sample_data(Genus.Rel.table1B)$sample
SampleVector = sample_data(Genus.Rel.table1B)$sample_c
sample_data(Genus.Rel.table1B)$sample_c

# Duplicate to create a color vector and replace value with color 
# Colorvector can only replace numbers! 
Colorvector <-SampleVector
Colorvector <- replace(Colorvector, which (Colorvector == "4"), "gold1")
Colorvector <- replace(Colorvector, which (Colorvector == "2"), "firebrick1")

# Cluster Bray Heatmap
# Cluster Genuses(row)
GenusData.Bray.dist <-vegdist(GenusData, method = "bray")
Genus.Bray.clus <-hclust(GenusData.Bray.dist, "aver")

# Cluster samples(Col) 
Bray.dist = distance(GenusData, method="bray")  
cluster.Bray = hclust(Bray.dist, "aver")

# Change the names for genuses that are labelled as "g__" --> Come back to this
tax_table(Genus.Rel.table1B)
Genus.Rel.table1B.New.Names = prune_taxa(tail(names(sort(taxa_sums(Genus.Rel.table1B))), ntaxa(Genus.Rel.table1B)), Genus.Rel.table1B)

# Add a new rank, Strain, with the Genus ids
tax_table(Genus.Rel.table1B.New.Names) <- cbind(tax_table(Genus.Rel.table1B.New.Names), Strain=taxa_names(Genus.Rel.table1B.New.Names))

# Define the ranks you want to include
myranks = c("Family", "Genus")
mylabels = apply(tax_table(Genus.Rel.table1B.New.Names)[, myranks], 1, paste, sep="", collapse="_")

# Add concatenated labels as a new rank after strain
tax_table(Genus.Rel.table1B.New.Names) <- cbind(tax_table(Genus.Rel.table1B.New.Names), catglab=mylabels)
tax_table(Genus.Rel.table1B.New.Names)  

# Plot heatmap with dendograms
mypalette <- colorRampPalette(c('#ffffff','#4169E1','#0000CD'))
pdf("heatmap.sputum.vs.saline.pdf", height = 15, width = 30)
heatmap.2(GenusData, margins = c(10,20), 
          density.info = "none",
          trace = "none",
          keysize = 0.75,
          key.title = "Relative abundance",
          offsetRow = 1, offsetCol = 1,
          dendrogram = "both",
          Rowv = as.dendrogram(Genus.Bray.clus),
          Colv = as.dendrogram(cluster.Bray),
          labRow=tax_table(Genus.Rel.table1B.New.Names)[,"catglab"],
          cexRow = 1.5,
          labCol = sample_data(Genus.Rel.table1B)$heatmap,
          cexCol = 1.5,
          col = mypalette(17),
          symm=F,symkey=F,symbreaks=T, scale="none",
          breaks =c(seq(0,.1,length=10),seq(.11,0.3,length=4),seq(0.31,.7,length=4)),
          ColSideColors=Colorvector,
          main = "Induced sputum vs saline",
)
dev.off()

#########################################################
#########################################################
#########################################################

### DIFFERENTIAL ABUNDANCE ANALYSIS WITH DESEQ2

# Set Theme For Figures
theme <-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),
              axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), 
              legend.position="none")

# Set alpha
alpha <- 0.20

## Oral wash vs water

# Subset otu.table to oral washes and water
subset.otu.table = subset_samples(otu.table, subset_paired %in% c(1))

# Subset otu table to genus level
subset.genus.table = tax_glom(subset.otu.table, taxrank = "Genus")

# Prune data to less than 100 genera remaining
filtered.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.20 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(filtered.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

# To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# Convert Phyloseq Object to DESEq, correcting for any potential confounders
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ sample)

# Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$sample <- droplevels(diagdds$sample)

# Choose which is the 'control' in your comparison variable 
diagdds$sample <- relevel(diagdds$sample, ref ="water")

# Run the differential Analysis
diagdds<- DESeq(diagdds)

# Output the results from DESEQ into a table
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 Significant Genes
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

# Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with Taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$OTU,rownames(res))

# Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

# Convert Resuts table into a data.frame
res <- as.data.frame(res)

# Set Names of Results Table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))

# Make the full trail the First Column
res$names <- res$Gene.symbol

res$Gene.symbol <- res$row2

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
oral.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, sample_c %in% c("1"))
water.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, sample_c %in% c("5"))

# From relative table we should get the mean across the row of the otu table
oral.pruned.genus.rel.table.df <- data.frame(otu_table(oral.pruned.genus.rel.table))
oral.pruned.genus.rel.table.df.meanRA <- rowMeans(oral.pruned.genus.rel.table.df)

water.pruned.genus.rel.table.df <- data.frame(otu_table(water.pruned.genus.rel.table))
water.pruned.genus.rel.table.df.meanRA <- rowMeans(water.pruned.genus.rel.table.df)

# Subset AND reorder just the otus that we have 
oral.pruned.genus.rel.table.df.meanRA.save <- oral.pruned.genus.rel.table.df.meanRA[otu.to.save]
water.pruned.genus.rel.table.df.meanRA.save <- water.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundnace data for the res dataframe
res$abundance.oral <- oral.pruned.genus.rel.table.df.meanRA.save
res$abundance.water <- water.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need for pathway analysis
res.1 <- res[,c("Gene.symbol", "abundance.oral", "abundance.water", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Write Tables to TXT file
write.table(res.1,file="taxa.oral.wash.vs.water.abundance.txt.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Compute FDR in a log scales
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "mediumseagreen"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "darkviolet"

# Create a Variable for the size of the dots in the Volcano Plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="oral.wash.vs.water.volcano.fdr.0.2.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 100 * res$abundance.oral, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 100 * res$abundance.water,2)),alpha=0.7) + #Chose Colors and size for dots
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=3,force=25, segment.colour="grey",segment.alpha=0.2) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)") + #label Y Axis
  theme #Set Theme
dev.off()

## Induced sputum vs saline

# Subset otu.table to induced sputum and saline
subset.otu.table = subset_samples(otu.table, subset_paired %in% c(2))

# Subset otu table to genus level
subset.genus.table = tax_glom(subset.otu.table, taxrank = "Genus")

# Prune data to less than 100 genera remaining
filtered.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.20 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(filtered.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

# To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# Convert Phyloseq Object to DESEq, correcting for any potential confounders
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ sample)

# Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$sample <- droplevels(diagdds$sample)

# Choose which is the 'control' in your comparison variable 
diagdds$sample <- relevel(diagdds$sample, ref ="saline")

# Run the differential Analysis
diagdds<- DESeq(diagdds)

# Output the results from DESEQ into a table
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 Significant Genes
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

# Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with Taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$OTU,rownames(res))

# Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

# Convert Resuts table into a data.frame
res <- as.data.frame(res)

# Set Names of Results Table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))

# Make the full trail the First Column
res$names <- res$Gene.symbol

res$Gene.symbol <- res$row2

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
sputum.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, sample_c %in% c("2"))
saline.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, sample_c %in% c("4"))

# From relative table we should get the mean across the row of the otu table
sputum.pruned.genus.rel.table.df <- data.frame(otu_table(sputum.pruned.genus.rel.table))
sputum.pruned.genus.rel.table.df.meanRA <- rowMeans(sputum.pruned.genus.rel.table.df)

saline.pruned.genus.rel.table.df <- data.frame(otu_table(saline.pruned.genus.rel.table))
saline.pruned.genus.rel.table.df.meanRA <- rowMeans(saline.pruned.genus.rel.table.df)

# Subset AND reorder just the otus that we have 
sputum.pruned.genus.rel.table.df.meanRA.save <- sputum.pruned.genus.rel.table.df.meanRA[otu.to.save]
saline.pruned.genus.rel.table.df.meanRA.save <- saline.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundnace data for the res dataframe
res$abundance.sputum <- sputum.pruned.genus.rel.table.df.meanRA.save
res$abundance.saline <- saline.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need for pathway analysis
res.1 <- res[,c("Gene.symbol", "abundance.sputum", "abundance.saline", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Write Tables to TXT file
write.table(res.1,file="taxa.induced.sputum.vs.saline.abundance.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Compute FDR in a log scales
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "firebrick1"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "gold1"

# Create a Variable for the size of the dots in the Volcano Plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="induced.sputum.vs.saline.volcano.fdr.0.2.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 100 * res$abundance.sputum, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 100 * res$abundance.saline,2)),alpha=0.7) + #Chose Colors and size for dots
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=3,force=25, segment.colour="grey",segment.alpha=0.2) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)") + #label Y Axis
  theme #Set Theme
dev.off()

#########################################################
#########################################################
#########################################################

### HEATMAP: ORAL WASHES VS INDUCED SPUTUM

# Subset oral washes and induced sputum
oral_sputum_Genus.rel.table = subset_samples(Genus.rel.table, sample_c %in% c(1,2))
rownames(sample_data(oral_sputum_Genus.rel.table))

# Select for genera present in >5% relative abundance in 0.5% of the samples (this approach brings in > 75% of the data in almost all samples) WE LIKE THIS APPROACH
Genus.Rel.wh1 = genefilter_sample(oral_sputum_Genus.rel.table, filterfun_sample(function(x) x > 0.05), A = 0.005 * nsamples(oral_sputum_Genus.rel.table))
Genus.Rel.table1B = prune_taxa(Genus.Rel.wh1, oral_sputum_Genus.rel.table)
colnames(sample_data(Genus.Rel.table1B))
plot_bar(Genus.Rel.table1B, fill="Genus")

# Set data tables  
GenusData <-otu_table(Genus.Rel.table1B) 

# Create vector to label by sample type
sample_data(Genus.Rel.table1B)$sample
SampleVector = sample_data(Genus.Rel.table1B)$sample_c
sample_data(Genus.Rel.table1B)$sample_c

# Duplicate to create a color vector and replace value w/ color 
# Colorvector can only replace numbers! 
Colorvector <-SampleVector
Colorvector <- replace(Colorvector, which (Colorvector == "1"), "mediumseagreen")
Colorvector <- replace(Colorvector, which (Colorvector == "2"), "firebrick1")

# Cluster Bray Heatmap
# Cluster Genuses(row)
GenusData.Bray.dist <-vegdist(GenusData, method = "bray")
Genus.Bray.clus <-hclust(GenusData.Bray.dist, "aver")

# Cluster samples(Col) 
Bray.dist = distance(GenusData, method="bray")  
cluster.Bray = hclust(Bray.dist, "aver")

# Change the names for genuses that are labelled as "g__"
tax_table(Genus.Rel.table1B)
Genus.Rel.table1B.New.Names = prune_taxa(tail(names(sort(taxa_sums(Genus.Rel.table1B))), ntaxa(Genus.Rel.table1B)), Genus.Rel.table1B)

# Add a new rank, Strain, with the Genus ids
tax_table(Genus.Rel.table1B.New.Names) <- cbind(tax_table(Genus.Rel.table1B.New.Names), Strain=taxa_names(Genus.Rel.table1B.New.Names))

# Define the ranks you want to include
myranks = c("Family", "Genus")
mylabels = apply(tax_table(Genus.Rel.table1B.New.Names)[, myranks], 1, paste, sep="", collapse="_")

# Add concatenated labels as a new rank after strain
tax_table(Genus.Rel.table1B.New.Names) <- cbind(tax_table(Genus.Rel.table1B.New.Names), catglab=mylabels)
tax_table(Genus.Rel.table1B.New.Names)   ###come back to this -> try output to csv

# Plot Heat map with dendograms
mypalette <- colorRampPalette(c('#ffffff','#4169E1','#0000CD'))

pdf("heatmap.oral.wash.vs.induced.sputum.pdf", height = 15, width = 30)
heatmap.2(GenusData, margins = c(10,20), 
          density.info = "none",
          trace = "none",
          keysize = 0.75,
          key.title = "Relative abundance",
          offsetRow = 1, offsetCol = 1,
          dendrogram = "both",
          Rowv = as.dendrogram(Genus.Bray.clus),
          Colv = as.dendrogram(cluster.Bray),
          labRow=tax_table(Genus.Rel.table1B.New.Names)[,"catglab"],
          cexRow = 1.5,
          labCol = sample_data(Genus.Rel.table1B)$heatmap,
          cexCol = 1.5,
          col = mypalette(17),
          symm=F,symkey=F,symbreaks=T, scale="none",
          breaks =c(seq(0,.1,length=10),seq(.11,0.3,length=4),seq(0.31,.7,length=4)),
          ColSideColors=Colorvector,
          main = "Oral wash vs. Induced sputum",
)
dev.off()


#############################################################################################
#############################################################################################
#############################################################################################

### INDUCED SPUTUM

# Subset induced sputum
sputum.otu.relative.table = subset_samples(otu.relative.table, sample_c %in% c(2))

## BETA DIVERSITY

# Create Distance Matrix
vegdist = distance(sputum.otu.relative.table, "wunifrac")

# Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

# Calculate Sample variance for each PC
vars <- apply(CmdScale, 2, var)

# Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

# Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(sputum.otu.relative.table), by = "row.names", all.x = TRUE)

# Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ group,data= newResults, mean)

# Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="group",suffixes=c("",".centroid"))

pdf("sputum.wuniF.pdf", height = 10, width = 15)
ggplot(newResults, aes(PC1, PC2, color= group)) +
  geom_point(size=3,alpha=0.7) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values=c("red", "purple", "cyan", "forestgreen", "grey")) + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= group), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= group, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=group, size=10)) +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Subset cases and SCs
patients.sputum.otu.relative.table = subset_samples(sputum.otu.relative.table, group_c %in% c(0,1))
# Create distance object
patients = distance(patients.sputum.otu.relative.table, "wunifrac")
# Statistics
adonis(patients ~ sample_data(patients.sputum.otu.relative.table)$group)

# Subset cases and CCCs
cases.ccc.sputum.otu.relative.table = subset_samples(sputum.otu.relative.table, group_c %in% c(1,2))
# Create distance object
cases = distance(cases.ccc.sputum.otu.relative.table, "wunifrac")
# Statistics
adonis(cases ~ sample_data(cases.ccc.sputum.otu.relative.table)$group)

# Subset SCs and CCSCs
sc.ccs.sputum.otu.relative.table = subset_samples(sputum.otu.relative.table, group_c %in% c(0,3))
# Create distance object
sc = distance(sc.ccs.sputum.otu.relative.table, "wunifrac")
# Statistics
adonis(sc ~ sample_data(sc.ccs.sputum.otu.relative.table)$group)

#########################################################
#########################################################
#########################################################

### BETA DIVERSITY IN EACH GROUP BY HIV STATUS

## Cases
# Subset cases
cases.sputum.otu.relative.table = subset_samples(sputum.otu.relative.table, group_c %in% c(1))

# Create Distance Matrix
vegdist = distance(cases.sputum.otu.relative.table, "wunifrac")

# Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

# Calculate Sample variance for each PC
vars <- apply(CmdScale, 2, var)

# Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

# Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(cases.sputum.otu.relative.table), by = "row.names", all.x = TRUE)

# Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ hiv,data= newResults, mean)

# Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="hiv",suffixes=c("",".centroid"))

pdf("sputum.wunif.hiv.cases.pdf", height = 10, width = 10)
ggplot(newResults, aes(PC1, PC2, color= hiv)) + # Graph PC1 and PC2
  geom_point(size=3,alpha=0.7) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  scale_color_manual(values=c("orange","red", "grey")) + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= hiv), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= hiv, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=hiv, size=10)) +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Statistics
adonis(vegdist ~ sample_data(cases.sputum.otu.relative.table)$hiv)


## Symptomatic controls (SCs)
# Subset SCs
sc.sputum.otu.relative.table = subset_samples(sputum.otu.relative.table, group_c %in% c(0))

# Create Distance Matrix
vegdist = distance(sc.sputum.otu.relative.table, "wunifrac")

# Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

# Calculate Sample variance for each PC
vars <- apply(CmdScale, 2, var)

# Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

# Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(sc.sputum.otu.relative.table), by = "row.names", all.x = TRUE)

# Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ hiv,data= newResults, mean)

# Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="hiv",suffixes=c("",".centroid"))

# Plot
pdf("sputum.wunif.hiv.SCs.pdf", height = 10, width = 10)
ggplot(newResults, aes(PC1, PC2, color= hiv)) + # Graph PC1 and PC2
  geom_point(size=3,alpha=0.7) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  scale_color_manual(values=c("orange","red", "grey")) + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= hiv), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= hiv, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=hiv, size=10)) +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Statistics
adonis(vegdist ~ sample_data(sc.sputum.otu.relative.table)$hiv)

## Close contacts of cases (CCCs)
# Subset CCCs
ccc.sputum.otu.relative.table = subset_samples(sputum.otu.relative.table, group_c %in% c(2))

#subset CCCs with only known hiv status
ccc.sputum.otu.relative.table = subset_samples(ccc.sputum.otu.relative.table, known_hiv_status %in% c(1))

# Create Distance Matrix
vegdist = distance(ccc.sputum.otu.relative.table, "wunifrac")

# Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

# Calculate Sample variance for each PC
vars <- apply(CmdScale, 2, var)

# Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

# Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(ccc.sputum.otu.relative.table), by = "row.names", all.x = TRUE)

# Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ hiv,data= newResults, mean)

# Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="hiv",suffixes=c("",".centroid"))

# Plot
pdf("sputum.wunif.hiv.CCCs.pdf", height = 10, width = 10)
ggplot(newResults, aes(PC1, PC2, color= hiv)) + # Graph PC1 and PC2
  geom_point(size=3,alpha=0.7) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  scale_color_manual(values=c("orange","red", "grey")) + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= hiv), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= hiv, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=hiv, size=10)) +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Statistics
adonis(vegdist ~ sample_data(ccc.sputum.otu.relative.table)$hiv)

## Close contacts of symptomatic controls (CCSCs)
# Subset CCSCs
ccsc.sputum.otu.relative.table = subset_samples(sputum.otu.relative.table, group_c %in% c(3))

#subset CCSCs with only known hiv status
ccsc.sputum.otu.relative.table = subset_samples(ccsc.sputum.otu.relative.table, known_hiv_status %in% c(1))

# Create Distance Matrix
vegdist = distance(ccsc.sputum.otu.relative.table, "wunifrac")

# Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

# Calculate Sample variance for each PC
vars <- apply(CmdScale, 2, var)

# Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

# Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(ccsc.sputum.otu.relative.table), by = "row.names", all.x = TRUE)

# Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ hiv,data= newResults, mean)

# Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="hiv",suffixes=c("",".centroid"))

# Plot
pdf("sputum.wunif.hiv.CCSCs.pdf", height = 10, width = 10)
ggplot(newResults, aes(PC1, PC2, color= hiv)) + # Graph PC1 and PC2
  geom_point(size=3,alpha=0.7) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  scale_color_manual(values=c("orange","red", "grey")) + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= hiv), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= hiv, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=hiv, size=10)) +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

#Statistics not done because HIV-positive group only had one patient

#########################################################
#########################################################
#########################################################

### BETA DIVERSITY IN CLOSE CONTACTS OF CASES (CCCs) VS. CLOSE CONTACTS OF SYMPTOMATIC CONTROLS (CCSCs)

# Subset CCCs and CCSCs
contacts.sputum.otu.relative.table = subset_samples(sputum.otu.relative.table, group_c %in% c(2,3))
rownames(sample_data(contacts.sputum.otu.relative.table))

# Create Distance Matrix
vegdist = distance(contacts.sputum.otu.relative.table, "wunifrac")

# Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

# Calculate Sample variance for each PC
vars <- apply(CmdScale, 2, var)

# Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#  PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(contacts.sputum.otu.relative.table), by = "row.names", all.x = TRUE)

# Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ group,data= newResults, mean)

# Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="group",suffixes=c("",".centroid"))

# Plot
pdf("sputum.wuniF.contacts.pdf", height = 10, width = 15)
ggplot(newResults, aes(PC1, PC2, color= group)) + # Graph PC1 and PC2
  geom_point(size=3,alpha=0.7) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  scale_color_manual(values=c("purple", "cyan")) + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= group), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= group, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=group, size=10)) + 
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Statistics
adonis(vegdist ~ sample_data(contacts.sputum.otu.relative.table)$group)

#########################################################
#########################################################
#########################################################

### DIFFERENTIAL ABUNDANCE ANALYSIS WITH DESEQ2

# Set Theme For Figures
theme <-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),
              axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), 
              legend.position="none")

# Load file
#load(file="all.samples.RData")

# Set alpha
alpha <- 0.20


## CASES VS SYMPTOMATIC CONTROLS (SCs)

# Subset OTU table to induced sputum
subset.otu.table = subset_samples(otu.table, sample_c %in% c(2))

# Subset OTU table to cases and SCs
subset.otu.table = subset_samples(subset.otu.table, group_c %in% c(0,1))

# Subset OTU table to genus level
subset.genus.table = tax_glom(subset.otu.table, taxrank = "Genus")

# Prune data to less than 100 genera remaining
filtered.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.25 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(filtered.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

# Create relative abundance table
# To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# genus table: pruned.genus.table
# genus relative table: pruned.genus.rel.table
# variable: group

# Convert Phyloseq Object to DESEq object, correcting for any potential confounders
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ group)

# Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$group <- droplevels(diagdds$group)

# Choose the 'control' in your comparison variable --> positive is upregulated in cases, negative is down-regulated
diagdds$group <- relevel(diagdds$group, ref ="SCs")

# Run the differential Analysis
diagdds<- DESeq(diagdds)

# Output the results from DESeq into a table
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant genes
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

# Tables

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$OTU,rownames(res))
# Replace spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

# Convert resuts table into a data.frame
res <- as.data.frame(res)

#Set names of results table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))

# Make the full trail the first column
res$names <- res$Gene.symbol

res$Gene.symbol <- res$row2

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
cases.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("1"))
sc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("0"))

# From relative table we should get the mean across the row of the otu table
cases.pruned.genus.rel.table.df <- data.frame(otu_table(cases.pruned.genus.rel.table))
cases.pruned.genus.rel.table.df.meanRA <- rowMeans(cases.pruned.genus.rel.table.df)

sc.pruned.genus.rel.table.df <- data.frame(otu_table(sc.pruned.genus.rel.table))
sc.pruned.genus.rel.table.df.meanRA <- rowMeans(sc.pruned.genus.rel.table.df)

# Subset AND reorder just the OTUs that we have 
cases.pruned.genus.rel.table.df.meanRA.save <- cases.pruned.genus.rel.table.df.meanRA[otu.to.save]
sc.pruned.genus.rel.table.df.meanRA.save <- sc.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res dataframe
res$abundance.cases <- cases.pruned.genus.rel.table.df.meanRA.save
res$abundance.sc <- sc.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need for pathway analysis
res.IPA <- res[,c("Gene.symbol", "logFC", "pvalue", "adj.P.Val")]
res.1 <- res[,c("Gene.symbol", "abundance.cases", "abundance.sc", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Write tables to TXT file
write.table(res.IPA,file="taxa.sputum.patients.IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
write.table(res.1,file="taxa.sputum.patients.abundance.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Compute FDR in a log scales
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "red"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "darkgreen"

# Create a Variable for the size of the dots in the Volcano Plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="sputum.patients.volcano.fdr.0.2.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 1000 * res$abundance.cases, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 1000 * res$abundance.sc,2)),alpha=0.7) + #Chose colors and size for dots
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=3,force=25, segment.colour="grey",segment.alpha=0.2) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)") + #label Y Axis
  theme #Set Theme
dev.off()


## CASES VS THEIR CLOSE CONTACTS (CCCs)

# Subset OTU table to induced sputum
subset.otu.table = subset_samples(otu.table, sample_c %in% c(2))

# Subset OTU table to cases and CCCs
subset.otu.table = subset_samples(subset.otu.table, group_c %in% c(1,2))

# Subset OTU table to genus level
subset.genus.table = tax_glom(subset.otu.table, taxrank = "Genus")

# Prune data to less than 100 genera remaining
filtered.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.25 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(filtered.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

# Create relative abundance table
# To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# genus table: pruned.genus.table
# genus relative table: pruned.genus.rel.table
# variable: group

# Convert Phyloseq Object to DESEq object, correcting for any potential confounders
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ pairs + group)

# Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$group <- droplevels(diagdds$group)

# Choose the 'control' in your comparison variable 
diagdds$group <- relevel(diagdds$group, ref ="CCCs")

# Run the differential Analysis
diagdds<- DESeq(diagdds)

# Output the results from DESeq into a table
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant genes
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

## Tables

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$OTU,rownames(res))
# Replace spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

# Convert resuts table into a data.frame
res <- as.data.frame(res)

#Set names of results table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))

# Make the full trail the first column
res$names <- res$Gene.symbol

res$Gene.symbol <- res$row2

# Get abundance data - use sputum.otu.relative.table - THIS CODE ADDS REL ABUNDANCE AS THE DOT SIZE

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
cases.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("1"))
ccc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("2"))

# From relative table we should get the mean across the row of the otu table
cases.pruned.genus.rel.table.df <- data.frame(otu_table(cases.pruned.genus.rel.table))
cases.pruned.genus.rel.table.df.meanRA <- rowMeans(cases.pruned.genus.rel.table.df)

ccc.pruned.genus.rel.table.df <- data.frame(otu_table(ccc.pruned.genus.rel.table))
ccc.pruned.genus.rel.table.df.meanRA <- rowMeans(ccc.pruned.genus.rel.table.df)

# Subset AND reorder just the OTUs that we have 
cases.pruned.genus.rel.table.df.meanRA.save <- cases.pruned.genus.rel.table.df.meanRA[otu.to.save]
ccc.pruned.genus.rel.table.df.meanRA.save <- ccc.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res dataframe
res$abundance.cases <- cases.pruned.genus.rel.table.df.meanRA.save
res$abundance.ccc <- ccc.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need for pathway analysis
res.IPA <- res[,c("Gene.symbol", "logFC", "pvalue", "adj.P.Val")]
res.1 <- res[,c("Gene.symbol", "abundance.cases", "abundance.ccc", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Write tables to TXT file
write.table(res.IPA,file="taxa.sputum.cases.vs.CCCs.IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
write.table(res.1,file="taxa.sputum.cases.vs.CCCs.abundance.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Volcano plot

# Compute FDR in a log scales
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "red"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "darkgreen"

# Create a Variable for the size of the dots in the Volcano Plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="sputum.cases.vs.CCCs.volcano.fdr.0.2.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 1000 * res$abundance.cases, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 1000 * res$abundance.ccc,2)),alpha=0.7) + #Chose colors and size for dots
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=3,force=25, segment.colour="grey",segment.alpha=0.2) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)") + #label Y Axis
  theme #Set Theme
dev.off()

## SYMPTOMATIC CONTROLS (SCs) VS THEIR CLOSE CONTACTS (CCCs)

# Subset OTU table to induced sputum
subset.otu.table = subset_samples(otu.table, sample_c %in% c(2))

# Subset OTU table to SCs and CCSCs
subset.otu.table = subset_samples(subset.otu.table, group_c %in% c(0,3))

# Subset OTU table to genus level
subset.genus.table = tax_glom(subset.otu.table, taxrank = "Genus")

# Prune data to less than 100 genera remaining
filtered.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.25 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(filtered.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

# Create relative abundance table
# To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# genus table: pruned.genus.table
# genus relative table: pruned.genus.rel.table
# variable: group

# Convert Phyloseq Object to DESEq object, correcting for any potential confounders
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ pairs + group)

# Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$group <- droplevels(diagdds$group)

# Choose the 'control' in your comparison variable 
diagdds$group <- relevel(diagdds$group, ref ="CCSCs")

# Run the differential Analysis
diagdds<- DESeq(diagdds)

# Output the results from DESeq into a table
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant genes
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

## TABLES

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$OTU,rownames(res))
# Replace spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

# Convert resuts table into a data.frame
res <- as.data.frame(res)

#Set names of results table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))

# Make the full trail the first column
res$names <- res$Gene.symbol

res$Gene.symbol <- res$row2

# Get abundance data - use sputum.otu.relative.table - THIS CODE ADDS REL ABUNDANCE AS THE DOT SIZE

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
sc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("0"))
ccsc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("3"))

# From relative table we should get the mean across the row of the otu table
sc.pruned.genus.rel.table.df <- data.frame(otu_table(sc.pruned.genus.rel.table))
sc.pruned.genus.rel.table.df.meanRA <- rowMeans(sc.pruned.genus.rel.table.df)

ccsc.pruned.genus.rel.table.df <- data.frame(otu_table(ccsc.pruned.genus.rel.table))
ccsc.pruned.genus.rel.table.df.meanRA <- rowMeans(ccsc.pruned.genus.rel.table.df)

# Subset AND reorder just the OTUs that we have 
sc.pruned.genus.rel.table.df.meanRA.save <- sc.pruned.genus.rel.table.df.meanRA[otu.to.save]
ccsc.pruned.genus.rel.table.df.meanRA.save <- ccsc.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res dataframe
res$abundance.sc <- sc.pruned.genus.rel.table.df.meanRA.save
res$abundance.ccsc <- ccsc.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need for pathway analysis
res.IPA <- res[,c("Gene.symbol", "logFC", "pvalue", "adj.P.Val")]
res.1 <- res[,c("Gene.symbol", "abundance.sc", "abundance.ccsc", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Write tables to TXT file
write.table(res.IPA,file="taxa.sputum.SCs.vs.CCSCs.IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
write.table(res.1,file="taxa.sputum.SCs.vs.CCSCs.abundance.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

## Volcano plot

# Compute FDR in a log scales
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "red"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "darkgreen"

# Create a Variable for the size of the dots in the Volcano Plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="sputum.SCs.vs.CCSCs.volcano.fdr.0.2.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 500 * res$abundance.sc, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 500 * res$abundance.ccsc,2)),alpha=0.7) + #Chose colors and size for dots
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=3,force=25, segment.colour="grey",segment.alpha=0.2) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)") + #label Y Axis
  theme #Set Theme
dev.off()


## ADJUSTING FOR HIV: CASES VS SYMPTOMATIC CONTROLS

# Subset OTU table to induced sputum
subset.otu.table = subset_samples(otu.table, sample_c %in% c(2))

# Subset OTU table to cases and SCs
subset.otu.table = subset_samples(subset.otu.table, group_c %in% c(0,1))

# Subset OTU table to genus level
subset.genus.table = tax_glom(subset.otu.table, taxrank = "Genus")

# Prune data to less than 100 genera remaining
filtered.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.25 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(filtered.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

# Create relative abundance table
# To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# genus table: pruned.genus.table
# genus relative table: pruned.genus.rel.table
# variable: group

# Convert Phyloseq Object to DESEq object, correcting for any potential confounders
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ hiv + group)

# Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$group <- droplevels(diagdds$group)

# Choose the 'control' in your comparison variable --> positive is upregulated in cases, negative is down-regulated
diagdds$group <- relevel(diagdds$group, ref ="SCs")

# Run the differential Analysis
diagdds<- DESeq(diagdds)

# Output the results from DESeq into a table
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant genes
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

## Tables

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$OTU,rownames(res))
# Replace spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

# Convert resuts table into a data.frame
res <- as.data.frame(res)

#Set names of results table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))

# Make the full trail the first column
res$names <- res$Gene.symbol

res$Gene.symbol <- res$row2

# Get abundance data - use sputum.otu.relative.table - THIS CODE ADDS REL ABUNDANCE AS THE DOT SIZE

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
cases.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("1"))
sc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("0"))

# From relative table we should get the mean across the row of the otu table
cases.pruned.genus.rel.table.df <- data.frame(otu_table(cases.pruned.genus.rel.table))
cases.pruned.genus.rel.table.df.meanRA <- rowMeans(cases.pruned.genus.rel.table.df)

sc.pruned.genus.rel.table.df <- data.frame(otu_table(sc.pruned.genus.rel.table))
sc.pruned.genus.rel.table.df.meanRA <- rowMeans(sc.pruned.genus.rel.table.df)

# Subset AND reorder just the OTUs that we have 
cases.pruned.genus.rel.table.df.meanRA.save <- cases.pruned.genus.rel.table.df.meanRA[otu.to.save]
sc.pruned.genus.rel.table.df.meanRA.save <- sc.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res dataframe
res$abundance.cases <- cases.pruned.genus.rel.table.df.meanRA.save
res$abundance.sc <- sc.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need for pathway analysis
res.IPA <- res[,c("Gene.symbol", "logFC", "pvalue", "adj.P.Val")]
res.1 <- res[,c("Gene.symbol", "abundance.cases", "abundance.sc", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Write tables to TXT file
write.table(res.IPA,file="taxa.sputum.patients.hiv.IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
write.table(res.1,file="taxa.sputum.patients.hiv.abundance.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Volcano plot

# Compute FDR in a log scales
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "red"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "darkgreen"

# Create a Variable for the size of the dots in the Volcano Plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="sputum.patients.hiv.volcano.fdr.0.2.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 1000 * res$abundance.cases, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 1000 * res$abundance.sc,2)),alpha=0.7) + #Chose colors and size for dots
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=3,force=25, segment.colour="grey",segment.alpha=0.2) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)") + #label Y Axis
  theme #Set Theme
dev.off()


## ADJUSTING FOR HIV: CASES VS THEIR CLOSE CONTACTS (CCCs)

# Subset OTU table to induced sputum
subset.otu.table = subset_samples(otu.table, sample_c %in% c(2))

# Subset OTU table to cases and CCCs
subset.otu.table = subset_samples(subset.otu.table, group_c %in% c(1,2))

# Subset OTU table to genus level
subset.genus.table = tax_glom(subset.otu.table, taxrank = "Genus")

# Prune data to less than 100 genera remaining
filtered.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.25 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(filtered.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

# Create relative abundance table
# To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# genus table: pruned.genus.table
# genus relative table: pruned.genus.rel.table
# variable: group

# Convert Phyloseq Object to DESEq object, correcting for any potential confounders
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ hiv + pairs + group)

# Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$group <- droplevels(diagdds$group)

# Choose the 'control' in your comparison variable 
diagdds$group <- relevel(diagdds$group, ref ="CCCs")

# Run the differential Analysis
diagdds<- DESeq(diagdds)

# Output the results from DESeq into a table
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant genes
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

## Tables

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$OTU,rownames(res))
# Replace spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

# Convert resuts table into a data.frame
res <- as.data.frame(res)

#Set names of results table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))

# Make the full trail the first column
res$names <- res$Gene.symbol

res$Gene.symbol <- res$row2

# Get abundance data - use sputum.otu.relative.table - THIS CODE ADDS REL ABUNDANCE AS THE DOT SIZE

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
cases.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("1"))
ccc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("2"))

# From relative table we should get the mean across the row of the otu table
cases.pruned.genus.rel.table.df <- data.frame(otu_table(cases.pruned.genus.rel.table))
cases.pruned.genus.rel.table.df.meanRA <- rowMeans(cases.pruned.genus.rel.table.df)

ccc.pruned.genus.rel.table.df <- data.frame(otu_table(ccc.pruned.genus.rel.table))
ccc.pruned.genus.rel.table.df.meanRA <- rowMeans(ccc.pruned.genus.rel.table.df)

# Subset AND reorder just the OTUs that we have 
cases.pruned.genus.rel.table.df.meanRA.save <- cases.pruned.genus.rel.table.df.meanRA[otu.to.save]
ccc.pruned.genus.rel.table.df.meanRA.save <- ccc.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res dataframe
res$abundance.cases <- cases.pruned.genus.rel.table.df.meanRA.save
res$abundance.ccc <- ccc.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need for pathway analysis
res.IPA <- res[,c("Gene.symbol", "logFC", "pvalue", "adj.P.Val")]
res.1 <- res[,c("Gene.symbol", "abundance.cases", "abundance.ccc", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Write tables to TXT file
write.table(res.IPA,file="taxa.sputum.cases.vs.CCCs.hiv.IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
write.table(res.1,file="taxa.sputum.cases.vs.CCCs.hiv.abundance.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Volcano plot

# Compute FDR in a log scales
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "red"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "darkgreen"

# Create a Variable for the size of the dots in the Volcano Plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="sputum.cases.vs.CCCs.hiv.volcano.fdr.0.2.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 1000 * res$abundance.cases, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 1000 * res$abundance.ccc,2)),alpha=0.7) + #Chose colors and size for dots
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=3,force=25, segment.colour="grey",segment.alpha=0.2) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)") + #label Y Axis
  theme #Set Theme
dev.off()



## ADJUSTING FOR HIV: SYMPTOMATIC CONTROLS (SCs) VS THEIR CLOSE CONTACTS (CCCs)

# Subset OTU table to induced sputum
subset.otu.table = subset_samples(otu.table, sample_c %in% c(2))

# Subset OTU table to SCs and CCSCs
subset.otu.table = subset_samples(subset.otu.table, group_c %in% c(0,3))

# Subset OTU table to genus level
subset.genus.table = tax_glom(subset.otu.table, taxrank = "Genus")

# Prune data to less than 100 genera remaining
filtered.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.25 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(filtered.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

# Create relative abundance table
# To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# genus table: pruned.genus.table
# genus relative table: pruned.genus.rel.table
# variable: group

# Convert Phyloseq Object to DESEq object, correcting for any potential confounders
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ hiv + pairs + group)

# Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$group <- droplevels(diagdds$group)

# Choose the 'control' in your comparison variable 
diagdds$group <- relevel(diagdds$group, ref ="CCSCs")

# Run the differential Analysis
diagdds<- DESeq(diagdds)

# Output the results from DESeq into a table
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant genes
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

## TABLES

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$OTU,rownames(res))
# Replace spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

# Convert resuts table into a data.frame
res <- as.data.frame(res)

#Set names of results table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))

# Make the full trail the first column
res$names <- res$Gene.symbol

res$Gene.symbol <- res$row2

# Get abundance data - use sputum.otu.relative.table - THIS CODE ADDS REL ABUNDANCE AS THE DOT SIZE

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
sc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("0"))
ccsc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("3"))

# From relative table we should get the mean across the row of the otu table
sc.pruned.genus.rel.table.df <- data.frame(otu_table(sc.pruned.genus.rel.table))
sc.pruned.genus.rel.table.df.meanRA <- rowMeans(sc.pruned.genus.rel.table.df)

ccsc.pruned.genus.rel.table.df <- data.frame(otu_table(ccsc.pruned.genus.rel.table))
ccsc.pruned.genus.rel.table.df.meanRA <- rowMeans(ccsc.pruned.genus.rel.table.df)

# Subset AND reorder just the OTUs that we have 
sc.pruned.genus.rel.table.df.meanRA.save <- sc.pruned.genus.rel.table.df.meanRA[otu.to.save]
ccsc.pruned.genus.rel.table.df.meanRA.save <- ccsc.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res dataframe
res$abundance.sc <- sc.pruned.genus.rel.table.df.meanRA.save
res$abundance.ccsc <- ccsc.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need for pathway analysis
res.IPA <- res[,c("Gene.symbol", "logFC", "pvalue", "adj.P.Val")]
res.1 <- res[,c("Gene.symbol", "abundance.sc", "abundance.ccsc", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Write tables to TXT file
write.table(res.IPA,file="taxa.sputum.SCs.vs.CCSCs.hiv.IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
write.table(res.1,file="taxa.sputum.SCs.vs.CCSCs.hiv.abundance.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Volcano plot

# Compute FDR in a log scales
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "red"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "darkgreen"

# Create a Variable for the size of the dots in the Volcano Plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="sputum.SCs.vs.CCSCs.hiv.volcano.fdr.0.2.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 500 * res$abundance.sc, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 500 * res$abundance.ccsc,2)),alpha=0.7) + #Chose colors and size for dots
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=3,force=25, segment.colour="grey",segment.alpha=0.2) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)") + #label Y Axis
  theme #Set Theme
dev.off()


## CLOSE CONTACTS OF CASES (CCCs) VS. CLOSE CONTACTS OF SYMPTOMATIC CONTROLS (CCSCs)

# Subset OTU table to induced sputum
subset.otu.table = subset_samples(otu.table, sample_c %in% c(2))

# Subset OTU table to CCCs and CCSCs
subset.otu.table = subset_samples(subset.otu.table, group_c %in% c(2,3))

# Subset OTU table to genus level
subset.genus.table = tax_glom(subset.otu.table, taxrank = "Genus")

# Prune data to less than 100 genera remaining
filtered.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.25 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(filtered.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

# Create relative abundance table
# To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# genus table: pruned.genus.table
# genus relative table: pruned.genus.rel.table
# variable: group

# Convert Phyloseq Object to DESEq object, correcting for any potential confounders
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ group)

# Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$group <- droplevels(diagdds$group)

# Choose the 'control' in your comparison variable 
diagdds$group <- relevel(diagdds$group, ref ="CCSCs")

# Run the differential Analysis
diagdds<- DESeq(diagdds)

# Output the results from DESeq into a table
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant genes
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

## TABLES

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$OTU,rownames(res))
# Replace spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

# Convert resuts table into a data.frame
res <- as.data.frame(res)

#Set names of results table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))

# Make the full trail the first column
res$names <- res$Gene.symbol

res$Gene.symbol <- res$row2

# Get abundance data - use sputum.otu.relative.table - THIS CODE ADDS REL ABUNDANCE AS THE DOT SIZE

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
ccc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("2"))
ccsc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("3"))

# From relative table we should get the mean across the row of the otu table
ccc.pruned.genus.rel.table.df <- data.frame(otu_table(ccc.pruned.genus.rel.table))
ccc.pruned.genus.rel.table.df.meanRA <- rowMeans(ccc.pruned.genus.rel.table.df)

ccsc.pruned.genus.rel.table.df <- data.frame(otu_table(ccsc.pruned.genus.rel.table))
ccsc.pruned.genus.rel.table.df.meanRA <- rowMeans(ccsc.pruned.genus.rel.table.df)

# Subset AND reorder just the OTUs that we have 
ccc.pruned.genus.rel.table.df.meanRA.save <- ccc.pruned.genus.rel.table.df.meanRA[otu.to.save]
ccsc.pruned.genus.rel.table.df.meanRA.save <- ccsc.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res dataframe
res$abundance.ccc <- ccc.pruned.genus.rel.table.df.meanRA.save
res$abundance.ccsc <- ccsc.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need for pathway analysis
res.IPA <- res[,c("Gene.symbol", "logFC", "pvalue", "adj.P.Val")]
res.1 <- res[,c("Gene.symbol", "abundance.ccc", "abundance.ccsc", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Write tables to TXT file
write.table(res.IPA,file="taxa.sputum.CCCs.vs.CCSCs.IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
write.table(res.1,file="taxa.sputum.CCCs.vs.CCSCs.abundance.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Volcano plot

# Compute FDR in a log scales
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "red"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "darkgreen"

# Create a Variable for the size of the dots in the Volcano Plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="sputum.CCCs.vs.CCSCs.volcano.fdr.0.2.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 1000 * res$abundance.ccc, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 1000 * res$abundance.ccsc,2)),alpha=0.7) + #Chose colors and size for dots
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=3,force=25, segment.colour="grey",segment.alpha=0.2) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)") + #label Y Axis
  theme #Set Theme
dev.off()

#############################################################################################
#############################################################################################
#############################################################################################

### STOOL

# Subset to stool
stool.otu.relative.table = subset_samples(otu.relative.table, sample_c %in% c(3))

## BETA DIVERSITY

# Create Distance Matrix
vegdist = distance(stool.otu.relative.table, "wunifrac")

# Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

# Calculate Sample variance for each PC
vars <- apply(CmdScale, 2, var)

# Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

# Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(stool.otu.relative.table), by = "row.names", all.x = TRUE)

# Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ group,data= newResults, mean)

# Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="group",suffixes=c("",".centroid"))

pdf("stool.wuniF.pdf", height = 10, width = 15)
ggplot(newResults, aes(PC1, PC2, color= group)) +
  geom_point(size=3,alpha=0.7) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values=c("red", "purple", "cyan", "forestgreen", "grey")) + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= group), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= group, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=group, size=10)) +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Subset cases and SCs
patients.stool.otu.relative.table = subset_samples(stool.otu.relative.table, group_c %in% c(0,1))
# Create distance object
patients = distance(patients.stool.otu.relative.table, "wunifrac")
# Statistics
adonis(patients ~ sample_data(patients.stool.otu.relative.table)$group)

# Statistics for other variables
adonis(patients ~ sample_data(patients.stool.otu.relative.table)$age)
adonis(patients ~ sample_data(patients.stool.otu.relative.table)$sex)
adonis(patients ~ sample_data(patients.stool.otu.relative.table)$ethnicity)
adonis(patients ~ sample_data(patients.stool.otu.relative.table)$previous_tb)
adonis(patients ~ sample_data(patients.stool.otu.relative.table)$smokes_cigarettes)
adonis(patients ~ sample_data(patients.stool.otu.relative.table)$hiv)
adonis(patients ~ sample_data(patients.stool.otu.relative.table)$tb_score_category)
adonis(patients ~ sample_data(patients.stool.otu.relative.table)$consumes_alcohol)

# Subset cases and CCCs
cases.ccc.stool.otu.relative.table = subset_samples(stool.otu.relative.table, group_c %in% c(1,2))
# Create distance object
cases = distance(cases.ccc.stool.otu.relative.table, "wunifrac")
# Statistics
adonis(cases ~ sample_data(cases.ccc.stool.otu.relative.table)$group)

# Subset SCs and CCSCs
sc.ccs.stool.otu.relative.table = subset_samples(stool.otu.relative.table, group_c %in% c(0,3))
# Create distance object
sc = distance(sc.ccs.stool.otu.relative.table, "wunifrac")
# Statistics
adonis(sc ~ sample_data(sc.ccs.stool.otu.relative.table)$group)

#########################################################
#########################################################
#########################################################

### BETA DIVERSITY IN EACH GROUP BY HIV STATUS

## Cases
# Subset cases
cases.stool.otu.relative.table = subset_samples(stool.otu.relative.table, group_c %in% c(1))

# Create Distance Matrix
vegdist = distance(cases.stool.otu.relative.table, "wunifrac")

# Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

# Calculate Sample variance for each PC
vars <- apply(CmdScale, 2, var)

# Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

# Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(cases.stool.otu.relative.table), by = "row.names", all.x = TRUE)

# Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ hiv,data= newResults, mean)

# Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="hiv",suffixes=c("",".centroid"))

pdf("stool.wunif.hiv.cases.pdf", height = 10, width = 10)
ggplot(newResults, aes(PC1, PC2, color= hiv)) + # Graph PC1 and PC2
  geom_point(size=3,alpha=0.7) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  scale_color_manual(values=c("orange","red", "grey")) + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= hiv), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= hiv, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=hiv, size=10)) +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Statistics
adonis(vegdist ~ sample_data(cases.stool.otu.relative.table)$hiv)


## Symptomatic controls (SCs)
# Subset SCs
sc.stool.otu.relative.table = subset_samples(stool.otu.relative.table, group_c %in% c(0))

# Create Distance Matrix
vegdist = distance(sc.stool.otu.relative.table, "wunifrac")

# Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

# Calculate Sample variance for each PC
vars <- apply(CmdScale, 2, var)

# Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

# Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(sc.stool.otu.relative.table), by = "row.names", all.x = TRUE)

# Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ hiv,data= newResults, mean)

# Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="hiv",suffixes=c("",".centroid"))

# Plot
pdf("stool.wunif.hiv.SCs.pdf", height = 10, width = 10)
ggplot(newResults, aes(PC1, PC2, color= hiv)) + # Graph PC1 and PC2
  geom_point(size=3,alpha=0.7) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  scale_color_manual(values=c("orange","red", "grey")) + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= hiv), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= hiv, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=hiv, size=10)) +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Statistics
adonis(vegdist ~ sample_data(sc.stool.otu.relative.table)$hiv)

## Close contacts of cases (CCCs)
# Subset CCCs
ccc.stool.otu.relative.table = subset_samples(stool.otu.relative.table, group_c %in% c(2))

#subset CCCs with only known hiv status
ccc.stool.otu.relative.table = subset_samples(ccc.stool.otu.relative.table, known_hiv_status %in% c(1))

# Create Distance Matrix
vegdist = distance(ccc.stool.otu.relative.table, "wunifrac")

# Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

# Calculate Sample variance for each PC
vars <- apply(CmdScale, 2, var)

# Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

# Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(ccc.stool.otu.relative.table), by = "row.names", all.x = TRUE)

# Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ hiv,data= newResults, mean)

# Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="hiv",suffixes=c("",".centroid"))

# Plot
pdf("stool.wunif.hiv.CCCs.pdf", height = 10, width = 10)
ggplot(newResults, aes(PC1, PC2, color= hiv)) + # Graph PC1 and PC2
  geom_point(size=3,alpha=0.7) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  scale_color_manual(values=c("orange","red", "grey")) + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= hiv), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= hiv, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=hiv, size=10)) +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Statistics
adonis(vegdist ~ sample_data(ccc.stool.otu.relative.table)$hiv)

## Close contacts of symptomatic controls (CCSCs)
# Subset CCSCs
ccsc.stool.otu.relative.table = subset_samples(stool.otu.relative.table, group_c %in% c(3))

#subset CCSCs with only known hiv status
ccsc.stool.otu.relative.table = subset_samples(ccsc.stool.otu.relative.table, known_hiv_status %in% c(1))

# Create Distance Matrix
vegdist = distance(ccsc.stool.otu.relative.table, "wunifrac")

# Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

# Calculate Sample variance for each PC
vars <- apply(CmdScale, 2, var)

# Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

# Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(ccsc.stool.otu.relative.table), by = "row.names", all.x = TRUE)

# Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ hiv,data= newResults, mean)

# Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="hiv",suffixes=c("",".centroid"))

# Plot
pdf("stool.wunif.hiv.CCSCs.pdf", height = 10, width = 10)
ggplot(newResults, aes(PC1, PC2, color= hiv)) + # Graph PC1 and PC2
  geom_point(size=3,alpha=0.7) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  scale_color_manual(values=c("orange","red", "grey")) + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= hiv), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= hiv, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=hiv, size=10)) +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

#Statistics not done because HIV-positive group only had one patient

#########################################################
#########################################################
#########################################################

### BETA DIVERSITY IN CLOSE CONTACTS OF CASES (CCCs) VS. CLOSE CONTACTS OF SYMPTOMATIC CONTROLS (CCSCs)

# Subset CCCs and CCSCs
contacts.stool.otu.relative.table = subset_samples(stool.otu.relative.table, group_c %in% c(2,3))
rownames(sample_data(contacts.stool.otu.relative.table))

# Create Distance Matrix
vegdist = distance(contacts.stool.otu.relative.table, "wunifrac")

# Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

# Calculate Sample variance for each PC
vars <- apply(CmdScale, 2, var)

# Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#  PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(contacts.stool.otu.relative.table), by = "row.names", all.x = TRUE)

# Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ group,data= newResults, mean)

# Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="group",suffixes=c("",".centroid"))

# Plot
pdf("stool.wuniF.contacts.pdf", height = 10, width = 15)
ggplot(newResults, aes(PC1, PC2, color= group)) + # Graph PC1 and PC2
  geom_point(size=3,alpha=0.7) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  scale_color_manual(values=c("purple", "cyan")) + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= group), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= group, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=group, size=10)) + 
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Statistics
adonis(vegdist ~ sample_data(contacts.stool.otu.relative.table)$group)

#########################################################
#########################################################
#########################################################

### DIFFERENTIAL ABUNDANCE ANALYSIS WITH DESEQ2

# Set Theme For Figures
theme <-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),
              axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), 
              legend.position="none")

# Load file
#load(file="all.samples.RData")

# Set alpha
alpha <- 0.20

## CASES VS SYMPTOMATIC CONTROLS (SCs)

# Subset OTU table to stool
subset.otu.table = subset_samples(otu.table, sample_c %in% c(3))

# Subset OTU table to cases and SCs
subset.otu.table = subset_samples(subset.otu.table, group_c %in% c(0,1))

# Subset OTU table to genus level
subset.genus.table = tax_glom(subset.otu.table, taxrank = "Genus")

# Prune data to less than 100 genera remaining
filtered.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.25 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(filtered.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

# Create relative abundance table
# To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# genus table: pruned.genus.table
# genus relative table: pruned.genus.rel.table
# variable: group

# Convert Phyloseq Object to DESEq object, correcting for any potential confounders
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ group)

# Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$group <- droplevels(diagdds$group)

# Choose the 'control' in your comparison variable --> positive is upregulated in cases, negative is down-regulated
diagdds$group <- relevel(diagdds$group, ref ="SCs")

# Run the differential Analysis
diagdds<- DESeq(diagdds)

# Output the results from DESeq into a table
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant genes
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

# Tables

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$OTU,rownames(res))
# Replace spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

# Convert resuts table into a data.frame
res <- as.data.frame(res)

#Set names of results table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))

# Make the full trail the first column
res$names <- res$Gene.symbol

res$Gene.symbol <- res$row2

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
cases.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("1"))
sc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("0"))

# From relative table we should get the mean across the row of the otu table
cases.pruned.genus.rel.table.df <- data.frame(otu_table(cases.pruned.genus.rel.table))
cases.pruned.genus.rel.table.df.meanRA <- rowMeans(cases.pruned.genus.rel.table.df)

sc.pruned.genus.rel.table.df <- data.frame(otu_table(sc.pruned.genus.rel.table))
sc.pruned.genus.rel.table.df.meanRA <- rowMeans(sc.pruned.genus.rel.table.df)

# Subset AND reorder just the OTUs that we have 
cases.pruned.genus.rel.table.df.meanRA.save <- cases.pruned.genus.rel.table.df.meanRA[otu.to.save]
sc.pruned.genus.rel.table.df.meanRA.save <- sc.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res dataframe
res$abundance.cases <- cases.pruned.genus.rel.table.df.meanRA.save
res$abundance.sc <- sc.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need for pathway analysis
res.IPA <- res[,c("Gene.symbol", "logFC", "pvalue", "adj.P.Val")]
res.1 <- res[,c("Gene.symbol", "abundance.cases", "abundance.sc", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Write tables to TXT file
write.table(res.IPA,file="taxa.stool.patients.IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
write.table(res.1,file="taxa.stool.patients.abundance.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Compute FDR in a log scales
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "red"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "darkgreen"

# Create a Variable for the size of the dots in the Volcano Plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="stool.patients.volcano.fdr.0.2.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 1000 * res$abundance.cases, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 1000 * res$abundance.sc,2)),alpha=0.7) + #Chose colors and size for dots
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=3,force=25, segment.colour="grey",segment.alpha=0.2) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)") + #label Y Axis
  theme #Set Theme
dev.off()


## CASES VS THEIR CLOSE CONTACTS (CCCs)

# Subset OTU table to stool
subset.otu.table = subset_samples(otu.table, sample_c %in% c(3))

# Subset OTU table to cases and CCCs
subset.otu.table = subset_samples(subset.otu.table, group_c %in% c(1,2))

# Subset OTU table to genus level
subset.genus.table = tax_glom(subset.otu.table, taxrank = "Genus")

# Prune data to less than 100 genera remaining
filtered.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.20 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(filtered.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

# Create relative abundance table
# To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# genus table: pruned.genus.table
# genus relative table: pruned.genus.rel.table
# variable: group

# Convert Phyloseq Object to DESEq object, correcting for any potential confounders
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ pairs + group)

# Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$group <- droplevels(diagdds$group)

# Choose the 'control' in your comparison variable 
diagdds$group <- relevel(diagdds$group, ref ="CCCs")

# Run the differential Analysis
diagdds<- DESeq(diagdds)

# Output the results from DESeq into a table
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant genes
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

## Tables

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$OTU,rownames(res))
# Replace spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

# Convert resuts table into a data.frame
res <- as.data.frame(res)

#Set names of results table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))

# Make the full trail the first column
res$names <- res$Gene.symbol

res$Gene.symbol <- res$row2

# Get abundance data - use stool.otu.relative.table - THIS CODE ADDS REL ABUNDANCE AS THE DOT SIZE

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
cases.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("1"))
ccc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("2"))

# From relative table we should get the mean across the row of the otu table
cases.pruned.genus.rel.table.df <- data.frame(otu_table(cases.pruned.genus.rel.table))
cases.pruned.genus.rel.table.df.meanRA <- rowMeans(cases.pruned.genus.rel.table.df)

ccc.pruned.genus.rel.table.df <- data.frame(otu_table(ccc.pruned.genus.rel.table))
ccc.pruned.genus.rel.table.df.meanRA <- rowMeans(ccc.pruned.genus.rel.table.df)

# Subset AND reorder just the OTUs that we have 
cases.pruned.genus.rel.table.df.meanRA.save <- cases.pruned.genus.rel.table.df.meanRA[otu.to.save]
ccc.pruned.genus.rel.table.df.meanRA.save <- ccc.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res dataframe
res$abundance.cases <- cases.pruned.genus.rel.table.df.meanRA.save
res$abundance.ccc <- ccc.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need for pathway analysis
res.IPA <- res[,c("Gene.symbol", "logFC", "pvalue", "adj.P.Val")]
res.1 <- res[,c("Gene.symbol", "abundance.cases", "abundance.ccc", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Write tables to TXT file
write.table(res.IPA,file="taxa.stool.cases.vs.CCCs.IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
write.table(res.1,file="taxa.stool.cases.vs.CCCs.abundance.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Volcano plot

# Compute FDR in a log scales
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "red"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "darkgreen"

# Create a Variable for the size of the dots in the Volcano Plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="stool.cases.vs.CCCs.volcano.fdr.0.2.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 1000 * res$abundance.cases, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 1000 * res$abundance.ccc,2)),alpha=0.7) + #Chose colors and size for dots
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=3,force=25, segment.colour="grey",segment.alpha=0.2) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)") + #label Y Axis
  theme #Set Theme
dev.off()

## SYMPTOMATIC CONTROLS (SCs) VS THEIR CLOSE CONTACTS (CCCs)

# Subset OTU table to stool
subset.otu.table = subset_samples(otu.table, sample_c %in% c(3))

# Subset OTU table to SCs and CCSCs
subset.otu.table = subset_samples(subset.otu.table, group_c %in% c(0,3))

# Subset OTU table to genus level
subset.genus.table = tax_glom(subset.otu.table, taxrank = "Genus")

# Prune data to less than 100 genera remaining
filtered.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.20 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(filtered.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

# Create relative abundance table
# To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# genus table: pruned.genus.table
# genus relative table: pruned.genus.rel.table
# variable: group

# Convert Phyloseq Object to DESEq object, correcting for any potential confounders
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ pairs + group)

# Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$group <- droplevels(diagdds$group)

# Choose the 'control' in your comparison variable 
diagdds$group <- relevel(diagdds$group, ref ="CCSCs")

# Run the differential Analysis
diagdds<- DESeq(diagdds)

# Output the results from DESeq into a table
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant genes
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

## TABLES

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$OTU,rownames(res))
# Replace spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

# Convert resuts table into a data.frame
res <- as.data.frame(res)

#Set names of results table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))

# Make the full trail the first column
res$names <- res$Gene.symbol

res$Gene.symbol <- res$row2

# Get abundance data - use stool.otu.relative.table - THIS CODE ADDS REL ABUNDANCE AS THE DOT SIZE

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
sc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("0"))
ccsc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("3"))

# From relative table we should get the mean across the row of the otu table
sc.pruned.genus.rel.table.df <- data.frame(otu_table(sc.pruned.genus.rel.table))
sc.pruned.genus.rel.table.df.meanRA <- rowMeans(sc.pruned.genus.rel.table.df)

ccsc.pruned.genus.rel.table.df <- data.frame(otu_table(ccsc.pruned.genus.rel.table))
ccsc.pruned.genus.rel.table.df.meanRA <- rowMeans(ccsc.pruned.genus.rel.table.df)

# Subset AND reorder just the OTUs that we have 
sc.pruned.genus.rel.table.df.meanRA.save <- sc.pruned.genus.rel.table.df.meanRA[otu.to.save]
ccsc.pruned.genus.rel.table.df.meanRA.save <- ccsc.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res dataframe
res$abundance.sc <- sc.pruned.genus.rel.table.df.meanRA.save
res$abundance.ccsc <- ccsc.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need for pathway analysis
res.IPA <- res[,c("Gene.symbol", "logFC", "pvalue", "adj.P.Val")]
res.1 <- res[,c("Gene.symbol", "abundance.sc", "abundance.ccsc", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Write tables to TXT file
write.table(res.IPA,file="taxa.stool.SCs.vs.CCSCs.IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
write.table(res.1,file="taxa.stool.SCs.vs.CCSCs.abundance.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

## Volcano plot

# Compute FDR in a log scales
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "red"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "darkgreen"

# Create a Variable for the size of the dots in the Volcano Plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="stool.SCs.vs.CCSCs.volcano.fdr.0.2.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 500 * res$abundance.sc, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 500 * res$abundance.ccsc,2)),alpha=0.7) + #Chose colors and size for dots
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=3,force=25, segment.colour="grey",segment.alpha=0.2) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)") + #label Y Axis
  theme #Set Theme
dev.off()

## ADJUSTING FOR HIV: CASES VS SYMPTOMATIC CONTROLS

# Subset OTU table to stool
subset.otu.table = subset_samples(otu.table, sample_c %in% c(3))

# Subset OTU table to cases and SCs
subset.otu.table = subset_samples(subset.otu.table, group_c %in% c(0,1))

# Subset OTU table to genus level
subset.genus.table = tax_glom(subset.otu.table, taxrank = "Genus")

# Prune data to less than 100 genera remaining
filtered.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.25 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(filtered.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

# Create relative abundance table
# To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# genus table: pruned.genus.table
# genus relative table: pruned.genus.rel.table
# variable: group

# Convert Phyloseq Object to DESEq object, correcting for any potential confounders
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ hiv + group)

# Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$group <- droplevels(diagdds$group)

# Choose the 'control' in your comparison variable --> positive is upregulated in cases, negative is down-regulated
diagdds$group <- relevel(diagdds$group, ref ="SCs")

# Run the differential Analysis
diagdds<- DESeq(diagdds)

# Output the results from DESeq into a table
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant genes
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

## Tables

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$OTU,rownames(res))
# Replace spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

# Convert resuts table into a data.frame
res <- as.data.frame(res)

#Set names of results table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))

# Make the full trail the first column
res$names <- res$Gene.symbol

res$Gene.symbol <- res$row2

# Get abundance data - use stool.otu.relative.table - THIS CODE ADDS REL ABUNDANCE AS THE DOT SIZE

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
cases.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("1"))
sc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("0"))

# From relative table we should get the mean across the row of the otu table
cases.pruned.genus.rel.table.df <- data.frame(otu_table(cases.pruned.genus.rel.table))
cases.pruned.genus.rel.table.df.meanRA <- rowMeans(cases.pruned.genus.rel.table.df)

sc.pruned.genus.rel.table.df <- data.frame(otu_table(sc.pruned.genus.rel.table))
sc.pruned.genus.rel.table.df.meanRA <- rowMeans(sc.pruned.genus.rel.table.df)

# Subset AND reorder just the OTUs that we have 
cases.pruned.genus.rel.table.df.meanRA.save <- cases.pruned.genus.rel.table.df.meanRA[otu.to.save]
sc.pruned.genus.rel.table.df.meanRA.save <- sc.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res dataframe
res$abundance.cases <- cases.pruned.genus.rel.table.df.meanRA.save
res$abundance.sc <- sc.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need for pathway analysis
res.IPA <- res[,c("Gene.symbol", "logFC", "pvalue", "adj.P.Val")]
res.1 <- res[,c("Gene.symbol", "abundance.cases", "abundance.sc", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Write tables to TXT file
write.table(res.IPA,file="taxa.stool.patients.hiv.IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
write.table(res.1,file="taxa.stool.patients.hiv.abundance.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Volcano plot

# Compute FDR in a log scales
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "red"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "darkgreen"

# Create a Variable for the size of the dots in the Volcano Plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="stool.patients.hiv.volcano.fdr.0.2.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 1000 * res$abundance.cases, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 1000 * res$abundance.sc,2)),alpha=0.7) + #Chose colors and size for dots
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=3,force=25, segment.colour="grey",segment.alpha=0.2) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)") + #label Y Axis
  theme #Set Theme
dev.off()

## ADJUSTING FOR HIV: CASES VS THEIR CLOSE CONTACTS (CCCs)

# Subset OTU table to stool
subset.otu.table = subset_samples(otu.table, sample_c %in% c(3))

# Subset OTU table to cases and CCCs
subset.otu.table = subset_samples(subset.otu.table, group_c %in% c(1,2))

# Subset OTU table to genus level
subset.genus.table = tax_glom(subset.otu.table, taxrank = "Genus")

# Prune data to less than 100 genera remaining
filtered.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.20 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(filtered.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

# Create relative abundance table
# To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# genus table: pruned.genus.table
# genus relative table: pruned.genus.rel.table
# variable: group

# Convert Phyloseq Object to DESEq object, correcting for any potential confounders
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ hiv + pairs + group)

# Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$group <- droplevels(diagdds$group)

# Choose the 'control' in your comparison variable 
diagdds$group <- relevel(diagdds$group, ref ="CCCs")

# Run the differential Analysis
diagdds<- DESeq(diagdds)

# Output the results from DESeq into a table
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant genes
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

## Tables

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$OTU,rownames(res))
# Replace spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

# Convert resuts table into a data.frame
res <- as.data.frame(res)

#Set names of results table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))

# Make the full trail the first column
res$names <- res$Gene.symbol

res$Gene.symbol <- res$row2

# Get abundance data - use stool.otu.relative.table - THIS CODE ADDS REL ABUNDANCE AS THE DOT SIZE

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
cases.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("1"))
ccc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("2"))

# From relative table we should get the mean across the row of the otu table
cases.pruned.genus.rel.table.df <- data.frame(otu_table(cases.pruned.genus.rel.table))
cases.pruned.genus.rel.table.df.meanRA <- rowMeans(cases.pruned.genus.rel.table.df)

ccc.pruned.genus.rel.table.df <- data.frame(otu_table(ccc.pruned.genus.rel.table))
ccc.pruned.genus.rel.table.df.meanRA <- rowMeans(ccc.pruned.genus.rel.table.df)

# Subset AND reorder just the OTUs that we have 
cases.pruned.genus.rel.table.df.meanRA.save <- cases.pruned.genus.rel.table.df.meanRA[otu.to.save]
ccc.pruned.genus.rel.table.df.meanRA.save <- ccc.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res dataframe
res$abundance.cases <- cases.pruned.genus.rel.table.df.meanRA.save
res$abundance.ccc <- ccc.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need for pathway analysis
res.IPA <- res[,c("Gene.symbol", "logFC", "pvalue", "adj.P.Val")]
res.1 <- res[,c("Gene.symbol", "abundance.cases", "abundance.ccc", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Write tables to TXT file
write.table(res.IPA,file="taxa.stool.cases.vs.CCCs.hiv.IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
write.table(res.1,file="taxa.stool.cases.vs.CCCs.hiv.abundance.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Volcano plot

# Compute FDR in a log scales
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "red"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "darkgreen"

# Create a Variable for the size of the dots in the Volcano Plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="stool.cases.vs.CCCs.hiv.volcano.fdr.0.2.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 1000 * res$abundance.cases, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 1000 * res$abundance.ccc,2)),alpha=0.7) + #Chose colors and size for dots
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=3,force=25, segment.colour="grey",segment.alpha=0.2) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)") + #label Y Axis
  theme #Set Theme
dev.off()

## ADJUSTING FOR HIV: SYMPTOMATIC CONTROLS (SCs) VS THEIR CLOSE CONTACTS (CCCs)

# Subset OTU table to stool
subset.otu.table = subset_samples(otu.table, sample_c %in% c(3))

# Subset OTU table to SCs and CCSCs
subset.otu.table = subset_samples(subset.otu.table, group_c %in% c(0,3))

# Subset OTU table to genus level
subset.genus.table = tax_glom(subset.otu.table, taxrank = "Genus")

# Prune data to less than 100 genera remaining
filtered.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.20 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(filtered.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

# Create relative abundance table
# To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# genus table: pruned.genus.table
# genus relative table: pruned.genus.rel.table
# variable: group

# Convert Phyloseq Object to DESEq object, correcting for any potential confounders
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ hiv + pairs + group)

# Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$group <- droplevels(diagdds$group)

# Choose the 'control' in your comparison variable 
diagdds$group <- relevel(diagdds$group, ref ="CCSCs")

# Run the differential Analysis
diagdds<- DESeq(diagdds)

# Output the results from DESeq into a table
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant genes
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

## TABLES

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$OTU,rownames(res))
# Replace spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

# Convert resuts table into a data.frame
res <- as.data.frame(res)

#Set names of results table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))

# Make the full trail the first column
res$names <- res$Gene.symbol

res$Gene.symbol <- res$row2

# Get abundance data - use stool.otu.relative.table - THIS CODE ADDS REL ABUNDANCE AS THE DOT SIZE

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
sc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("0"))
ccsc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("3"))

# From relative table we should get the mean across the row of the otu table
sc.pruned.genus.rel.table.df <- data.frame(otu_table(sc.pruned.genus.rel.table))
sc.pruned.genus.rel.table.df.meanRA <- rowMeans(sc.pruned.genus.rel.table.df)

ccsc.pruned.genus.rel.table.df <- data.frame(otu_table(ccsc.pruned.genus.rel.table))
ccsc.pruned.genus.rel.table.df.meanRA <- rowMeans(ccsc.pruned.genus.rel.table.df)

# Subset AND reorder just the OTUs that we have 
sc.pruned.genus.rel.table.df.meanRA.save <- sc.pruned.genus.rel.table.df.meanRA[otu.to.save]
ccsc.pruned.genus.rel.table.df.meanRA.save <- ccsc.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res dataframe
res$abundance.sc <- sc.pruned.genus.rel.table.df.meanRA.save
res$abundance.ccsc <- ccsc.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need for pathway analysis
res.IPA <- res[,c("Gene.symbol", "logFC", "pvalue", "adj.P.Val")]
res.1 <- res[,c("Gene.symbol", "abundance.sc", "abundance.ccsc", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Write tables to TXT file
write.table(res.IPA,file="taxa.stool.SCs.vs.CCSCs.hiv.IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
write.table(res.1,file="taxa.stool.SCs.vs.CCSCs.hiv.abundance.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Volcano plot

# Compute FDR in a log scales
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "red"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "darkgreen"

# Create a Variable for the size of the dots in the Volcano Plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="stool.SCs.vs.CCSCs.hiv.volcano.fdr.0.2.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 1000 * res$abundance.sc, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 1000 * res$abundance.ccsc,2)),alpha=0.7) + #Chose colors and size for dots
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=3,force=25, segment.colour="grey",segment.alpha=0.2) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)") + #label Y Axis
  theme #Set Theme
dev.off()

## CLOSE CONTACTS OF CASES (CCCs) VS. CLOSE CONTACTS OF SYMPTOMATIC CONTROLS (CCSCs)

# Subset OTU table to stool
subset.otu.table = subset_samples(otu.table, sample_c %in% c(3))

# Subset OTU table to CCCs and CCSCs
subset.otu.table = subset_samples(subset.otu.table, group_c %in% c(2,3))

# Subset OTU table to genus level
subset.genus.table = tax_glom(subset.otu.table, taxrank = "Genus")

# Prune data to less than 100 genera remaining
filtered.genus.table = genefilter_sample(subset.genus.table, filterfun_sample(function(x) x > 0.05), A = 0.15 * nsamples(subset.genus.table))
pruned.genus.table = prune_taxa(filtered.genus.table, subset.genus.table)
colnames(sample_data(pruned.genus.table))
rownames(sample_data(pruned.genus.table))
ntaxa(pruned.genus.table)

# Create relative abundance table
# To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# genus table: pruned.genus.table
# genus relative table: pruned.genus.rel.table
# variable: group

# Convert Phyloseq Object to DESEq object, correcting for any potential confounders
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ group)

# Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$group <- droplevels(diagdds$group)

# Choose the 'control' in your comparison variable 
diagdds$group <- relevel(diagdds$group, ref ="CCSCs")

# Run the differential Analysis
diagdds<- DESeq(diagdds)

# Output the results from DESeq into a table
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant genes
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

## TABLES

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$OTU,rownames(res))
# Replace spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

# Convert resuts table into a data.frame
res <- as.data.frame(res)

#Set names of results table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))

# Make the full trail the first column
res$names <- res$Gene.symbol

res$Gene.symbol <- res$row2

# Get abundance data - use stool.otu.relative.table - THIS CODE ADDS REL ABUNDANCE AS THE DOT SIZE

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
ccc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("2"))
ccsc.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, group_c %in% c("3"))

# From relative table we should get the mean across the row of the otu table
ccc.pruned.genus.rel.table.df <- data.frame(otu_table(ccc.pruned.genus.rel.table))
ccc.pruned.genus.rel.table.df.meanRA <- rowMeans(ccc.pruned.genus.rel.table.df)

ccsc.pruned.genus.rel.table.df <- data.frame(otu_table(ccsc.pruned.genus.rel.table))
ccsc.pruned.genus.rel.table.df.meanRA <- rowMeans(ccsc.pruned.genus.rel.table.df)

# Subset AND reorder just the OTUs that we have 
ccc.pruned.genus.rel.table.df.meanRA.save <- ccc.pruned.genus.rel.table.df.meanRA[otu.to.save]
ccsc.pruned.genus.rel.table.df.meanRA.save <- ccsc.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res dataframe
res$abundance.ccc <- ccc.pruned.genus.rel.table.df.meanRA.save
res$abundance.ccsc <- ccsc.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need for pathway analysis
res.IPA <- res[,c("Gene.symbol", "logFC", "pvalue", "adj.P.Val")]
res.1 <- res[,c("Gene.symbol", "abundance.ccc", "abundance.ccsc", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Write tables to TXT file
write.table(res.IPA,file="taxa.stool.CCCs.vs.CCSCs.IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
write.table(res.1,file="taxa.stool.CCCs.vs.CCSCs.abundance.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Volcano plot

# Compute FDR in a log scales
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "red"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "darkgreen"

# Create a Variable for the size of the dots in the Volcano Plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="stool.CCCs.vs.CCSCs.volcano.fdr.0.2.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 2000 * res$abundance.ccc, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 2000 * res$abundance.ccsc,2)),alpha=0.7) + #Chose colors and size for dots
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=3,force=25, segment.colour="grey",segment.alpha=0.2) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)") + #label Y Axis
  theme #Set Theme
dev.off()
