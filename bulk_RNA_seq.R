---
title: "Cambi bulk RNA-seq"
author: "AnnaS"
date: "Last edited `r format (Sys.time(), '%d %B %Y')`"
output: 
  html_document:
    toc: true
    toc_float: true
    number_sections: true
    code_folding: hide
    toc_depth: 4
    fig_path: figure-html/
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r setseed}
set.seed(1)
```

Goal: To study the mRNA expression between the genotypes


# Experimental set-up and preparing the data

Samples(genotypes):  
- NS1-C2d  
- NS1-TGFb2d  
- NS1-IL1a2d  
- NS1-TplusI2d  
- NS1-C4d  
- NS1-TGFb4d  
- NS1-IL1a4d  
- NS1-TplusI4d  
- NS1-TGFbswap  
- NS1-TGFbwash  
- NS1-IL1aswap  
- NS1-IL1awash  
- NS2-C2d   
- NS2-TGFb2d  
- NS2-IL1a2d  
- NS2-TplusI2d  
- NS2-C4d  
- NS2-TGFb4d  
- NS2-IL1a4d  
- NS2-TplusI4d  
- NS2-TGFbswap  
- NS2-TGFbwash  
- NS2-IL1aswap  
- NS2-IL1awash  
- NS3-C2d  
- NS3-TGFb2d  
- NS3-IL1a2d  
- NS3-TplusI2d  
- NS3-C4d  
- NS3-TGFb4d  
- NS3-IL1a4d  
- NS3-TplusI4d  
- NS3-TGFbswap  
- NS3-TGFbwash  
- NS3-IL1aswap  
- NS3-IL1awash  



Load the libraries
```{r load libs, results='hide', message=FALSE}
library(readr)
library(dplyr)
library(tidyverse)
library(BiocManager)
library(DESeq2) 
library(vsn)
library(pheatmap)
library(ComplexHeatmap)
library(ggplot2)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(hexbin)
library(ggrepel)
library(clusterProfiler)
library(edgeR) #pca
library(factoextra) #pca
library(RColorBrewer) #heatmap
library(progeny)


experiment_id <- "Cambi_bulk_RNA-seq"
```

## Load the counts dataset from seq2science results directory
```{r loading counts, message=FALSE}
#read in counts and filter out low counts
data.counts <- read.table("../GRCm39-counts.tsv", header=TRUE, row.names=1)
data.counts[,c(1,2,3,4,13,14,15,16,25,26,27,28)] <- NULL
colnames(data.counts) <- gsub("4d", "", gsub("TplusI", "mix", gsub("TGFbswap", "TGFb->IL1a", gsub("IL1aswap", "IL1a->TGFb", colnames(data.counts)))))

#reorder the columns so the replicates are next to each other
data.counts.reord <- data.counts[,c(1,9,17,6,14,22,2,10,18,4,12,20,5,13,21,7,15,23,3,11,19,8,16,24)]

#total number of mapped reads in the pool
data.reads.tot <- sum(data.counts.reord)
data.reads.tot

#check the number of mapped reads per sample based on the gene count table
data.reads.sum <-colSums(data.counts.reord)
data.reads.sum
barplot(data.reads.sum, main = "Number of mapped reads per sample", col = c("seashell1", "seashell1", "seashell1", "salmon", "salmon", "salmon", "indianred3", "indianred3", "indianred3", "indianred4", "indianred4", "indianred4", "seashell3", "seashell3", "seashell3", "seashell4", "seashell4", "seashell4", "cyan2", "cyan2", "cyan2", "lightseagreen", "lightseagreen", "lightseagreen", "turquoise4", "turquoise4", "turquoise4"), las=2, cex.names = 0.7)

#check the number of detected genes per sample
data.gene.sum <- colSums(data.counts.reord !=0)
data.gene.sum
barplot(data.gene.sum, main = "Number of detected genes per sample", col = c("seashell1", "seashell1", "seashell1", "salmon", "salmon", "salmon", "indianred3", "indianred3", "indianred3", "indianred4", "indianred4", "indianred4", "seashell3", "seashell3", "seashell3", "seashell4", "seashell4", "seashell4", "cyan2", "cyan2", "cyan2", "lightseagreen", "lightseagreen", "lightseagreen", "turquoise4", "turquoise4", "turquoise4"), las=2, cex.names = 0.7)



#Violin plot and boxplot of the log data counts

#filter out rows containing only 0
data.wo.0 <- data.counts.reord[rowSums(data.counts.reord != 0) > 0,]
#log transform the data
data.log <- as.data.frame(log2(data.counts.reord))
#gather data in two columns
data.log.gather <- gather(data.log, Sample, Value)


#plot the gather data
##data.log
ggplot(data.log.gather, aes(x=Sample, y=Value)) + geom_violin() + ggtitle("Violin plot of log data.counts") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(data.log.gather, aes(x=Sample, y=Value)) + geom_boxplot(alpha=0.2) + ggtitle("Boxplot of log data.counts") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```

## Gene annotation
```{r gene ann,fig.keep='all', message=FALSE}
#filter the data, more than 10 reads inmore than 3 columns
data.counts <- data.counts.reord[rowSums(data.counts<10) <3,]

#rename the genes from ENSMUSG00000040952.17 to ENSMUSG00000040952
rownames(data.counts.reord) <- sub("\\..*", "", rownames(data.counts.reord))

#annotate genes in counts, if there is NA gene name, keep the gene symbol
gene_names <- mapIds(org.Mm.eg.db, keys = rownames(data.counts.reord), keytype = "ENSEMBL", column = "SYMBOL")
data.counts.reord$genes <- ifelse(is.na(gene_names), rownames(data.counts.reord), gene_names)

#change the gene column to rownames
data.merge.duplicates <- data.counts.reord %>% group_by(genes) %>% summarise_all(sum)
data.counts.final <- column_to_rownames(data.merge.duplicates, var="genes")

```


# DESeq

## Prepare colData for DESeq
```{r coldata, message=FALSE}
#Preparing the metadata table for DESeq2

metadata <- data.frame(sample = colnames(data.counts.final))

metadata$stimulation <- ifelse(grepl("TGFb->IL1a", metadata$sample), "TGFb->IL1a",
                              ifelse(grepl("IL1a->TGFb", metadata$sample), "IL1a->TGFb",
                                     ifelse(grepl("TGFbwash", metadata$sample), "TGFbwash",
                                            ifelse(grepl("IL1awash", metadata$sample), "IL1awash",
                                                   ifelse(grepl("C", metadata$sample), "ctrl",
                                                          ifelse(grepl("mix", metadata$sample), "mix",
                                                                 ifelse(grepl("TGFb", metadata$sample), "TGFb",
                                                                        ifelse(grepl("IL1a", metadata$sample), "IL1a",
                                           "other"))))))))

metadata$batch <- ifelse(grepl("NS1", metadata$sample), "1+2",
                         ifelse(grepl("NS2", metadata$sample), "1+2",
                                ifelse(grepl("NS3", metadata$sample), "3",
                                       "other")))

rownames(metadata) <- metadata$sample
metadata$sample <- NULL
metadata
```


## Sample distance plot without the outlier
Run the actual DESeq() function (estimates size factors, dispersion and fits the negative binomial GLM) and run the variance stabilizing transformation vst().

```{r sample dist, warning=FALSE, message=FALSE}
#make sure row names in colData match column names in counts and are in correct order
all(colnames(data.counts.final) %in% rownames(metadata))
all(colnames(data.counts.final) == rownames(metadata))

#construct a DEseq object
dds <- DESeqDataSetFromMatrix(countData = data.counts.final,
                       colData = metadata,
                       design = ~ batch + stimulation)

#remove rows with low count, keep rows with at least 10 reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#run DESEq
dds <- DESeq(dds)
res <- results(dds)
res

#export data
dds.as.df <- counts(dds, normalized=TRUE)
write.csv2(dds.as.df, "../output/4days/dds4days0.1.csv")

#data transformation
vsd <- vst(dds, blind = TRUE) # blind=TRUE to calculate across-all-samples variability, blind=FALSE to calculate within-group variability 
assay(vsd) <- limma::removeBatchEffect(assay(vsd), batch=vsd$batch, design=model.matrix(~stimulation, colData(vsd)))
plotPCA(vsd, "stimulation") + ggtitle("PCA on dds data")
vsd.assay <- assay(vsd)
write.csv2(vsd.assay, "../output/4days/vsd4days.csv")

#heatmap
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,"stimulation"])
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

#heatmap of sample distance
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(vsd))
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
```

### Heatmap of all DEG
Comment this block when knitting.
Heatmap based on the DEG in contrasts: ctrl vs TGFb+IL1a+mix, mix vs swaps+TGFb, ctrl vs mix+swaps+TGFb, ctrl vs TGFb+TGFbwash.
Total 373 genes.

```{r all DEG heatmap, warning=FALSE, message=FALSE}
# #merge all DEG
# genes.DE.all <- unique(c(genes.DE1, genes.DE2, genes.DE3, genes.DE4))
# 
# # Subset the vsd object to include only the genes of interest
# heatmap <- assay(vsd)[genes.DE.all,]
# 
# # Scale the counts for better visualization in the heatmap
# heatmap <- t(scale(t(heatmap)))
# 
# # Plot the heatmap
# pdf("../output/4days/4days-0.1-Heatmap-of-DEG.pdf")
# ComplexHeatmap::Heatmap(heatmap,
#                         cluster_rows = TRUE,
#                         cluster_columns = FALSE,
#                         column_names_gp = gpar(fontsize = 6),
#                         row_names_gp = gpar(fontsize = 1),
#                         col = colprgn,
#                         column_title = "Heatmap of all DEG")
# dev.off()
# 
# # almost 1GB file
# tiff("../output/4days/4days-0.1-Heatmap-of-DEG.tiff", width = 9, height = 12, units = "in", res = 1000)
# ComplexHeatmap::Heatmap(heatmap,
#                         cluster_rows = TRUE,
#                         cluster_columns = FALSE,
#                         column_names_gp = gpar(fontsize = 6),
#                         row_names_gp = gpar(fontsize = 2,5),
#                         col = colprgn,
#                         column_title = "Heatmap of all DEG")
# dev.off()

```



## PROGENy
Compute PROGENy scores for every sample (with the replicates) using the DESeq2 normalised counts. 100 most responsive genes per pathway are used.
Inside PROGENy, one can find gene signatures for 14 different pathways:  
- **Androgen:** involved in the growth and development of the male reproductive organs.  
- **EGFR:** regulates growth, survival, migration, apoptosis, proliferation, and differentiation in mammalian cells  
- **Estrogen:** promotes the growth and development of the female reproductive organs.  
- **Hypoxia:** promotes angiogenesis and metabolic reprogramming when O2 levels are low.  
- **JAK-STAT:** involved in immunity, cell division, cell death, and tumor formation.  
- **MAPK:** integrates external signals and promotes cell growth and proliferation.  
- **NFkB:** regulates immune response, cytokine production and cell survival.  
- **p53:** regulates cell cycle, apoptosis, DNA repair and tumor suppression.  
- **PI3K:** promotes growth and proliferation.  
- **TGFb:** involved in development, homeostasis, and repair of most tissues.  
- **TNFa:** mediates haematopoiesis, immune surveillance, tumour regression and protection from infection.  
- **Trail:** induces apoptosis.  
- **VEGF:** mediates angiogenesis, vascular permeability, and cell migration.  
- **WNT:** regulates organ morphogenesis during development and tissue repair.  
```{r progeny, warning=FALSE, message=FALSE}
dds.as.df <- counts(dds, normalized=TRUE)
colprgn = colorRampPalette(brewer.pal(9, "PRGn"))(255)

#vsd data
vsd.as.df <- as.matrix(assay(vsd))
progeny.scores <- progeny(vsd.as.df, scale = TRUE, organism = "Mouse", top = 100)
progeny.vector <- as.vector(progeny.scores)
paletteLength <- 100

progenyBreaks <- c(seq(min(progeny.vector), 0, 
    length.out=ceiling(paletteLength/2) + 1),
    seq(max(progeny.vector)/paletteLength, 
    max(progeny.vector), 
    length.out=floor(paletteLength/2)))

pheatmap.plot <- pheatmap(t(progeny.scores),fontsize=12, 
    fontsize_row = 10, fontsize_col = 10, 
    cluster_rows = TRUE, cluster_cols = FALSE, 
    color=colprgn, breaks = progenyBreaks, 
    main = "PROGENy (100) 4 days",
    treeheight_col = 0,  border_color = NA)
pheatmap.plot
tiff("../output/4days/progeny.tiff", width = 12, height = 7, units = "in", res = 300)
plot(pheatmap.plot)
dev.off()
```

## Different contrasts 
In the `results()` function, the comparison of interest can be defined using contrasts.  
To explore the output of `results()` function, I transformed the dds data, ordered results by significance (p-value), scaled the data, and plotted top 20 differently expressed genes in heatmap. The heatmaps compute the dispersion again, so they do not acknowledge the batch effect defined in the DESeq design as a covariate. Therefore, just for the heatmaps, I ran also function limma::removeBatchEffect()  
Volcano plot shows significantly (y-axis p adj. value) up/downregulated (x-axis log2FoldChange) genes.  
In GO enrichment analysis, DE genes below padj 0.05 and log2FoldChange over 1.5 were chosen. The last plot summarizes GO terms into clusters based on their semantic similarity of the functions.  


### ctrl vs mix, TGFb, TGFb->IL1a, IL1a->TGFb {.tabset .tabset-fade .tabset-pills}
```{r ctrl vs mix+TGFb+TGFbswap+IL1aswap, message=FALSE, warning=FALSE}
#reorder the columns and metadata for a given contrast
data.counts.sep <- data.counts.final[,c(1:3, 7:17)]
metadata.sep <- metadata[c(1:3, 7:17),]

#make sure row names in colData match column names in counts and are in correct order
all(colnames(data.counts.sep) %in% rownames(metadata.sep))
all(colnames(data.counts.sep) == rownames(metadata.sep))

#construct a DEseq object
dds <- DESeqDataSetFromMatrix(countData = data.counts.sep,
                       colData = metadata.sep,
                       design = ~ batch + stimulation)

#remove rows with low count, keep rows with at least 10 reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$stimulation <- relevel(dds$stimulation, ref = "ctrl")

#run DESeq
dds <- DESeq(dds)
res <- results(dds)
res

#explore results
summary(res)
colData(dds)

#data transformation
vsd <- vst(dds, blind = TRUE) # blind=TRUE to calculate across-all-samples variability (when we do not expect most genes to change between treatment groups), blind=FALSE to calculate within-group variability. Blind=TRUE calculates the dispersion again, so any batch effect can be visible again. It does not use the design to remove variation in the data -> it does not remove variation that can be associated with batch covariate.

#effect of data transformation on the variance
ntd <- normTransform(dds)

#remove batch effect with limma
assay(vsd) <- limma::removeBatchEffect(assay(vsd), batch=vsd$batch, design=model.matrix(~stimulation, colData(vsd)))

plotPCA(vsd, intgroup="stimulation") + geom_label_repel(aes(label=name)) + ggtitle("PCA after removing batch effect")

select_heatmap <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:25]
dds_as_df <- as.data.frame(colData(dds)[,c("batch", "stimulation")])
pheatmap(assay(vsd)[select_heatmap,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=dds_as_df, main = "vsd data after removing batch effect")

#Complex heatmap
#get top DEG
genes <- res[order(res$padj, -abs(res$log2FoldChange)),] %>%
  head(20) %>%
  rownames()
heatmap <- assay(vsd)[genes,]

#scale counts for visualization
heatmap <- t(scale(t(heatmap)))

# Plot as heatmap
heatmat.top20 <- ComplexHeatmap::Heatmap(heatmap,
                        cluster_rows = TRUE, 
                        cluster_columns = FALSE, 
                        col = colprgn,
                        column_title = "Heatmap top 20 significant genes")
heatmat.top20
tiff("../output/4days/heatmat.top20-ctrl-vs-swaps.tiff", width = 8, height = 7, units = "in", res = 300)
plot(heatmat.top20)
dev.off()

# Heatmap of specific genes
# List of specific genes you want to visualize
genes_of_interest <- c("Col1a1", "Col1a2","Col6a1", "Col6a2", "Col6a3", "Col6a4", "Col6a5", "Col6a6", "Acta2", "Ccl2", "Tagln", "Il6", "Lif", "Fn1")

# Ensure these genes are present in the vsd object
genes_present <- genes_of_interest[genes_of_interest %in% rownames(assay(vsd))]

# Subset the vsd object to include only the genes of interest
heatmap <- assay(vsd)[genes_present,]

# Scale the counts for better visualization in the heatmap
heatmap <- t(scale(t(heatmap)))

# Plot the heatmap
ComplexHeatmap::Heatmap(heatmap,
                        cluster_rows = TRUE, 
                        cluster_columns = FALSE, 
                        column_title = "Heatmap for Selected Genes")

## Heatmap of all DEG
res1 <- results(dds, contrast = c("stimulation", "TGFb->IL1a", "ctrl")) %>%
  as.data.frame() %>%
  na.omit() %>%
  filter(padj < 0.1 & abs(log2FoldChange) > 0.6)

res2 <- results(dds, contrast = c("stimulation", "IL1a->TGFb", "ctrl")) %>%
  as.data.frame() %>%
  na.omit() %>%
  filter(padj < 0.1 & abs(log2FoldChange) > 0.6)

res3 <- results(dds, contrast = c("stimulation", "TGFb", "ctrl")) %>%
  as.data.frame() %>%
  na.omit() %>%
  filter(padj < 0.1 & abs(log2FoldChange) > 0.6)

res4 <- results(dds, contrast = c("stimulation", "mix", "ctrl")) %>%
  as.data.frame() %>%
  na.omit() %>%
  filter(padj < 0.1 & abs(log2FoldChange) > 0.6)

genes.DE <- rbind(res1, res2, res3, res4)
genes.DE.names <- rownames(genes.DE)

# Remove duplicate rows (for overlapping genes across contrasts)
genes.DE <- genes.DE.names[genes.DE.names %in% rownames(assay(vsd))]
genes.DE6 <- genes.DE

# Subset the vsd object to include only the genes of interest
heatmap <- assay(vsd)[genes.DE,]

# Scale the counts for better visualization in the heatmap
heatmap <- t(scale(t(heatmap)))

# Plot the heatmap
pdf("../output/4days/4days_0.1_ctrl-mix-TGFbswap-IL1aswap-TGFb_Heatmap-of-DEG.pdf")
ComplexHeatmap::Heatmap(heatmap,
                        cluster_rows = TRUE, 
                        cluster_columns = FALSE, 
                        column_names_gp = gpar(fontsize = 6),
                        row_names_gp = gpar(fontsize = 1.5),
                        col = colprgn,
                        column_title = "ctrl vs TGFb, swaps and mix - Heatmap of DEG")
dev.off()


## GO condition per column dotplot

TGFbswap <- res1
IL1aswap <- res2
TGFb <- res3
mix <- res4

# GGPLOTS
#data wrangling for compareCluster()
TGFb$othergroup <- "TGFb"
TGFb <- rownames_to_column(TGFb, "Entrez")
mix$othergroup <- "mix"
mix <- rownames_to_column(mix, "Entrez")
TGFbswap$othergroup <- "TGFb->IL1a"
TGFbswap <- rownames_to_column(TGFbswap, "Entrez")
IL1aswap$othergroup <- "IL1a->TGFb"
IL1aswap <- rownames_to_column(IL1aswap, "Entrez")
df <- TGFb %>%
  full_join(mix) %>%
  full_join(TGFbswap) %>%
  full_join(IL1aswap) %>%
  mutate(diff.expr. = ifelse(log2FoldChange > 0, "up", "down"))

#BP
ccgo.bp <- compareCluster(Entrez~diff.expr.+othergroup, data = df, fun = "enrichGO", OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont="BP")
#data wrangling and plotting with ggplot
ccgo.bp.df <- as.data.frame(ccgo.bp) 

#MF
ccgo.mf <- compareCluster(Entrez~diff.expr.+othergroup, data = df, fun = "enrichGO", OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont="MF")
ccgo.mf.df <- as.data.frame(ccgo.mf) %>%
  filter(Description %in% c("protein-macromolecule adaptor activity", "glycosaminoglycan binding", "growth factor activity", "protein tyrosine kinase activator activity", "extracellular matrix structural constituent", "extracellular matrix structural constituent conferring compression resistance", "ATPase binding", "Wnt-protein binding", "proteoglycan binding", "cytokine activity", "chemokine activity", "chemokine receptor binding", "platelet-derived growth factor binding", "growth factor binding", "extracellular matrix structural constituent conferring tensile strength", "filamin binding", "MMP activity", "fibronectin binding", "integrin binding", "cell adhesion molecule binding", "collagen binding", "SMAD binding", "cytokine receptor binding"))

#instead of the GeneRatio, use a continuous value and change order of columns by a factor to define levels (orders of columns)
ccgo.mf.df.dec <- ccgo.mf.df %>%
  separate(GeneRatio, into = c("numerator", "denominator"), sep = "/") %>%
  mutate(GeneRatioDec = as.numeric(numerator) / as.numeric(denominator)) %>%
  select(-numerator, -denominator) %>%
  mutate(othergroup = factor(othergroup, levels = c("TGFb", "mix", "TGFb->IL1a", "IL1a->TGFb")))

#ggplot
ccgo.mg.ggplot <- ggplot(
  ccgo.mf.df.dec, 
  aes(
    x = diff.expr., 
    y = reorder(str_wrap(Description, width = 30), -p.adjust), 
    size = GeneRatioDec, 
    color = diff.expr., 
    alpha = -log10(p.adjust)  # Intensity based on -log10(p.adjust)
  )) +
  geom_point() +  # Color by diff.expr., alpha by p.adjust for intensity
  facet_wrap(~othergroup, ncol = 4) + 
  scale_color_manual(values = c("down" = "#762A83", "up" = "#1B7837")) +  # PRGn purple and green
  scale_alpha_continuous(range = c(0.3, 1)) +  # Controls dot intensity based on p.adjust, higher is more intense
  guides(alpha = guide_legend(override.aes = list(size = 6)),  # bigger dots in the legend
         size = guide_legend(override.aes = list(alpha = 1)),   # opacity of dots in the legend
         color = "none") +  # Adjust color legend dots
  theme_minimal() + #frames
    theme(
    axis.text.x = element_text(hjust = 1),
    strip.background = element_rect(color = "black", fill = "white", size = 1),  # Frame around facets
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Border around each facet
    panel.spacing = unit(1, "lines")  # Adds space between facets
  ) +
  ggtitle("Molecular function") +
  labs(size = "Gene Ratio", alpha = "-log10(p.adjust)") +
  xlab("Differential Expression") +
  ylab("GO Term")

ccgo.mg.ggplot

```

# Session info
```{r sesison info}
sessionInfo()
```


