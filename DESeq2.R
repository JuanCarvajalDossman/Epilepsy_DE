# Libraries
library(ggplot2)
library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(apeglm)
library(msigdbr)
library(pathview)
library(svglite)

# Load Data: Count matrix and metadata
counts <- read.csv("../Data/combat_corrected_counts.csv", row.names = 1)  # Your raw count data
head(counts)
coldata <- read.csv("../Data/expDes.csv", row.names = 1)      # Metadata with conditions/samples info
head(coldata)
# Ensure row names of coldata match column names of counts
all(rownames(coldata) == colnames(counts))  # Should return TRUE

# Define Count Filter
# row_Sum <- read.csv("../Data/rowsum.csv")
# df<-data.frame(row_Sum)
# df <- data.frame(lapply(df, log2))
# Create histogram
# hist(as.numeric(unlist(df)), breaks = 200, main = "Custom Bin Histogram", xlab = "Values", ylab = "Frequency", col = "lightblue")

# Experiments Design
# experiments <- read.csv("../Data/experiments.csv", header = TRUE)      

###### Experiments Design ######
exp_name <- "comb"

# loop trough experiments
# for (i in 1:nrow(experiments)) 
# {

# cat("Processing iteration:", i, "\n")

# Get Experiment info and names
# exp_name <- experiments[i, 1]  # Esperimental Desing
# A_name <- experiments[i, 2]  # Control
# B_name <- experiments[i, 3]  # Treatment

# Create DESeq2 object
if (exp_name == "neurons") 
{
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ neurons)
} 
if (exp_name == "condition") 
{
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
} 
if (exp_name == "comb")
{
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ comb)
}
dds <- DESeq(dds)

# Coeficients
resultsNames(dds)

###### Levels #######
A_name <- "inref"        
B_name <- "intemp"

# Create Folder and suffix
suffix <- paste(exp_name, A_name, B_name, sep = "_")
path=paste("../Output/Results", suffix, sep = "_")
dir.create(path)

# Differential expression analysis
res <- results(dds, contrast = c(exp_name, A_name, B_name))
name1 <- paste(exp_name, A_name, B_name, "Result.csv", sep="_")
write.csv(as.data.frame(res), file = paste(path, name1, sep="/"))

#Shrinkage
coef_n <- paste(exp_name, B_name, "vs", A_name, sep="_") 
res <- lfcShrink(dds, coef=coef_n, type="apeglm")
name2 <- paste(exp_name, A_name, B_name, "Shrink.csv", sep="_")
write.csv(as.data.frame(res), file = paste(path, name2, sep="/"))

# List of # Genes
genes_info <- list()
# Number of genes
genes_info[["Total"]] <- nrow(res)

# Filter significant results (Adjusted p-value < 0.05 and |log2FoldChange| > 2)
res_sig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
name3 <- paste(exp_name, A_name, B_name, "Significant.csv", sep="_")
write.csv(as.data.frame(res_sig), file = paste(path, name3, sep="/"))

# Number of significant genes
genes_info[["Significant"]] <- nrow(res_sig)

# Upregulated genes
res_up <- subset(res_sig, log2FoldChange > 0)
name4 <- paste(exp_name, A_name, B_name, "Upregulated.csv", sep="_")
write.csv(as.data.frame(res_up), file = paste(path, name4, sep="/"))
genes_info[["Up"]] <- nrow(res_up)

# Downregulated genes
res_down <- subset(res_sig, log2FoldChange < 0)
name5 <- paste(exp_name, A_name, B_name, "Downregulated.csv", sep="_")
write.csv(as.data.frame(res_down), file = paste(path, name5, sep="/"))
genes_info[["down"]] <- nrow(res_down)

# Volcano Plot
name6 <- paste(exp_name, A_name, B_name, "Volcano.svg", sep="_")
svg(paste(path, name6, sep="/"), width = 6, height = 6)
EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'pvalue', pCutoff = 0.05, FCcutoff = 2)
dev.off() 

#MA Plot
name7 <- paste(exp_name, A_name, B_name, "MA.svg", sep="_")
svg(paste(path, name7, sep="/"), width = 6, height = 6)
title <- paste(exp_name, A_name, "vs", B_name, sep=" ")
plotMA(res, main = title, ylim = c(-10, 10))
points(
  res$baseMean[abs(res$log2FoldChange) > 2 & res$padj < 0.05], 
  res$log2FoldChange[abs(res$log2FoldChange) > 2 & res$padj < 0.05],
  col = "red",
  pch = 20
)
abline(h = c(-2, 2), col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Significant", "Not significant (log2FoldChange)","Not significant (padj)"), col = c("red", "blue", "gray"), pch = 16)
dev.off()

# Perform variance-stabilizing transformation
vsd <- vst(dds, blind = FALSE)
# Select the top 100 most significant genes for the heatmap
top_genes <- head(order(res$padj), 100)

# Heatmap
name8 <- paste(exp_name, A_name, B_name, "Heatmap.pdf", sep="_")
pdf(paste(path, name8, sep="/"), width = 20, height = 18)
pheatmap(assay(vsd)[top_genes, ], cluster_cols = TRUE, 
         main = "Top 100 Significant Genes Heatmap",
         fontsize_row = 10,         # Font size for row labels
         fontsize_col = 10,         # Font size for column labels
         labels_row = rownames(assay(vsd)[top_genes, ]),  # Row labels
         labels_col = colnames(assay(vsd)))
dev.off()

# GO Enrichment Analysis - Upregulated Genes

# Gene List Upregulated
gene_list_up <- rownames(res_up)

# Perform GO enrichment - MF - Molecular Function
ego1 <- enrichGO(gene = gene_list_up, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "MF",  pAdjustMethod = "BH", pvalueCutoff = 0.05)
go_df1 <- as.data.frame(ego1)
name9 <- paste(exp_name, A_name, B_name, "GO_MF_UP.csv", sep="_")
write.csv(go_df1, file = paste(path, name9, sep="/"), row.names = FALSE)
genes_info[["GO_UP_MF"]] <- nrow(go_df1)

# BarPlot
name10 <- paste(exp_name, A_name, B_name, "GO_MF_UP_barplot.svg", sep="_")
svg(filename = paste(path, name10, sep="/"), width = 6, height = 6)
barplot(ego1, showCategory = 10, title = "GO Enrichment UP - Barplot")
dev.off()

# DotPlot
name11 <- paste(exp_name, A_name, B_name, "GO_MF_UP_dotplot.svg", sep="_")
svg(filename = paste(path, name11, sep="/"), width = 6, height = 6)
dotplot(ego1, showCategory = 10, title = "GO Enrichment UP - Dotplot")
dev.off()

# Perform GO enrichment - BP - Biological Process
ego2 <- enrichGO(gene = gene_list_up, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",  pAdjustMethod = "BH", pvalueCutoff = 0.05)
go_df2 <- as.data.frame(ego2)
name12 <- paste(exp_name, A_name, B_name, "GO_UP_BP.csv", sep="_")
write.csv(go_df2, file = paste(path, name12, sep="/"), row.names = FALSE)
genes_info[["GO_UP_BP"]] <- nrow(go_df2)

# BarPlot
name13 <- paste(exp_name, A_name, B_name, "GO_BP_UP_barplot.svg", sep="_")
svg(filename = paste(path, name13, sep="/"), width = 6, height = 6)
barplot(ego2, showCategory = 10, title = "GO Enrichment UP - Barplot")
dev.off()

# DotPlot
name14 <- paste(exp_name, A_name, B_name, "GO_BP_UP_dotplot.svg", sep="_")
svg(filename = paste(path, name14, sep="/"), width = 6, height = 6)
dotplot(ego2, showCategory = 10, title = "GO Enrichment UP - Dotplot")
dev.off()

# Perform GO enrichment - CC - Cellular Component
ego3 <- enrichGO(gene = gene_list_up, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "CC",  pAdjustMethod = "BH", pvalueCutoff = 0.05)
go_df3 <- as.data.frame(ego3)
name15 <- paste(exp_name, A_name, B_name, "GO_CC_UP.csv", sep="_")
write.csv(go_df3, file = paste(path, name15, sep="/"), row.names = FALSE)
genes_info[["GO_UP_CC"]] <- nrow(go_df3)

# BarPlot
name16 <- paste(exp_name, A_name, B_name, "GO_CC_UP_barplot.svg", sep="_")
svg(filename = paste(path, name16, sep="/"), width = 6, height = 6)
barplot(ego3, showCategory = 10, title = "GO Enrichment UP - Barplot")
dev.off()

# DotPlot
name17 <- paste(exp_name, A_name, B_name, "GO_CC_UP_dotplot.svg", sep="_")
svg(filename = paste(path, name17, sep="/"), width = 6, height = 6)
dotplot(ego3, showCategory = 10, title = "GO Enrichment UP - Dotplot")
dev.off()

# GO Enrichment Analysis - Downregulated Genes

# Gene List Downregulated
gene_list_down <- rownames(res_down)

# Perform GO enrichment - MF - Molecular Function
ego4 <- enrichGO(gene = gene_list_down, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "MF",  pAdjustMethod = "BH", pvalueCutoff = 0.05)
go_df4 <- as.data.frame(ego4)
name18 <- paste(exp_name, A_name, B_name, "GO_MF_DOWN.csv", sep="_")
write.csv(go_df4, file = paste(path, name18, sep="/"), row.names = FALSE)
genes_info[["GO_DOWN_MF"]] <- nrow(go_df4)

# BarPlot
name19 <- paste(exp_name, A_name, B_name, "GO_MF_DOWN_barplot.svg", sep="_")
svg(filename = paste(path, name19, sep="/"), width = 6, height = 6)
barplot(ego4, showCategory = 10, title = "GO Enrichment DOWN - Barplot")
dev.off()

# DotPlot
name20 <- paste(exp_name, A_name, B_name, "GO_MF_DOWN_dotplot.svg", sep="_")
svg(filename = paste(path, name20, sep="/"), width = 6, height = 6)
dotplot(ego4, showCategory = 10, title = "GO Enrichment DOWN - Dotplot")
dev.off()

# Perform GO enrichment - BP - Biological Process
ego5 <- enrichGO(gene = gene_list_down, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",  pAdjustMethod = "BH", pvalueCutoff = 0.05)
go_df5 <- as.data.frame(ego5)
name21 <- paste(exp_name, A_name, B_name, "GO_DOWN_BP.csv", sep="_")
write.csv(go_df5, file = paste(path, name21, sep="/"), row.names = FALSE)
genes_info[["GO_DOWN_BP"]] <- nrow(go_df5)

# BarPlot
name22 <- paste(exp_name, A_name, B_name, "GO_BP_DOWN_barplot.svg", sep="_")
svg(filename = paste(path, name22, sep="/"), width = 6, height = 6)
barplot(ego5, showCategory = 10, title = "GO Enrichment DOWN - Barplot")
dev.off()

# DotPlot
name23 <- paste(exp_name, A_name, B_name, "GO_BP_DOWN_dotplot.svg", sep="_")
svg(filename = paste(path, name23, sep="/"), width = 6, height = 6)
dotplot(ego5, showCategory = 10, title = "GO Enrichment DOWN - Dotplot")
dev.off()

# Perform GO enrichment - CC - Cellular Component
ego6 <- enrichGO(gene = gene_list_down, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "CC",  pAdjustMethod = "BH", pvalueCutoff = 0.05)
go_df6 <- as.data.frame(ego6)
name24 <- paste(exp_name, A_name, B_name, "GO_CC_DOWN.csv", sep="_")
write.csv(go_df6, file = paste(path, name24, sep="/"), row.names = FALSE)
genes_info[["GO_DOWN_CC"]] <- nrow(go_df6)

# BarPlot
name25 <- paste(exp_name, A_name, B_name, "GO_CC_DOWN_barplot.svg", sep="_")
svg(filename = paste(path, name25, sep="/"), width = 6, height = 6)
barplot(ego6, showCategory = 10, title = "GO Enrichment DOWN - Barplot")
dev.off()

# DotPlot
name26 <- paste(exp_name, A_name, B_name, "GO_CC_DOWN_dotplot.svg", sep="_")
svg(filename = paste(path, name26, sep="/"), width = 6, height = 6)
dotplot(ego6, showCategory = 10, title = "GO Enrichment DOWN - Dotplot")
dev.off()

# GO Enrichment Analysis - All Genes

# Gene List All
gene_list_all <- rownames(res_sig)

# Perform GO enrichment - MF - Molecular Function
ego7 <- enrichGO(gene = gene_list_all, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "MF",  pAdjustMethod = "BH", pvalueCutoff = 0.05)
go_df7 <- as.data.frame(ego7)
name27 <- paste(exp_name, A_name, B_name, "GO_MF_ALL.csv", sep="_")
write.csv(go_df7, file = paste(path, name27, sep="/"), row.names = FALSE)
genes_info[["GO_ALL_MF"]] <- nrow(go_df7)

# BarPlot
name28 <- paste(exp_name, A_name, B_name, "GO_MF_ALL_barplot.svg", sep="_")
svg(filename = paste(path, name28, sep="/"), width = 6, height = 6)
barplot(ego7, showCategory = 10, title = "GO Enrichment ALL - Barplot")
dev.off()

# DotPlot
name29 <- paste(exp_name, A_name, B_name, "GO_MF_ALL_dotplot.svg", sep="_")
svg(filename = paste(path, name29, sep="/"), width = 6, height = 6)
dotplot(ego7, showCategory = 10, title = "GO Enrichment ALL - Dotplot")
dev.off()

# Perform GO enrichment - BP - Biological Process
ego8 <- enrichGO(gene = gene_list_all, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",  pAdjustMethod = "BH", pvalueCutoff = 0.05)
go_df8 <- as.data.frame(ego8)
name30 <- paste(exp_name, A_name, B_name, "GO_BP_ALL.csv", sep="_")
write.csv(go_df8, file = paste(path, name30, sep="/"), row.names = FALSE)

# BarPlot
name31 <- paste(exp_name, A_name, B_name, "GO_BP_ALL_barplot.svg", sep="_")
svg(filename = paste(path, name31, sep="/"), width = 6, height = 6)
barplot(ego8, showCategory = 10, title = "GO Enrichment ALL - Barplot")
dev.off()

# DotPlot
name32 <- paste(exp_name, A_name, B_name, "GO_BP_ALL_dotplot.svg", sep="_")
svg(filename = paste(path, name32, sep="/"), width = 6, height = 6)
dotplot(ego8, showCategory = 10, title = "GO Enrichment ALL - Dotplot")
dev.off()

# Perform GO enrichment - CC - Cellular Component
ego9 <- enrichGO(gene = gene_list_all, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "CC",  pAdjustMethod = "BH", pvalueCutoff = 0.05)
go_df9 <- as.data.frame(ego9)
name33 <- paste(exp_name, A_name, B_name, "GO_CC_ALL.csv", sep="_")
write.csv(go_df9, file = paste(path, name33, sep="/"), row.names = FALSE)
genes_info[["GO_ALL_CC"]] <- nrow(go_df9)

# BarPlot
name34 <- paste(exp_name, A_name, B_name, "GO_CC_ALL_barplot.svg", sep="_")
svg(filename = paste(path, name34, sep="/"), width = 6, height = 6)
barplot(ego9, showCategory = 10, title = "GO Enrichment ALL - Barplot")
dev.off()

# DotPlot
name35 <- paste(exp_name, A_name, B_name, "GO_CC_ALL_dotplot.svg", sep="_")
svg(filename = paste(path, name35, sep="/"), width = 6, height = 6)
dotplot(ego9, showCategory = 10, title = "GO Enrichment ALL - Dotplot")
dev.off()

# KEGG Enrichment Analysis - Upregulated Genes

# Convert gene symbols to Entrez IDs
gene_entrez_up <- bitr(gene_list_up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Perform KEGG enrichment
kegg_enrich_up <- enrichKEGG(gene = gene_entrez_up$ENTREZID, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff  = 0.2)
kegg_df_up <- as.data.frame(kegg_enrich_up)
name36 <- paste(exp_name, A_name, B_name, "KEGG_UP.csv", sep="_")
write.csv(kegg_df_up, file = paste(path, name36, sep="/"), row.names = FALSE)
genes_info[["KEGG_UP"]] <- nrow(kegg_df_up)

# BarPlot
name37 <- paste(exp_name, A_name, B_name, "KEGG_UP_barplot.svg", sep="_")
svg(filename = paste(path, name37, sep="/"), width = 6, height = 6)
barplot(kegg_enrich_up, showCategory = 10, title = "KEGG Enrichment UP - Barplot")
dev.off()

# DotPlot
name38 <- paste(exp_name, A_name, B_name, "KEGG_UP_dotplot.svg", sep="_")
svg(filename = paste(path, name38, sep="/"), width = 6, height = 6)
dotplot(kegg_enrich_up, showCategory = 10, title = "KEGG Enrichment UP - Dotplot")
dev.off()

# KEGG Enrichment Analysis - Downregulated Genes

# Convert gene symbols to Entrez IDs
gene_entrez_down <- bitr(gene_list_down, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Perform KEGG enrichment
kegg_enrich_down <- enrichKEGG(gene = gene_entrez_down$ENTREZID, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff  = 0.2)
kegg_df_down <- as.data.frame(kegg_enrich_down)
name39 <- paste(exp_name, A_name, B_name, "KEGG_DOWN.csv", sep="_")
write.csv(kegg_df_down, file = paste(path, name39, sep="/"), row.names = FALSE)
genes_info[["KEGG_DOWN"]] <- nrow(kegg_df_down)

# BarPlot
name40 <- paste(exp_name, A_name, B_name, "KEGG_DOWN_barplot.svg", sep="_")
svg(filename = paste(path, name40, sep="/"), width = 6, height = 6)
barplot(kegg_enrich_down, showCategory = 10, title = "KEGG Enrichment DOWN - Barplot")
dev.off()

# DotPlot
name41 <- paste(exp_name, A_name, B_name, "KEGG_DOWN_dotplot.svg", sep="_")
svg(filename = paste(path, name41, sep="/"), width = 6, height = 6)
dotplot(kegg_enrich_down, showCategory = 10, title = "KEGG Enrichment DOWN - Dotplot")
dev.off()

# KEGG Enrichment Analysis - All Genes

# Convert gene symbols to Entrez IDs
gene_entrez_all <- bitr(gene_list_all, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Perform KEGG enrichment
kegg_enrich_all <- enrichKEGG(gene = gene_entrez_all$ENTREZID, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff  = 0.2)
kegg_df_all <- as.data.frame(kegg_enrich_all)
name42 <- paste(exp_name, A_name, B_name, "KEGG_ALL.csv", sep="_")
write.csv(kegg_df_all, file = paste(path, name42, sep="/"), row.names = FALSE)
genes_info[["KEGG_ALL"]] <- nrow(kegg_df_all)

# BarPlot
name43 <- paste(exp_name, A_name, B_name, "KEGG_ALL_barplot.svg", sep="_")
svg(filename = paste(path, name43, sep="/"), width = 6, height = 6)
barplot(kegg_enrich_all, showCategory = 10, title = "KEGG Enrichment ALL - Barplot")
dev.off()

# DotPlot
name44 <- paste(exp_name, A_name, B_name, "KEGG_ALL_dotplot.svg", sep="_")
svg(filename = paste(path, name44, sep="/"), width = 6, height = 6)
dotplot(kegg_enrich_all, showCategory = 10, title = "KEGG Enrichment ALL - Dotplot")
dev.off()

# Save gene stats
df <- data.frame(Key = names(genes_info), Value = unlist(genes_info), stringsAsFactors = FALSE)
name45 <- paste(exp_name, A_name, B_name, "GeneStats.csv", sep="_")
write.csv(df, file = paste(path, name45, sep="/"), row.names = FALSE)

# }
  
  




