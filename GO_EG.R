# Libraries
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(pathview)
library(svglite)

directory <- "../Output/EpilepsyGenes/GO_KEGG"

file <- "extemp_intemp.csv"

# Load Data
data <- read.csv(paste(directory, file, sep = "/"), header = TRUE)
head(data)
df_data <- as.data.frame(data)
head(df_data)
name <- strsplit(file, split = "\\.")[[1]][1]
name

# Create Folder
path <- paste("../Output/EpilepsyGenes/GO_KEGG", name, sep = "/")
dir.create(path)

# Dictionary
genes_info <- list()

# Gene List
gene_list <- df_data$Genes
head(gene_list)
genes_info[["Total"]] <- nrow(df_data)

# GO Enrichment Analysis - Upregulated Genes

# Perform GO enrichment - MF - Molecular Function
ego1 <- enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "MF",  pAdjustMethod = "BH", pvalueCutoff = 0.05)
go_df1 <- as.data.frame(ego1)
name1 <- paste(name, "GO_MF_UP.csv", sep="_")
write.csv(go_df1, file = paste(path, name1, sep="/"), row.names = FALSE)
genes_info[["GO_MF"]] <- nrow(go_df1)

# BarPlot
name2 <- paste(name, "GO_MF_barplot.svg", sep="_")
svg(filename = paste(path, name2, sep="/"), width = 6, height = 6)
barplot(ego1, showCategory = 10, title = "GO Enrichment - Barplot")
dev.off()

# DotPlot
name3 <- paste(name, "GO_MF_dotplot.svg", sep="_")
svg(filename = paste(path, name3, sep="/"), width = 6, height = 6)
dotplot(ego1, showCategory = 10, title = "GO Enrichment - Dotplot")
dev.off()

# Perform GO enrichment - BP - Biological Process
ego2 <- enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",  pAdjustMethod = "BH", pvalueCutoff = 0.05)
go_df2 <- as.data.frame(ego2)
name4 <- paste(name, "GO_BP.csv", sep="_")
write.csv(go_df2, file = paste(path, name4, sep="/"), row.names = FALSE)
genes_info[["GO_BP"]] <- nrow(go_df2)

# BarPlot
name5 <- paste(name, "GO_BP_barplot.svg", sep="_")
svg(filename = paste(path, name5, sep="/"), width = 6, height = 6)
barplot(ego2, showCategory = 10, title = "GO Enrichment - Barplot")
dev.off()

# DotPlot
name6 <- paste(name, "GO_BP_dotplot.svg", sep="_")
svg(filename = paste(path, name6, sep="/"), width = 6, height = 6)
dotplot(ego2, showCategory = 10, title = "GO Enrichment - Dotplot")
dev.off()

# Perform GO enrichment - CC - Cellular Component
ego3 <- enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "CC",  pAdjustMethod = "BH", pvalueCutoff = 0.05)
go_df3 <- as.data.frame(ego3)
name7 <- paste(name, "GO_CC.csv", sep="_")
write.csv(go_df3, file = paste(path, name7, sep="/"), row.names = FALSE)
genes_info[["GO_CC"]] <- nrow(go_df3)

# BarPlot
name8 <- paste(name, "GO_CC_barplot.svg", sep="_")
svg(filename = paste(path, name8, sep="/"), width = 6, height = 6)
barplot(ego3, showCategory = 10, title = "GO Enrichment - Barplot")
dev.off()

# DotPlot
name9 <- paste(name, "GO_CC_dotplot.svg", sep="_")
svg(filename = paste(path, name9, sep="/"), width = 6, height = 6)
dotplot(ego3, showCategory = 10, title = "GO Enrichment - Dotplot")
dev.off()

# KEGG Enrichment Analysis

# Convert gene symbols to Entrez IDs
gene_entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Perform KEGG enrichment
kegg_enrich <- enrichKEGG(gene = gene_entrez$ENTREZID, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff  = 0.2)
kegg_df <- as.data.frame(kegg_enrich)
name10 <- paste(name, "KEGG.csv", sep="_")
write.csv(kegg_df, file = paste(path, name10, sep="/"), row.names = FALSE)
genes_info[["KEGG"]] <- nrow(kegg_df)

# BarPlot
name11 <- paste(name, "KEGG_barplot.svg", sep="_")
svg(filename = paste(path, name11, sep="/"), width = 6, height = 6)
barplot(kegg_enrich, showCategory = 10, title = "KEGG Enrichment - Barplot")
dev.off()

# DotPlot
name12 <- paste(name, "KEGG_dotplot.svg", sep="_")
svg(filename = paste(path, name12, sep="/"), width = 6, height = 6)
dotplot(kegg_enrich, showCategory = 10, title = "KEGG Enrichment - Dotplot")
dev.off()

# Save gene stats
df <- data.frame(Key = names(genes_info), Value = unlist(genes_info), stringsAsFactors = FALSE)
name13 <- paste(name, "GeneStats.csv", sep="_")
write.csv(df, file = paste(path, name13, sep="/"), row.names = FALSE)


  




