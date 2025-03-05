# Libraries
library(ggplot2)
library(sva)

# Load Data: Count matrix and metadata
counts <- read.csv("../Data/final_counts_matrix_bulk.csv", row.names = 1)  # Your raw count data
head(counts)
coldata <- read.csv("../Data/expDes_bulk.csv", row.names = 1)      # Metadata with conditions/samples info
head(coldata)
# Ensure row names of coldata match column names of counts
all(rownames(coldata) == colnames(counts))  # Should return TRUE

#Apply Filters & Add 1 to all counts
counts <- counts[rowSums(counts) >= 40, ]
counts <- counts + 1

# PCA Raw
log_raw_counts <- log2(counts)

# Perform PCA before normalization
pca_raw <- prcomp(t(log_raw_counts))
percent_variance <- round(100 * pca_raw$sdev^2 / sum(pca_raw$sdev^2), 1)

# Create PCA DataFrame
pca_df_raw <- as.data.frame(pca_raw$x)
pca_df_raw$batch <- coldata$batch
pca_df_raw$neurons <- coldata$neurons
# pca_df_raw$condition <- coldata$condition -bulk
# pca_df_raw$comb <- coldata$comb -bulk 

# Plot PCA 
# pca_bef_plot <- ggplot(pca_df_raw, aes(x = PC1, y = PC2)) +
#   geom_point(aes(color = condition, shape = batch), size = 3) +
#   geom_text(aes(label = substr(pca_df_raw$neurons, 1, 1)), vjust = -0.5, hjust = -0.2, check_overlap = FALSE) +
#   xlab(paste0("PC1 (", percent_variance[1], "%)")) +
#   ylab(paste0("PC2 (", percent_variance[2], "%)")) +
#   ggtitle("PCA on Raw Counts") +
#   theme_minimal()

# Bulk
pca_bef_plot <- ggplot(pca_df_raw, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = neurons, shape = batch), size = 3) +
  xlab(paste0("PC1 (", percent_variance[1], "%)")) +
  ylab(paste0("PC2 (", percent_variance[2], "%)")) +
  ggtitle("PCA on Raw Counts") +
  theme_minimal()

# Image
ggsave("../Output/BatchCorrection/PCA_Raw_bulk.svg", pca_bef_plot, width = 6, height = 6)


# Normalized bath-corrected counts
combatoutput <- ComBat_seq(counts = as.matrix(counts), batch = coldata$batch, group = coldata$neurons) 
dim(combatoutput)
head(combatoutput)

log_norm_counts <- log2(combatoutput + 1)
head(log_norm_counts)

# Perform PCA before normalization
pca_norm <- prcomp(t(log_norm_counts))
head(pca_norm$x)
percent_variance_norm <- round(100 * pca_norm$sdev^2 / sum(pca_norm$sdev^2), 1)

# Create PCA DataFrame
pca_df_norm <- as.data.frame(pca_norm$x)
pca_df_norm$batch <- coldata$batch
pca_df_norm$neurons <- coldata$neurons
# pca_df_norm$condition <- coldata$condition -bulk
# pca_df_norm$comb <- coldata$comb -bulk
# PCA Normalized

# Plot PCA
# pca_aft_plot <- ggplot(pca_df_norm, aes(x = PC1, y = PC2)) +
#   geom_point(aes(color = condition, shape = batch), size = 3) +
#   geom_text(aes(label = substr(pca_df_norm$neurons, 1, 1)), vjust = -0.5, hjust = -0.2, check_overlap = FALSE) +
#   xlab(paste0("PC1 (", percent_variance_norm[1], "%)")) +
#   ylab(paste0("PC2 (", percent_variance_norm[2], "%)")) +
#   ggtitle("PCA on Corrected Counts") +
#   theme_minimal()

#bulk 
pca_aft_plot <- ggplot(pca_df_norm, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = neurons, shape = batch), size = 3) +
  xlab(paste0("PC1 (", percent_variance_norm[1], "%)")) +
  ylab(paste0("PC2 (", percent_variance_norm[2], "%)")) +
  ggtitle("PCA on Corrected Counts") +
  theme_minimal()

# Image
ggsave("../Output/BatchCorrection/PCA_Norm_bulk.svg", pca_aft_plot, width = 6, height = 6)

combat_df <- as.data.frame(combatoutput)
write.csv(combat_df, "../Data/combat_corrected_counts_bulk.csv", row.names = TRUE)
