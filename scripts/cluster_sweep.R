#!/usr/bin/env Rscript

# Load libraries
library(reshape2)
library(dplyr)


setwd('/Users/juanjovel/Library/CloudStorage/OneDrive-UniversityofCalgary/jj/UofC/data_analysis/dongyanNiu/mawra/clustering_phylogenetics_human_Ecoli')

# ---- INPUT ----
# Replace with your SNP distance TSV file (from snp-dists)
dist_file <- "snpdist_DIST.tsv"

# ---- READ DISTANCE MATRIX ----
dist_mat <- as.matrix(read.table(dist_file, header=TRUE, row.names=1, check.names=FALSE))

dim(dist_mat)        # should be 114 x 114
length(rownames(dist_mat))
length(unique(rownames(dist_mat)))

# Convert to "dist" object

# Clean names (strip extensions, spaces)
clean_names <- gsub("\\.fasta$", "", rownames(dist_mat))
clean_names <- gsub("\\s+", "_", clean_names)
rownames(dist_mat) <- clean_names
colnames(dist_mat) <- clean_names

dist_obj <- as.dist(dist_mat)

# ---- HIERARCHICAL CLUSTERING (complete linkage) ----
hc <- hclust(dist_obj, method="complete")

# ---- SWEEP THRESHOLDS ----
thresholds <- seq(5, 100, by=5)  # SNP thresholds to test
results <- data.frame()

for (t in thresholds) {
  clusters <- cutree(hc, h=t)
  cluster_sizes <- table(clusters)

  res <- data.frame(
    threshold = t,
    n_clusters = length(cluster_sizes),
    mean_size = mean(cluster_sizes),
    median_size = median(cluster_sizes),
    max_size = max(cluster_sizes),
    singletons = sum(cluster_sizes == 1)
  )
  results <- rbind(results, res)
}

# ---- SAVE SUMMARY ----
write.table(results, "cluster_threshold_summary.tsv", sep="\t", row.names=FALSE, quote=FALSE)

# ---- OPTIONAL: make plots ----
png("clusters_vs_threshold.png", width=800, height=600)
plot(results$threshold, results$n_clusters, type="b", pch=19,
     xlab="SNP threshold", ylab="Number of clusters",
     main="Number of clusters vs SNP threshold (complete linkage)")
dev.off()

png("avg_size_vs_threshold.png", width=800, height=600)
plot(results$threshold, results$mean_size, type="b", pch=19,
     xlab="SNP threshold", ylab="Average cluster size",
     main="Cluster size vs SNP threshold")
dev.off()

cat("✅ Done. See 'cluster_threshold_summary.tsv' and plots.\n")

