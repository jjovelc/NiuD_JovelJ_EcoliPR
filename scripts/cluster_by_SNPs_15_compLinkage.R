library(data.table)

setwd('/Users/juanjovel/Library/CloudStorage/OneDrive-UniversityofCalgary/jj/UofC/data_analysis/dongyanNiu/mawra/clustering_phylogenetics_human_Ecoli')

T <- 15                               # threshold
infile <- "snpdist_DIST.tsv"      # <-- set your filename

# Read distance matrix
dt <- fread(infile)
mat <- as.matrix(dt[, -1, with = FALSE])
rownames(mat) <- dt[[1]]
mode(mat) <- "numeric"

# Convert to 'dist' object
d <- as.dist(mat)

# Complete-linkage hierarchical clustering
hc <- hclust(d, method = "complete")

# Cut tree at SNP threshold
clusters <- cutree(hc, h = T)

# Save cluster membership
cl <- data.table(isolate = names(clusters),
                 cluster_complete = paste0("CL_", clusters))
fwrite(cl, sprintf("clusters_complete_T%d.tsv", T), sep = "\t")

# Summary table (cluster sizes)
cl[, .N, by = cluster_complete][order(-N)] |>
  fwrite(sprintf("clusters_complete_T%d_summary.tsv", T), sep="\t")

