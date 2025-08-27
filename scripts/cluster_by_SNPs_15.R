# files from snp-dists:
#   snpdist_<PREFIX>.tsv   (wide matrix with header, top-left blank due to -b)
#   snpdist_<PREFIX>_long.tsv  (optional long format)

library(data.table)
library(igraph)

setwd('/Users/juanjovel/Library/CloudStorage/OneDrive-UniversityofCalgary/jj/UofC/data_analysis/dongyanNiu/mawra/clustering_phylogenetics_human_Ecoli')

T <- 15                               # threshold
infile <- "snpdist_DIST.tsv      # <-- set your filename
dt <- fread(infile)

# Build numeric distance matrix
mat <- as.matrix(dt[, -1, with = FALSE])  # drop first col (row names)
rownames(mat) <- dt[[1]]
mode(mat) <- "numeric"

## -------- Option A: SINGLE-LINKAGE (graph components) --------
adj <- mat <= T
diag(adj) <- FALSE
g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
cc <- components(g)$membership
cl_single <- data.table(isolate = names(cc), cluster_single = paste0("SL_", cc))

# Edge list (pairs within threshold) – useful for Cytoscape/Gephi
edges <- which(adj & upper.tri(adj), arr.ind = TRUE)
edge_df <- data.table(
  from = rownames(mat)[edges[,1]],
  to   = colnames(mat)[edges[,2]],
  dist = mat[edges]
)

## -------- Option B: COMPLETE-LINKAGE (all pairs ≤ T) --------
d <- as.dist(mat)
hc <- hclust(d, method = "complete")
grp <- cutree(hc, h = T)
cl_complete <- data.table(isolate = names(grp), cluster_complete = paste0("CL_", grp))

## -------- Outputs --------
cl <- merge(cl_single, cl_complete, by = "isolate", all = TRUE)
fwrite(cl,      sprintf("clusters_T%d.tsv", T), sep = "\t")
fwrite(edge_df, sprintf("edges_T%d.tsv", T),    sep = "\t")

# Quick summaries
cl[, .N, by = cluster_single][order(-N)][, fwrite(.SD, sprintf("clusters_T%d_single_summary.tsv", T), sep="\t")]
cl[, .N, by = cluster_complete][order(-N)][, fwrite(.SD, sprintf("clusters_T%d_complete_summary.tsv", T), sep="\t")]

