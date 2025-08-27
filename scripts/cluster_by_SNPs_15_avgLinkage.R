# Load distance matrix
dist_mat <- as.matrix(read.table("snpdist_DIST.tsv",
                                 header=TRUE, row.names=1, check.names=FALSE))
dist_obj <- as.dist(dist_mat)

# Average-linkage hierarchical clustering
hc_avg <- hclust(dist_obj, method="average")

# Cut the tree at SNP threshold = 15
clusters_avg <- cutree(hc_avg, h=15)

# Save cluster assignments
write.table(data.frame(Sample=names(clusters_avg),
                       Cluster=clusters_avg),
            file="clusters_avg_15.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# Quick look at cluster sizes
table(clusters_avg)

