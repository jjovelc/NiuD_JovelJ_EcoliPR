import pandas as pd, numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

T = 15
df = pd.read_csv("snpdist.tsv", sep="\t")
samples = df.iloc[:,0].tolist()
mat = df.iloc[:,1:].to_numpy(dtype=float)

# Convert symmetric matrix to condensed form
dists = squareform(mat)
Z = linkage(dists, method="complete")

# Assign clusters
clusters = fcluster(Z, t=T, criterion="distance")
pd.DataFrame({"isolate": samples,
              "cluster_complete": ["CL_"+str(c) for c in clusters]}) \
  .to_csv(f"clusters_complete_T{T}.tsv", sep="\t", index=False)

