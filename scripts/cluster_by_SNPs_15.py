import pandas as pd, numpy as np, networkx as nx
T = 15
df = pd.read_csv("snpdist_<PREFIX>.tsv", sep="\t")
samples = df.iloc[:,0].tolist()
mat = df.iloc[:,1:].to_numpy(dtype=float)
G = nx.Graph()
G.add_nodes_from(samples)
for i in range(len(samples)):
    for j in range(i+1, len(samples)):
        if mat[i,j] <= T:
            G.add_edge(samples[i], samples[j], weight=mat[i,j])
clusters = {n: f"SL_{k}" for k, comp in enumerate(nx.connected_components(G), 1) for n in comp}
pd.DataFrame({"isolate": samples, "cluster_single": [clusters[s] for s in samples]}).to_csv(f"clusters_T{T}.tsv", sep="\t", index=False)

