# optimal transport to determine the distance between node sets
import time

import numpy as np

import kgfe

# timing ran on intel i7-9700 CPU @ 3.00GHz with 64gb memory

df = kgfe.load_graph('reactome_genes_chems.csv.gz')

t = time.time()
graph = kgfe.df_to_graph(df)
graph.simplify()
print('igraph graph from df time:', time.time() - t)
# timing: 12s

topic_ids = ['NCBIGene::5972',
        'NCBIGene::958',
        'NCBIGene::100', 'NCBIGene::8797', 'NCBIGene::26762']

# TODO: get some known sets of genes/outliers, and get some randomly chosen groups of genes as well
t = time.time()
ot_distances = []
for i in range(50):
    genes_1 = kgfe.graph_info.random_nodes_in_category(graph, 'Gene', 10)
    genes_2 = kgfe.graph_info.random_nodes_in_category(graph, 'Gene', 10)
    ot_matrix, emd = kgfe.ot_distance.ot_table(graph, genes_1, genes_2)
    ot_distances.append(emd)
print('time for 50 OT distance calculations:', time.time() - t)
print('mean EMD between two randomly selected lists of 10 genes, 50 iterations:', np.mean(ot_distances)) 
