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
# time: 3.5, mean distance: 36.52

t = time.time()
ot_distances = []
for i in range(50):
    genes_1 = kgfe.graph_info.random_nodes_in_category(graph, 'Gene', 10)
    genes_2 = kgfe.graph_info.random_nodes_in_category(graph, 'Gene', 20)
    emd = kgfe.ot_distance.ot_sinkhorn_unbalanced_cost(graph, genes_1, genes_2)
    ot_distances.append(emd)
print('time for 50 unbalanced OT-sinkhorn distance calculations:', time.time() - t)
print('mean EMD between two randomly selected lists of 10 genes and 20 genes, 50 iterations:', np.mean(ot_distances)) 
# time: 2.5s

t = time.time()
gene_lists = [kgfe.graph_info.random_nodes_in_category(graph, 'Gene', 20) for i in range(200)]
distance_matrix = kgfe.ot_distance.ot_balanced_emd_distance_matrix(graph, gene_lists, verbose=True)
print('time for OT exact-EMD distance matrix calculation for 200 lists of 20 genes:', time.time() - t)
print('mean EMD between gene sets:', np.mean(distance_matrix)) 
# time: 17.7s
# EMD mean: 2.442875



t = time.time()
gene_lists = [kgfe.graph_info.random_nodes_in_category(graph, 'Gene', 20) for i in range(200)]
distance_matrix = kgfe.ot_distance.ot_sinkhorn_distance_matrix(graph, gene_lists, verbose=True)
print('time for OT unbalanced Sinkhorn distance matrix calculation for 200 lists of 20 genes:', time.time() - t)
print('mean EMD between gene sets:', np.mean(distance_matrix)) 
# timing: 32s
# EMD: 7.12


t = time.time()
gene_lists = [kgfe.graph_info.random_nodes_in_category(graph, 'Gene', 20) for i in range(200)]
distance_matrix = kgfe.ot_distance.ot_balanced_sinkhorn_distance_matrix(graph, gene_lists, verbose=True)
print('time for OT balanced Sinkhorn distance matrix calculation for 200 lists of 20 genes, reg=1:', time.time() - t)
print('mean EMD between gene sets:', np.mean(distance_matrix)) 
# timing: 26s
# EMD: 3.02


"""
t = time.time()
gene_lists = [kgfe.graph_info.random_nodes_in_category(graph, 'Gene', 20) for i in range(200)]
distance_matrix = kgfe.ot_distance.ot_balanced_sinkhorn_distance_matrix(graph, gene_lists, verbose=True, reg=0.5)
print('time for OT balanced Sinkhorn distance matrix calculation for 200 lists of 20 genes, reg=0.5:', time.time() - t)
print('mean EMD between gene sets:', np.mean(distance_matrix)) 
# timing: 40s
# EMD: 2.66

t = time.time()
gene_lists = [kgfe.graph_info.random_nodes_in_category(graph, 'Gene', 20) for i in range(200)]
distance_matrix = kgfe.ot_distance.ot_balanced_sinkhorn_distance_matrix(graph, gene_lists, verbose=True, reg=0.25)
print('time for OT balanced Sinkhorn distance matrix calculation for 200 lists of 20 genes, reg=0.25:', time.time() - t)
print('mean EMD between gene sets:', np.mean(distance_matrix)) 
# timing: 109s
# EMD: 2.48

t = time.time()
gene_lists = [kgfe.graph_info.random_nodes_in_category(graph, 'Gene', 20) for i in range(200)]
distance_matrix = kgfe.ot_distance.ot_balanced_sinkhorn_distance_matrix(graph, gene_lists, verbose=True, method='greenkhorn', reg=0.25)
print('time for OT balanced Greenkhorn distance matrix calculation for 200 lists of 20 genes, reg=0.25:', time.time() - t)
print('mean EMD between gene sets:', np.mean(distance_matrix)) 
# timing: 381s
# EMD: 2.45

t = time.time()
gene_lists = [kgfe.graph_info.random_nodes_in_category(graph, 'Gene', 20) for i in range(200)]
distance_matrix = kgfe.ot_distance.ot_balanced_sinkhorn_distance_matrix(graph, gene_lists, verbose=True, reg=0.1)
print('time for OT balanced Sinkhorn distance matrix calculation for 200 lists of 20 genes, reg=0.1:', time.time() - t)
print('mean EMD between gene sets:', np.mean(distance_matrix)) 
# timing: 603s
# EMD: 2.41


# TODO: Screenkhorn is not implemented
t = time.time()
gene_lists = [kgfe.graph_info.random_nodes_in_category(graph, 'Gene', 20) for i in range(200)]
distance_matrix = kgfe.ot_distance.ot_balanced_sinkhorn_distance_matrix(graph, gene_lists, verbose=True, method='screenkhorn')
print('time for OT balanced Screenkhorn distance matrix calculation for 200 lists of 20 genes:', time.time() - t)
print('mean EMD between gene sets:', np.mean(distance_matrix)) 
"""
