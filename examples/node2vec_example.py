import time

import numpy as np

import kgfe

# node2vec - start out by computing random walks, and then running word2vec

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

# random walks, parallel random walks
t = time.time()
random_walks = kgfe.node2vec.random_walks_igraph(graph, 10, 25)
print('time for random walks:', time.time() - t)
print(len(random_walks), len(random_walks[0]))

t = time.time()
random_walks_parallel = kgfe.node2vec.random_walks_igraph_parallel(graph, 10, 25)
print('time for parallel random walks:', time.time() - t)
print(len(random_walks_parallel), len(random_walks_parallel[0]))

# use word2vec to generate the actual embeddings
t = time.time()
model1 = kgfe.node2vec.run_word2vec(random_walks)
print('time for word2vec:', time.time() - t)
