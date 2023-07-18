# basic usage example for kgfe
import time

import kgfe

# timing ran on intel i7-9700 CPU @ 3.00GHz with 64gb memory

df = kgfe.load_graph('reactome_genes_chems.csv.gz')
t = time.time()
graph = kgfe.df_to_graph(df)
print('igraph graph from df time:', time.time() - t)
# timing: 12s
topic_ids = ['NCBIGene::5972',
        'NCBIGene::958',
        'NCBIGene::100', 'NCBIGene::8797', 'NCBIGene::26762']
# graph size: 15957 nodes, 176595 edges
t = time.time()
pr_results, top_nodes = kgfe.explanations.topic_pagerank(graph, topic_ids)
print('igraph pagerank time:', time.time() - t)
# timing for pr: 76 ms (also includes postprocessing)

# %prun pr_results, top_nodes = kgfe.explanations.topic_pagerank(graph, topic_ids)

# TODO: time for random node id access
t0 = time.time()
for i in range(20000):
    for t in topic_ids:
        node = graph.vs.find(name=t)
print('igraph random access time (20000 iterations):', time.time() - t0)
# time: 183ms

#t = time.time()
#shortest_path_lengths = graph.distances()
#print('time for all-pairs shortest paths igraph:', time.time() - t)
# timing for shortest paths: 19s - 27s
t = time.time()
hypergeom_scores = kgfe.explanations.hypergeom_test(graph, topic_ids, query_category='Gene')
print('time for hypergeom:', time.time() - t)
# time: 46ms

t = time.time()
node_stats = kgfe.explanations.graph_node_stats(graph, topic_ids)
print('time for node stats:', time.time() - t)
# time: 2ms

# null graph statistics
t = time.time()
null_stats = kgfe.explanations.null_graph_stats(graph, 'Gene', 20, 100)
print('time for null stats:', time.time() - t)
# time: 1.4s

# steiner tree
t = time.time()
st = kgfe.explanations.steiner_tree(graph, topic_ids, method='takahashi')
print('time for steiner tree (Takahashi method):', time.time() - t)
print('nodes in steiner tree (Takahashi method):', len(st.vs))
# time: 3ms

for node in topic_ids:
    assert(st.vs.find(name=node))
    assert(len(st.vs.find(name=node).neighbors()) == 1)
assert(st.is_tree())

t = time.time()
st_mehlhorn = kgfe.explanations.steiner_tree(graph, topic_ids, method='mehlhorn')
print('time for steiner tree (Mehlhorn method):', time.time() - t)
print('nodes in steiner tree (Mehlhorn method):', len(st.vs))
# time: 3ms

for node in topic_ids:
    assert(st_mehlhorn.vs.find(name=node))
    assert(len(st_mehlhorn.vs.find(name=node).neighbors()) == 1)
assert(st_mehlhorn.is_tree())

print()

import networkx as nx
t = time.time()
nx_graph = kgfe.graph_info.df_to_networkx(df)
print('nx graph from df time:', time.time() - t)
# time: 12s
# this graph has 158959 edges (possibly removing duplicate edges between nodes?)
# note: networkx by default does not have duplicate edges between nodes (requres nx.MultiGraph class to have that)

t = time.time()
nx_pr_results = nx.pagerank(nx_graph, personalization={i: 1 for i in topic_ids})
print('nx pagerank time:', time.time() - t)
# timing for pr: 372 ms (does not include postprocessing)

t0 = time.time()
for i in range(20000):
    for t in topic_ids:
        node = nx_graph[t]
print('nx random access time (20000 iterations):', time.time() - t0)
# time: 59ms

# compare pr results
pr_diff_sum = 0
for k, v in pr_results.items():
    diff = v - nx_pr_results[k]
    pr_diff_sum += diff**2
print('sum of squared differences in PR results between networkx and igraph:', pr_diff_sum)

t = time.time()
nx_hypergeom_scores = kgfe.explanations.hypergeom_test_networkx(nx_graph, topic_ids, query_category='Gene')
print('time for hypergeom networkx:', time.time() - t)
hg_diff_sum = 0
for k, v in hypergeom_scores.items():
    diff = v[0] - nx_hypergeom_scores[k][0]
    hg_diff_sum += diff**2
print('sum of squared differences in hypergeom results between networkx and igraph:', hg_diff_sum)

# get node stats
t = time.time()
nx_node_stats = kgfe.explanations.graph_node_stats_networkx(nx_graph, topic_ids)
print('time for node stats networkx:', time.time() - t)
# time: 7ms

print('ig pairwise distance:', node_stats['average_pairwise_distance'], 'nx pairwise distance:', nx_node_stats['average_pairwise_distance'])
print('difference in pairwise distances:', node_stats['average_pairwise_distance'] - nx_node_stats['average_pairwise_distance'])

t = time.time()
nx_null_stats = kgfe.explanations.null_graph_stats_networkx(nx_graph, 'Gene', 20, 100)
print('time for null stats networkx:', time.time() - t)
# time: 6.6s

import pandas as pd
null_stats = pd.DataFrame(null_stats)
nx_null_stats = pd.DataFrame(nx_null_stats)
print('ig average pairwise distance for random genes:', null_stats.average_pairwise_distance.mean())
print('nx average pairwise distance for random genes:', nx_null_stats.average_pairwise_distance.mean())

t = time.time()
nx_steiner_tree = nx.approximation.steiner_tree(nx_graph, topic_ids, method='mehlhorn')
print('nx steiner tree time (mehlhorn):', time.time() - t)
print('size of nx steiner tree (mehlhorn):', len(nx_steiner_tree.nodes))
# time: 590ms

#shortest_path_lengths = list(nx.shortest_path_length(nx_graph))
# timing for shortest paths: gets killed on 16gb of memory (takes at least many minutes)

# also try doing a hypergeometric test
