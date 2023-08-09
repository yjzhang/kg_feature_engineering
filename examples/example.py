# basic usage example for kgfe
import time

import numpy as np
import pandas as pd

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
# graph size: 15957 nodes, 176595 edges
t = time.time()
pr_results, top_nodes = kgfe.explanations.topic_pagerank(graph, topic_ids)
print('igraph pagerank time:', time.time() - t)
base_pr_results, base_top_nodes = kgfe.explanations.topic_pagerank(graph)
# timing for pr: 76 ms (also includes postprocessing)

# %prun pr_results, top_nodes = kgfe.explanations.topic_pagerank(graph, topic_ids)

# test degree sampling
import scipy.stats
degrees = graph.degree(topic_ids)
dist = scipy.stats.gaussian_kde(degrees)
nodes = [x.index for x in kgfe.graph_info.nodes_in_category(graph, 'Gene')]
samples = kgfe.graph_info.degree_sample(graph, nodes, 10, dist)
print(samples)


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
null_stats = kgfe.explanations.null_graph_stats(graph, 'Gene', 20, 100, method='shortest_paths', target_nodes=[1])
print('time for null stats (shortest_paths method):', time.time() - t)
# time: 1.4s

# test null model degree sampling
t = time.time()
null_stats_degree_sample = kgfe.explanations.null_graph_stats(graph, 'Gene', 20, 100,
        use_degree_sampling=True, input_id_set=topic_ids, target_nodes=[1])
print('time for null stats (distances method), degree sampling:', time.time() - t)
null_stats_degree_sample = pd.DataFrame(null_stats_degree_sample)
print('mean degree (with degree sampling):', null_stats_degree_sample.degree_mean.mean())
print('mean degree for topic ids:', np.mean(graph.degree(topic_ids)))

t = time.time()
null_stats = kgfe.explanations.null_graph_stats(graph, 'Gene', 20, 100, method='distances', target_nodes=[1])
print('time for null stats (distances method):', time.time() - t)
# time: 1.4s (about the same)

# steiner tree
t = time.time()
st = kgfe.explanations.steiner_tree(graph, topic_ids, method='takahashi')
print('time for steiner tree (Takahashi method):', time.time() - t)
print('nodes in steiner tree (Takahashi method):', len(st.vs))
# time: 3ms

for node in topic_ids:
    assert(st.vs.find(name=node))
    assert(len(st.vs.find(name=node).neighbors()) > 0)
assert(st.is_tree())

t = time.time()
st_sp = kgfe.explanations.steiner_tree(graph, topic_ids, method='shortest_paths')
print('time for steiner tree (Shortest paths method):', time.time() - t)
print('nodes in steiner tree (Shortest paths method):', len(st.vs))
# time: 3ms

for node in topic_ids:
    assert(st_sp.vs.find(name=node))
    assert(len(st_sp.vs.find(name=node).neighbors()) > 0)
assert(st_sp.is_tree())


t = time.time()
st_mehlhorn = kgfe.explanations.steiner_tree(graph, topic_ids, method='mehlhorn')
print('time for steiner tree (Mehlhorn method):', time.time() - t)
print('nodes in steiner tree (Mehlhorn method):', len(st.vs))
# time: 3ms

for node in topic_ids:
    assert(st_mehlhorn.vs.find(name=node))
    assert(len(st_mehlhorn.vs.find(name=node).neighbors()) > 0)
assert(st_mehlhorn.is_tree())

# TODO: compare the size of different steiner tree methods
# takahashi, mehlhorn, shortest paths
takahashi_sizes = []
shortest_path_sizes = []
mehlhorn_sizes = []
takahashi_times = []
shortest_path_times = []
mehlhorn_times = []
for i in range(50):
    rand_topic_ids = kgfe.graph_info.random_nodes(graph, 20)
    t0 = time.time()
    st1 = kgfe.explanations.steiner_tree(graph, rand_topic_ids, method='takahashi')
    takahashi_sizes.append(len(st1.vs))
    takahashi_times.append(time.time() - t0)
    t0 = time.time()
    st2 = kgfe.explanations.steiner_tree(graph, rand_topic_ids, method='shortest_paths')
    shortest_path_sizes.append(len(st2.vs))
    shortest_path_times.append(time.time() - t0)
    t0 = time.time()
    st3 = kgfe.explanations.steiner_tree(graph, rand_topic_ids, method='mehlhorn')
    mehlhorn_sizes.append(len(st3.vs))
    mehlhorn_times.append(time.time() - t0)
print()
print('steiner tree approximation stats with 20 randomly selected nodes:')
print('mean size of takahashi steiner tree:', sum(takahashi_sizes)/50.)
print('mean size of shortest paths steiner tree:', sum(shortest_path_sizes)/50.)
print('mean size of mehlhorn steiner tree:', sum(mehlhorn_sizes)/50.)

print('mean running time of takahashi steiner tree:', sum(takahashi_times)/50.)
print('mean running time of shortest paths steiner tree:', sum(shortest_path_times)/50.)
print('mean running time of mehlhorn steiner tree:', sum(mehlhorn_times)/50.)

print()

import networkx as nx
t = time.time()
nx_graph = kgfe.graph_info.df_to_networkx(df)
print('nx graph from df time:', time.time() - t)
# time: 12s
# this graph has 158959 edges (possibly removing duplicate edges between nodes?)
# note: networkx by default does not have duplicate edges between nodes (requres nx.MultiGraph class to have that)

t = time.time()
nx_pr_results = nx.pagerank(nx_graph, personalization={i: 1 for i in topic_ids}, alpha=0.7)
print('nx pagerank time:', time.time() - t)
nx_base_pr_results = nx.pagerank(nx_graph, alpha=0.7)
# timing for pr: 372 ms (does not include postprocessing)

# compare pr results
pr_diff_sum = 0
for k, v in pr_results.items():
    diff = v - nx_pr_results[k]
    pr_diff_sum += diff**2
print('sum of squared differences in personalized PR results between networkx and igraph:', pr_diff_sum)

base_pr_diff_sum = 0
for k, v in base_pr_results.items():
    diff = v - nx_base_pr_results[k]
    base_pr_diff_sum += diff**2
print('sum of squared differences in base PR results between networkx and igraph:', base_pr_diff_sum)



# rank correlation?
top_pr_node_names = [x['name'] for x in top_nodes]
ig_node_names = set(top_pr_node_names)
top_nx_pr_nodes = sorted(nx_pr_results.items(), key=lambda x: x[1], reverse=True)
top_nx_pr_node_names = [x[0] for x in top_nx_pr_nodes if x[0] in ig_node_names]
from scipy.stats import spearmanr
print('rank correlation between ig and nx personalized pagerank results: ', spearmanr(top_pr_node_names, top_nx_pr_node_names))

top_base_pr_node_names = [x['name'] for x in base_top_nodes]
base_ig_node_names = set(top_base_pr_node_names)
top_nx_base_pr_nodes = sorted(nx_base_pr_results.items(), key=lambda x: x[1], reverse=True)
top_nx_base_pr_node_names = [x[0] for x in top_nx_base_pr_nodes if x[0] in base_ig_node_names]
print('rank correlation between ig and nx base pagerank results: ', spearmanr(top_base_pr_node_names, top_nx_base_pr_node_names))


print('rank correlation between personalized and base pagerank results for igraph:', spearmanr(top_pr_node_names, [x for x in top_base_pr_node_names if x in ig_node_names]))
print('rank correlation between personalized and base pagerank results for nx:', spearmanr(top_nx_pr_node_names, [x for x in top_nx_base_pr_node_names if x in ig_node_names]))

print()


t0 = time.time()
for i in range(20000):
    for t in topic_ids:
        node = nx_graph[t]
print('nx random access time (20000 iterations):', time.time() - t0)
# time: 59ms


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
