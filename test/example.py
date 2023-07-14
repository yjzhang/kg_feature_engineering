# basic usage example for kgfe
import time

import kgfe

df = kgfe.load_graph('reactome_genes_chems.csv.gz')
t = time.time()
graph = kgfe.df_to_graph(df)
print('igraph graph from df time:', time.time() - t)
topic_ids = ['NCBIGene::5972',
        'NCBIGene::958',
        'NCBIGene::100', 'NCBIGene::8797', 'NCBIGene::26762']
t = time.time()
pr_results, top_nodes = kgfe.explanations.topic_pagerank(graph, topic_ids)
print('igraph pagerank time:', time.time() - t)
# graph size: 15957 nodes, 176595 edges
# timing for pr: 53 ms (also includes postprocessing)
#t = time.time()
#shortest_path_lengths = graph.distances()
#print('time for all-pairs shortest paths igraph:', time.time() - t)
# timing for shortest paths: 19s - 27s
t = time.time()
hypergeom_scores = kgfe.explanations.hypergeom_test(graph, topic_ids, query_category='Gene')
print('time for hypergeom:', time.time() - t)
t = time.time()
node_stats = kgfe.explanations.graph_node_stats(graph, topic_ids)
print('time for node stats:', time.time() - t)

# null graph statistics
t = time.time()
null_stats = kgfe.explanations.null_graph_stats(graph, 'Gene', 10, 100)
print('time for null stats:', time.time() - t)

import networkx as nx
t = time.time()
nx_graph = kgfe.graph_info.df_to_networkx(df)
print('nx graph from df time:', time.time() - t)
t = time.time()
nx_pr_results = nx.pagerank(nx_graph, personalization={i: 1 for i in topic_ids})
print('nx pagerank time:', time.time() - t)

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

print('ig pairwise distance:', node_stats['average_pairwise_distance'], 'nx pairwise distance:', nx_node_stats['average_pairwise_distance'])
print('difference in pairwise distances:', node_stats['average_pairwise_distance'] - nx_node_stats['average_pairwise_distance'])

t = time.time()
nx_null_stats = kgfe.explanations.null_graph_stats_networkx(nx_graph, 'Gene', 10, 100)
print('time for null stats networkx:', time.time() - t)

# this graph has 158959 edges (possibly removing duplicate edges between nodes?)
# note: networkx by default does not have duplicate edges between nodes (requres nx.MultiGraph class to have that)
# timing for pr: 137 ms (does not include postprocessing)
#shortest_path_lengths = list(nx.shortest_path_length(nx_graph))
# timing for shortest paths: gets killed on 16gb of memory (takes at least many minutes)

# also try doing a hypergeometric test
