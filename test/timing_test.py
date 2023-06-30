import time
import cProfile

import networkx as nx
import kgfe

random_graph = nx.random_graphs.barabasi_albert_graph(100, 3)
t0 = time.time()
for i in range(100):
    paths_nx = dict(nx.shortest_path_length(random_graph))
print('shortest paths nx, BA, 100 nodes:', (time.time() - t0)/100)

t0 = time.time()
for i in range(100):
    paths_kgfe = kgfe.graph_utils.all_pairs_shortest_path_lengths(random_graph, list(random_graph.nodes))
print('shortest paths kgfe, BA, 100 nodes:', (time.time() - t0)/100)

t0 = time.time()
for i in range(100):
    paths_kgfe = kgfe.graph_utils.all_pairs_shortest_path_lengths_full(random_graph)
print('shortest paths-full kgfe, BA, 100 nodes:', (time.time() - t0)/100)



random_graph = nx.random_graphs.barabasi_albert_graph(1000, 3)
t0 = time.time()
for i in range(10):
    paths_nx = dict(nx.shortest_path_length(random_graph))
print('shortest paths nx, BA, 1000 nodes:', (time.time() - t0)/10)

t0 = time.time()
for i in range(10):
    paths_kgfe = kgfe.graph_utils.all_pairs_shortest_path_lengths(random_graph, list(random_graph.nodes))
print('shortest paths kgfe, BA, 1000 nodes:', (time.time() - t0)/10)

t0 = time.time()
for i in range(10):
    paths_kgfe = kgfe.graph_utils.all_pairs_shortest_path_lengths_full(random_graph)
print('shortest paths-full kgfe, BA, 1000 nodes:', (time.time() - t0)/10)


random_graph = nx.random_graphs.barabasi_albert_graph(10000, 3)
t0 = time.time()
paths_nx = dict(nx.shortest_path_length(random_graph))
print('shortest paths nx, BA, 10000 nodes:', (time.time() - t0))

t0 = time.time()
paths_kgfe = kgfe.graph_utils.all_pairs_shortest_path_lengths(random_graph, list(random_graph.nodes))
print('shortest paths kgfe, BA, 10000 nodes:', (time.time() - t0))



random_graph = nx.random_graphs.erdos_renyi_graph(1000, 0.5)
t0 = time.time()
paths_nx = dict(nx.shortest_path_length(random_graph))
print('shortest paths nx, ER, 1000 nodes:', (time.time() - t0))

t0 = time.time()
paths_kgfe = kgfe.graph_utils.all_pairs_shortest_path_lengths(random_graph, list(random_graph.nodes))
print('shortest paths kgfe, ER, 1000 nodes:', (time.time() - t0))

t0 = time.time()
paths_kgfe = kgfe.graph_utils.all_pairs_shortest_path_lengths_full(random_graph)
print('shortest paths-full kgfe, ER, 1000 nodes:', (time.time() - t0))


# TODO: try to implement a node subset?
import random
random_graph = nx.random_graphs.barabasi_albert_graph(10000, 3)
nodes_list = list(random_graph.nodes)
t0 = time.time()
for i in range(1000):
    nodes_set = random.sample(nodes_list, 20)
    nx_shortest_path_lengths = {}
    for n1 in nodes_set:
        nx_shortest_path_lengths[n1] = {}
        for n2 in nodes_set:
            dist = nx.shortest_path_length(random_graph, n1, n2)
            nx_shortest_path_lengths[n1][n2] = dist
print('shortest paths node subset nx, BA, 10000 nodes:', time.time() - t0)

from collections import defaultdict
t0 = time.time()
source_cache = defaultdict(lambda: {})
dest_cache = defaultdict(lambda: {})
for i in range(1000):
    nodes_set = random.sample(nodes_list, 20)
    paths_kgfe = kgfe.graph_utils.all_pairs_shortest_path_lengths(random_graph, nodes_set, source_cache, dest_cache)
print('shortest paths node subset kgfe, BA, 10000 nodes:', time.time() - t0)


