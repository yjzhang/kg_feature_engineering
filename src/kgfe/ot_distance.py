# TODO: optimal transport distance on graphs
# Given two node sets:
# 1. create a distance matrix using pairwise shortest paths between the two node sets
# 2. use some optimal transport solver
import numpy as np

def ot_table(graph, nodes_1, nodes_2):
    """
    Returns an optimal transport assignment table, and an optimal transport distance, between two node lists.
    """
    import ot
    shortest_paths_table = np.zeros((len(nodes_1), len(nodes_2)))
    for i, n in enumerate(nodes_1):
        path_lengths = graph.get_shortest_paths(n, nodes_2)
        for j, p in enumerate(path_lengths):
            shortest_paths_table[i, j] = len(p)
    # shortest_paths
    a = np.ones(len(nodes_1))
    b = np.ones(len(nodes_2))
    T = ot.emd(a, b, shortest_paths_table)
    cost = np.sum(T*shortest_paths_table)
    return T, cost
