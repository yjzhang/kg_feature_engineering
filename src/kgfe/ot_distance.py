# optimal transport distance on graphs
# Given two node sets:
# 1. create a distance matrix using pairwise shortest paths between the two node sets
# 2. use some optimal transport solver
import numpy as np

def ot_table(graph, nodes_1, nodes_2, **emd_params):
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
    T = ot.emd(a/a.sum(), b/b.sum(), shortest_paths_table, **emd_params)
    cost = np.sum(T*shortest_paths_table)
    return T, cost

def ot_cost(graph, nodes_1, nodes_2, **emd_params):
    "Only returns the OT cost (reduces compute time by removing duplicate nodes, so this will return a different result from ot_table.)"
    import ot
    s1 = set(nodes_1)
    s2 = set(nodes_2)
    nodes_1 = list(s1.difference(s2))
    nodes_2 = list(s2.difference(s1))
    shortest_paths_table = np.array(graph.distances(nodes_1, nodes_2))
    # shortest_paths
    a = np.ones(len(nodes_1))
    b = np.ones(len(nodes_2))
    cost = ot.emd2(a/a.sum(), b/b.sum(), shortest_paths_table, **emd_params)
    return cost

# for some applications (i.e. when computing many distances over a small subset of nodes), pre-computing a shortest paths table for a set of nodes can improve performance.
def precompute_shortest_paths(graph, nodes):
    "Returns a shortest paths dict that returns the path length given a pair of node ids."
    distances = graph.distances(nodes, nodes)
    distance_dict = {}
    for i, n1 in enumerate(nodes):
        for j, n2 in enumerate(nodes):
            distance_dict[n1, n2] = distances[i][j]
    return distance_dict

def ot_sinkhorn_unbalanced_cost(graph, nodes_1, nodes_2, **sinkhorn_params):
    "Returns the unbalanced sinkhorn cost between two node lists."
    # remove duplicate nodes between the two groups (necessary for graph.distances)
    import ot
    s1 = set(nodes_1)
    s2 = set(nodes_2)
    nodes_1 = list(s1.difference(s2))
    nodes_2 = list(s2.difference(s1))
    # compute shortest_paths
    shortest_paths_table = np.array(graph.distances(nodes_1, nodes_2))
    a = np.ones(len(nodes_1))
    b = np.ones(len(nodes_2))
    # regularization defaults
    if 'reg' not in sinkhorn_params:
        reg = 1
    else:
        reg = sinkhorn_params['reg']
        del sinkhorn_params['reg']
    if 'reg_m' not in sinkhorn_params:
        reg_m = 1
    else:
        reg_m = sinkhorn_params['reg_m']
        del sinkhorn_params['reg_m']
    cost = ot.unbalanced.sinkhorn_unbalanced2(a/a.sum(), b/b.sum(), shortest_paths_table, reg, reg_m, **sinkhorn_params)[0]
    return cost

def ot_sinkhorn_unbalanced_cost_distance_dict(graph, nodes_1, nodes_2, distance_dict, **sinkhorn_params):
    "Returns the unbalanced sinkhorn cost between two node lists, assuming that we have a dict that maps node pairs to distances."
    import ot
    if len(nodes_1) == 0 or len(nodes_2) == 0:
        return np.inf
    shortest_paths_table = np.zeros((len(nodes_1), len(nodes_2)))
    for i, n1 in enumerate(nodes_1):
        for j, n2 in enumerate(nodes_2):
            shortest_paths_table[i, j] = distance_dict[n1, n2]
    a = np.ones(len(nodes_1))
    b = np.ones(len(nodes_2))
    if 'reg' not in sinkhorn_params:
        reg = 1
    else:
        reg = sinkhorn_params['reg']
        del sinkhorn_params['reg']
    if 'reg_m' not in sinkhorn_params:
        reg_m = 1
    else:
        reg_m = sinkhorn_params['reg_m']
        del sinkhorn_params['reg_m']
    cost = ot.unbalanced.sinkhorn_unbalanced2(a/a.sum(), b/b.sum(), shortest_paths_table, reg, reg_m, **sinkhorn_params)
    return cost[0]

def ot_sinkhorn_distance_dict(graph, nodes_1, nodes_2, distance_dict, **sinkhorn_params):
    "Returns the balanced sinkhorn cost between two node lists, assuming that we have a dict that maps node pairs to distances."
    import ot
    if len(nodes_1) == 0 or len(nodes_2) == 0:
        return np.inf
    shortest_paths_table = np.zeros((len(nodes_1), len(nodes_2)))
    for i, n1 in enumerate(nodes_1):
        for j, n2 in enumerate(nodes_2):
            shortest_paths_table[i, j] = distance_dict[n1, n2]
    a = np.ones(len(nodes_1))
    b = np.ones(len(nodes_2))
    if 'reg' not in sinkhorn_params:
        reg = 0.1
    else:
        reg = sinkhorn_params['reg']
        del sinkhorn_params['reg']
    cost = ot.sinkhorn2(a/a.sum(), b/b.sum(), shortest_paths_table, reg, **sinkhorn_params)
    return cost[0]

def ot_sinkhorn_distance_matrix(graph, node_lists, distance_dict=None, verbose=False, parallel=False, **sinkhorn_params):
    """
    This creates a distance matrix between all pairs of node lists in node_lists, using unbalanced sinkhorn.

    Args:
        graph: an igraph graph
        node_lists: a list of lists of nodes
        distance_dict: optional, a dict mapping pairs (tuples) of node ids to shortest-paths distances. If this is not provided, one will be calculated. 
        verbose: boolean, whether or not to print progress
        parallel: not implemented
    """
    if distance_dict is None:
        node_set = set()
        for n in node_lists:
            node_set.update(n)
        distance_dict = precompute_shortest_paths(graph, list(node_set))
        if verbose:
            print('Finished computing node-pair distances')
    distance_matrix = np.zeros((len(node_lists), len(node_lists)))
    for i, n1 in enumerate(node_lists):
        for j in range(i+1, len(node_lists)):
            n2 = node_lists[j]
            distance_matrix[i, j] = ot_sinkhorn_unbalanced_cost_distance_dict(graph, n1, n2, distance_dict, **sinkhorn_params)
        if verbose and i%10 == 0:
            print('Nodes computed:', i)
    # convert matrix to full
    return distance_matrix + distance_matrix.T


def ot_balanced_sinkhorn_distance_matrix(graph, node_lists, distance_dict=None, verbose=False, parallel=False, **sinkhorn_params):
    """
    Assumes that all node lists in node_lists have the same number of nodes.
    This creates a distance matrix between all pairs of node lists in node_lists, using balanced sinkhorn.

    Args:
        graph: an igraph graph
        node_lists: a list of lists of nodes
        distance_dict: optional, a dict mapping pairs (tuples) of node ids to shortest-paths distances. If this is not provided, one will be calculated. 
        verbose: boolean, whether or not to print progress
        parallel: not implemented
    """
    if distance_dict is None:
        node_set = set()
        for n in node_lists:
            node_set.update(n)
        distance_dict = precompute_shortest_paths(graph, list(node_set))
        if verbose:
            print('Finished computing node-pair distances')
    distance_matrix = np.zeros((len(node_lists), len(node_lists)))
    for i, n1 in enumerate(node_lists):
        for j in range(i+1, len(node_lists)):
            n2 = node_lists[j]
            distance_matrix[i, j] = ot_sinkhorn_distance_dict(graph, n1, n2, distance_dict, **sinkhorn_params)
        if verbose and i%10 == 0:
            print('Nodes computed:', i)
    # convert matrix to full
    return distance_matrix + distance_matrix.T


