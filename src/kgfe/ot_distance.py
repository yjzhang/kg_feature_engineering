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

def precompute_shortest_paths_matrix(graph, nodes):
    "Returns a shortest paths matrix that returns the path length given a pair of indices in nodes."
    distances = graph.distances(nodes, nodes)
    return np.array(distances)


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

def _get_shortest_paths_table(nodes_1, nodes_2, distance_dict):
    "Generates a shortest paths table betwen the two sets of nodes, given a dict of distances."
    shortest_paths_table = np.zeros((len(nodes_1), len(nodes_2)))
    for i, n1 in enumerate(nodes_1):
        for j, n2 in enumerate(nodes_2):
            shortest_paths_table[i, j] = distance_dict[n1, n2]
    return shortest_paths_table

def ot_unbalanced_sinkhorn_distance_dict(graph, nodes_1, nodes_2, distance_dict, nodes_1_weights=None, nodes_2_weights=None, **sinkhorn_params):
    "Returns the unbalanced sinkhorn cost between two node lists, assuming that we have a dict that maps node pairs to distances."
    import ot
    if len(nodes_1) == 0 or len(nodes_2) == 0:
        return np.inf
    shortest_paths_table = _get_shortest_paths_table(nodes_1, nodes_2, distance_dict)
    a = np.ones(len(nodes_1))
    b = np.ones(len(nodes_2))
    if nodes_1_weights is not None:
        a = np.abs(nodes_1_weights)
    a = a/a.sum()
    if nodes_2_weights is not None:
        b = np.abs(nodes_2_weights)
    b = b/b.sum()
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
    cost = ot.unbalanced.sinkhorn_unbalanced2(a, b, shortest_paths_table, reg, reg_m, **sinkhorn_params)
    return cost[0]

def ot_balanced_sinkhorn_distance_dict(graph, nodes_1, nodes_2, distance_dict, nodes_1_weights=None, nodes_2_weights=None, **sinkhorn_params):
    "Returns the balanced sinkhorn cost between two node lists, assuming that we have a dict that maps node pairs to distances."
    import ot
    if len(nodes_1) == 0 or len(nodes_2) == 0:
        return np.inf
    shortest_paths_table = _get_shortest_paths_table(nodes_1, nodes_2, distance_dict)
    a = np.ones(len(nodes_1))
    b = np.ones(len(nodes_2))
    if 'reg' not in sinkhorn_params:
        reg = 1.0
    else:
        reg = sinkhorn_params['reg']
        del sinkhorn_params['reg']
    a = np.ones(len(nodes_1))
    b = np.ones(len(nodes_2))
    if nodes_1_weights is not None:
        a = np.abs(nodes_1_weights)
    a = a/a.sum()
    if nodes_2_weights is not None:
        b = np.abs(nodes_2_weights)
    b = b/b.sum()
    res = ot.sinkhorn(a, b, shortest_paths_table, reg, **sinkhorn_params)
    cost = np.sum(shortest_paths_table * res)
    return cost

def ot_balanced_exact_distance_dict(graph, nodes_1, nodes_2, distance_dict, nodes_1_weights=None, nodes_2_weights=None, **emd_params):
    "Uses the exact EMD solver to calculate the EMD between two node sets on a graph."
    import ot
    if len(nodes_1) == 0 or len(nodes_2) == 0:
        return np.inf
    shortest_paths_table = _get_shortest_paths_table(nodes_1, nodes_2, distance_dict)
    a = np.ones(len(nodes_1))
    b = np.ones(len(nodes_2))
    if nodes_1_weights is not None:
        a = np.abs(nodes_1_weights)
    a = a/a.sum()
    if nodes_2_weights is not None:
        b = np.abs(nodes_2_weights)
    b = b/b.sum()
    res = ot.emd(a, b, shortest_paths_table)
    cost = np.sum(shortest_paths_table * res)
    return cost

def ot_balanced_exact_distance(nodes_1_weights, nodes_2_weights, distances, **emd_params):
    "Uses the exact EMD solver to calculate the EMD between two node sets on a graph."
    import ot
    a = np.abs(nodes_1_weights)
    a = a/a.sum()
    b = np.abs(nodes_2_weights)
    b = b/b.sum()
    res = ot.emd(a, b, distances)
    cost = np.sum(distances * res)
    return cost

def ot_balanced_sinkhorn_distance(nodes_1_weights, nodes_2_weights, distances, **emd_params):
    "Uses the Sinkhorn solver to calculate the EMD between two node sets on a graph."
    import ot
    if 'reg' not in emd_params:
        reg = 1.0
    else:
        reg = emd_params['reg']
        del emd_params['reg']
    a = np.abs(nodes_1_weights)
    a = a/a.sum()
    b = np.abs(nodes_2_weights)
    b = b/b.sum()
    res = ot.sinkhorn(a, b, distances, reg, **emd_params)
    cost = np.sum(distances * res)
    return cost

def ot_sinkhorn_distance_matrix(graph, node_lists, distance_dict=None, node_weights=None, verbose=False, parallel=False, **sinkhorn_params):
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
            if node_weights is not None:
                distance_matrix[i, j] = ot_unbalanced_sinkhorn_distance_dict(graph, n1, n2, distance_dict, nodes_1_weights=node_weights[i], nodes_2_weights=node_weights[j], **sinkhorn_params)
            else:
                distance_matrix[i, j] = ot_unbalanced_sinkhorn_distance_dict(graph, n1, n2, distance_dict, **sinkhorn_params)
        if verbose and i%10 == 0:
            print('Nodes computed:', i)
    # convert matrix to full
    return distance_matrix + distance_matrix.T


def ot_balanced_sinkhorn_distance_matrix(graph, node_lists, distance_dict=None, node_weights=None, verbose=False, parallel=False, **sinkhorn_params):
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
            if node_weights is not None:
                distance_matrix[i, j] = ot_balanced_sinkhorn_distance_dict(graph, n1, n2, distance_dict, nodes_1_weights=node_weights[i], nodes_2_weights=node_weights[j], **sinkhorn_params)
            else:
                distance_matrix[i, j] = ot_balanced_sinkhorn_distance_dict(graph, n1, n2, distance_dict, **sinkhorn_params)
        if verbose and i%10 == 0:
            print('Nodes computed:', i)
    # convert matrix to full
    return distance_matrix + distance_matrix.T


def ot_balanced_emd_distance_matrix(graph, node_lists, distance_dict=None, node_weights=None, verbose=False, parallel=False, **emd_params):
    """
    Assumes that all node lists in node_lists have the same number of nodes.
    This creates a distance matrix between all pairs of node lists in node_lists, using the exact EMD algorithm.

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
            if node_weights is not None:
                distance_matrix[i, j] = ot_balanced_exact_distance_dict(graph, n1, n2, distance_dict, nodes_1_weights=node_weights[i], nodes_2_weights=node_weights[j], **emd_params)
            else:
                distance_matrix[i, j] = ot_balanced_exact_distance_dict(graph, n1, n2, distance_dict, **emd_params)
        if verbose and i%10 == 0:
            print('Nodes computed:', i)
    # convert matrix to full
    return distance_matrix + distance_matrix.T


def ot_data_distance_matrix(graph, data, nodes, distances=None, method='exact', verbose=False, parallel=False, **emd_params):
    """
    Calculates a distance matrix using OT on a dataset using a distance matrix.

    Args:
        graph: an igraph graph
        data: a n x k numpy array
        nodes: a list of node ids
        method: one of 'exact', 'sinkhorn'
    """
    if distances is None:
        distances = precompute_shortest_paths_matrix(graph, nodes)
        if verbose:
            print('Finished computing node-pair distances')
    distance_matrix = np.zeros((data.shape[0], data.shape[0]))
    for i in range(data.shape[0]):
        n1 = data[i, :]
        for j in range(i+1, data.shape[0]):
            n2 = data[j, :]
            if method == 'exact':
                distance_matrix[i, j] = ot_balanced_exact_distance(n1, n2, distances, **emd_params)
            elif method == 'sinkhorn':
                distance_matrix[i, j] = ot_balanced_sinkhorn_distance(n1, n2, distances, **emd_params)
        if verbose and i%10 == 0:
            print('Nodes computed:', i)
    # convert matrix to full
    return distance_matrix + distance_matrix.T


