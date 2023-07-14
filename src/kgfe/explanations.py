# Graph explanations
# how?

from collections import Counter

import igraph as ig
from scipy.stats import hypergeom


def topic_pagerank(graph, topic_ids, topic_category=None, topic_weights=None,
        topic_id_prefix=None,
        alpha=0.7, max_iter=50, nstart=None):
    """
    Params:
        graph - a networkx graph
        topic_ids - a list of topics for personalization (random restart in pagerank) - use None for regular PR
        topic_category: 'Gene', 'Drug', 'SmallMolecule', 'Pathway'
        topic_weights - a dict of topic_id : weight TODO: not implemented

    Returns:
        dict of node id : pagerank score
    """
    # alpha is set to 0.7 based on https://academic.oup.com/bioinformatics/article/35/3/497/5055408
    # all random restarts go to the topic nodes.
    # TODO: igraph pagerank
    graph_ids = [graph.vs.find(name=t).index for t in topic_ids]
    pr_results = graph.personalized_pagerank(reset_vertices=graph_ids, damping=alpha)
    # postprocessing
    top_nodes = []
    ids_set = set(topic_ids)
    pr_results = Counter({v['name']: pr_results[v.index] for v in graph.vs})
    for node_id, score in pr_results.most_common():
        if node_id in ids_set:
            continue
        node = graph.vs.find(name=node_id).attributes().copy()
        node['score'] = score
        top_nodes.append(node)
    return pr_results, top_nodes

def steiner_tree(graph, ids, method='takahashi', **params):
    """
    A thin wrapper around steiner_tree.

    Methods to implement:
    - Kou
    - Takahashi
    - Wang
    - Mehlhorn
    """
    from . import steiner_tree
    # TODO: implement steiner tree method in igraph
    indices = [graph.vs.find(name=i).index for i in ids]
    if method == 'takahashi':
        tree = steiner_tree.takahashi_matsuyama_steiner_tree(graph, indices)
    return tree

def create_shortest_path_lengths_cached(graph):
    import functools
    @functools.cache
    def shortest_paths(n1, n2):
        return len(graph.get_shortest_paths(n1, n2)[0]) - 1
    return shortest_paths

def create_shortest_paths_cached(graph):
    import functools
    @functools.cache
    def shortest_paths(n1, n2):
        return graph.get_shortest_paths(n1, n2)
    return shortest_paths

def create_shortest_paths_cached_networkx(graph):
    import functools
    import networkx as nx
    @functools.cache
    def shortest_paths(n1, n2):
        return nx.shortest_path_length(graph, n1, n2)
    return shortest_paths

def graph_node_stats(graph, ids, target_nodes=None,
        shortest_paths_cached_function=None,
        metrics=None):
    """
    This works with igraph graphs. Currently, this only gets average pairwise distance (I wasn't really using any of the other statistics.)

    Args:
        - graph - a networkx graph
        - ids - a list of graph node ids
        - target_nodes - a list of nodes of interest that we want to find the distances to.
        - shortest_paths_cached_function - the output of create_shortest_paths_cached(graph)
        - metrics - default: pairwise only

    Gets some summary statistics for a set of nodes:
    - average pairwise distance
    - average distance from a node in the set to target node(s)
    """
    # average pairwise distance
    all_path_lengths = []
    if not shortest_paths_cached_function:
        shortest_paths_cached_function = create_shortest_path_lengths_cached(graph)
    for i, n1 in enumerate(ids[:-1]):
        paths = graph.get_shortest_paths(n1, ids[i+1:])
        all_path_lengths.extend(len(x) - 1 for x in paths)
    average_pairwise_distance = sum(all_path_lengths)/(len(all_path_lengths))
    # should jaccard similarity or clustering coefficient be included?
    if target_nodes is not None:
        target_node_distances = []
        for n1 in ids:
            for n2 in target_nodes:
                target_node_distances.append(shortest_paths_cached_function(n1, n2))
        average_target_distance = sum(target_node_distances)/len(target_node_distances)
        return {'average_pairwise_distance': average_pairwise_distance,
                'average_target_distance': average_target_distance,
                }
    return {'average_pairwise_distance': average_pairwise_distance,
            }


def graph_node_stats_networkx(graph, ids, target_nodes=None, shortest_paths_cached_function=None):
    """
    Args:
        - graph - a networkx graph
        - ids - a list of graph node ids
        - target_nodes - a list of nodes of interest that we want to find the distances to.
        - shortest_paths_cached_function - the output of create_shortest_paths_cached(graph)

    Gets some summary statistics for a set of nodes?
    - average pairwise distance
    - average clustering score
    - average jaccard score
    - average distance from a node in the set to target node(s)
    """
    import networkx as nx
    # cliquishness - clustering score
    clustering = nx.average_clustering(graph, ids)
    # average pairwise distance
    all_path_lengths = []
    all_pairs = []
    if not shortest_paths_cached_function:
        shortest_paths_cached_function = create_shortest_paths_cached_networkx(graph)
    for i, n1 in enumerate(ids[:-1]):
        for j in range(i+1, len(ids)):
            n2 = ids[j]
            all_pairs.append((n1, n2))
            all_path_lengths.append(shortest_paths_cached_function(n1, n2))
    average_pairwise_distance = sum(all_path_lengths)/len(all_path_lengths)
    # jaccard similarity coefficient of all pairs in ids - average fraction of neighbors shared among pairs of nodes in the set.
    jaccard_coefficients = [x[2] for x in nx.jaccard_coefficient(graph, all_pairs)]
    average_jaccard = sum(jaccard_coefficients)/len(jaccard_coefficients)
    if target_nodes is not None:
        target_node_distances = []
        for n1 in ids:
            for n2 in target_nodes:
                target_node_distances.append(shortest_paths_cached_function(n1, n2))
        average_target_distance = sum(target_node_distances)/len(target_node_distances)
        return {'average_pairwise_distance': average_pairwise_distance,
                'clustering': clustering,
                'average_jaccard': average_jaccard,
                'average_target_distance': average_target_distance,
                }
    return {'average_pairwise_distance': average_pairwise_distance,
            'clustering': clustering,
            'average_jaccard': average_jaccard,
            }

def null_graph_stats(graph, category, n_ids, n_samples=100, ids_subset=None):
    """
    This generates node set statistics for n_samples random sets of nodes of size n_ids, where all nodes are either belonging to category, or are part of the ids_subset (if provided)

    Args:
        graph: an igraph graph
        category: a category of nodes - e.g. 'Gene', 'Protein', etc.
        n_ids: the number of nodes to be selected in each sample.
        n_samples: the number of random samples.
        ids_subset: If this is not None, then category will not be used.

    Returns:
        all_stats - a list of dicts as returned by graph_node_stats
    """
    import random
    from .graph_info import nodes_in_category
    if ids_subset is None:
        ids_subset = [x['name'] for x in nodes_in_category(graph, category)]
    all_stats = []
    shortest_paths_cached_function = create_shortest_path_lengths_cached(graph)
    for _ in range(n_samples):
        id_sample = random.sample(ids_subset, n_ids)
        stats = graph_node_stats(graph, id_sample, shortest_paths_cached_function=shortest_paths_cached_function)
        all_stats.append(stats)
    return all_stats

def null_graph_stats_networkx(graph, category, n_ids, n_samples=100, ids_subset=None):
    """
    This generates node set statistics for n_samples random sets of nodes of size n_ids, where all nodes are either belonging to category, or are part of the ids_subset (if provided)

    Args:
        graph: a networkx graph
        category: a category of nodes - e.g. 'Gene', 'Protein', etc.
        n_ids: the number of nodes to be selected in each sample.
        n_samples: the number of random samples.
        ids_subset: If this is not None, then category will not be used.

    Returns:
        all_stats - a list of dicts as returned by graph_node_stats
    """
    import random
    from .graph_info import nodes_in_category_networkx
    if ids_subset is None:
        ids_subset = nodes_in_category_networkx(graph, category)
    all_stats = []
    shortest_paths_cached_function = create_shortest_paths_cached_networkx(graph)
    for _ in range(n_samples):
        id_sample = random.sample(ids_subset, n_ids)
        stats = graph_node_stats_networkx(graph, id_sample, shortest_paths_cached_function=shortest_paths_cached_function)
        all_stats.append(stats)
    return all_stats


def hypergeom_test(graph, query_ids, query_category, query_universe=None):
    """
    Hypergeometric test:
    N = number of nodes in the query category (OR the number of nodes in the query universe, if that's available),
    n = number of nodes in the query set,
    K = number of nodes connected to each given target node in the target category
    k = number of nodes connected to the target node that are in the query set

    Params:
        graph - a networkx graph
        query_ids - a list of IDs (NCBI Gene or PubChem), or a list of lists of ids
        query_category: 'Gene', 'Drug', 'SmallMolecule', 'Pathway'
        query_universe: either the universe of the query, or None if it's all nodes of the category in the graph.

    Returns:
        either a dict of hypergeometric p-values for node ids, or a list of dicts of hypergeometric p-values.
    """
    # 1. get all nodes of the query category in the graph
    # 2. get all nodes in the graph that are connected to nodes in the query set
    # 2. compute the overlaps and the hypergeometric score
    if query_universe is None:
        category_nodes = set([n.index for n in graph.vs if n['category'] == query_category])
    else:
        category_nodes = set([graph.vs.find(name=n).index for n in query_universe])
    def single_hypergeom(ids):
        indices = [graph.vs.find(name=i).index for i in ids]
        neighbors = set()
        for i in indices:
            neighbors.update(graph.neighbors(i))
        neighbor_vals = {}
        for n in neighbors:
            K = set([m for m in graph.neighbors(n) if m in category_nodes])
            k = K.intersection(indices)
            neighbor_vals[graph.vs[n]['name']] = (1 - hypergeom.cdf(len(k) - 1, len(category_nodes), len(K), len(query_ids)), {graph.vs[k1]['name'] for k1 in k})
        return neighbor_vals
    if isinstance(query_ids[0], list):
        all_results = []
        for ids in query_ids:
            all_results.append(single_hypergeom(ids))
        return all_results
    else:
        return single_hypergeom(query_ids)


def hypergeom_test_networkx(graph, query_ids, query_category, query_universe=None):
    """
    Hypergeometric test:
    N = number of nodes in the query category (OR the number of nodes in the query universe, if that's available),
    n = number of nodes in the query set,
    K = number of nodes connected to each given target node in the target category
    k = number of nodes connected to the target node that are in the query set

    Params:
        graph - a networkx graph
        query_ids - a list of IDs (NCBI Gene or PubChem), or a list of lists of ids
        query_category: 'Gene', 'Drug', 'SmallMolecule', 'Pathway'
        query_universe: either the universe of the query, or None if it's all nodes of the category in the graph. 

    Returns:
        either a dict of hypergeometric p-values for node ids, or a list of dicts of hypergeometric p-values.
    """
    # 1. get all nodes of the query category in the graph
    # 2. get all nodes in the graph that are connected to nodes in the query set
    # 2. compute the overlaps and the hypergeometric score
    if query_universe is None:
        category_nodes = set([n for n in graph.nodes if graph.nodes[n]['category'] == query_category])
    else:
        category_nodes = query_universe
    def single_hypergeom(ids):
        neighbors = set()
        for i in ids:
            neighbors.update([n for n in graph.neighbors(i)])
        neighbor_vals = {}
        for n in neighbors:
            K = set([m for m in graph.neighbors(n) if m in category_nodes])
            k = K.intersection(ids)
            neighbor_vals[n] = (1 - hypergeom.cdf(len(k) - 1, len(category_nodes), len(K), len(query_ids)), k)
        return neighbor_vals
    if isinstance(query_ids[0], list):
        all_results = []
        for ids in query_ids:
            all_results.append(single_hypergeom(ids))
        return all_results
    else:
        return single_hypergeom(query_ids)
