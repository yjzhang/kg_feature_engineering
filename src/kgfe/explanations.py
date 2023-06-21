# Graph explanations
# how?

from collections import Counter

import networkx as nx
from scipy.stats import hypergeom


def topic_pagerank(graph, topic_ids, topic_category=None, topic_weights=None,
        topic_id_prefix=None,
        alpha=0.85, max_iter=50, nstart=None):
    """
    Params:
        graph - a networkx graph
        topic_ids - a list of topics
        topic_category: 'Gene', 'Drug', 'SmallMolecule', 'Pathway'
        topic_weights - a dict of topic_id : weight

    Returns:
        dict of node id : pagerank score
    """
    if topic_weights is None:
        topic_weights = {i: 1 for i in topic_ids}
    if topic_weights is not None and topic_category is not None:
        topic_weights = {i: t for i, t in topic_weights.items() \
                if i in graph.nodes and graph.nodes[i]['category'] == topic_category}
    pr_results = nx.pagerank(graph, alpha=alpha,
            personalization=topic_weights, max_iter=max_iter, nstart=nstart)
    # postprocessing
    top_nodes = []
    ids_set = set(topic_ids)
    pr_results = Counter(pr_results)
    for node_id, score in pr_results.most_common():
        if node_id in ids_set:
            continue
        node = graph.nodes[node_id].copy()
        node['score'] = score
        top_nodes.append(node)
    return pr_results, top_nodes

def steiner_tree_subgraph(graph, ids, method='mehlhorn', **params):
    """
    Just a thin wrapper around the steiner tree method in networkx.
    Returns both the generated tree and a subgraph.
    """
    # TODO: get a connected subgraph with the given ids?
    steiner_tree = nx.approximation.steiner_tree(graph, ids, method=method, **params)
    subgraph = nx.subgraph(graph, steiner_tree)
    return steiner_tree, subgraph

def graph_node_stats(graph, ids):
    """
    Gets some summary statistics for a set of nodes?
    - average pairwise distance
    - average clustering score
    - average jaccard score
    """
    # cliquishness - clustering score
    clustering = nx.average_clustering(graph, ids)
    # average pairwise distance
    all_path_lengths = []
    all_pairs = []
    for i, n1 in enumerate(ids[:-1]):
        for j in range(i+1, len(ids)):
            n2 = ids[j]
            all_pairs.append((n1, n2))
            all_path_lengths.append(nx.shortest_path_length(graph, n1, n2))
    average_pairwise_distance = sum(all_path_lengths)/len(all_path_lengths)
    # jaccard similarity coefficient of all pairs in ids - average fraction of neighbors shared among pairs of nodes in the set.
    jaccard_coefficients = [x[2] for x in nx.jaccard_coefficient(graph, all_pairs)]
    average_jaccard = sum(jaccard_coefficients)/len(jaccard_coefficients)
    return {'average_pairwise_distance': average_pairwise_distance,
            'clustering': clustering,
            'average_jaccard': average_jaccard,
            }

def hypgergeom_test(graph, query_ids, query_category, query_universe=None):
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

