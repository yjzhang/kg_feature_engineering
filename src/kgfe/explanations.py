# Graph explanations
# how?
import networkx as nx


def topic_pagerank(graph, topic_ids, topic_category, topic_weights=None, topic_id_prefix=None, alpha=0.85):
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
    topic_weights = {i: t for i, t in topic_weights.items() \
            if i in graph.nodes and graph.nodes[i]['category'] == topic_category}
    pr_results = nx.pagerank(graph, topic_weights, alpha=alpha)
    return pr_results


def hypgergeom_test(graph, query_ids, query_category):
    """
    Params:
        graph - a networkx graph
        query_ids - a list of IDs (NCBI Gene or PubChem), or a list of lists of ids
        query_category: 'Gene', 'Drug', 'SmallMolecule', 'Pathway'
        query_weights - a dict of topic_id : weight

    Returns:
        either a dict of hypergeometric p-values for node ids, or a list of dicts of hypergeometric p-values.
    """
    # 1. get all nodes of the query category in the graph
    # 2. get all nodes in the graph that are connected to nodes in the query set
    # 2. compute the overlaps and the hypergeometric score
    def single_hypergeom():
        pass
    category_nodes = [n for n in graph.nodes if graph.nodes[n]['category'] == query_category]
    if isinstance(query_ids[0], list):
        pass
    else:
        connected_nodes = []
    pass

