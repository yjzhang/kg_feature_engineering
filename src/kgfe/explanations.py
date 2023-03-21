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
    """
    if topic_weights is None:
        topic_weights = {i: 1 for i in topic_weights}
    topic_weights = {i: t for i, t in topic_weights.items() \
            if i in graph.nodes and graph.nodes[i]['category'] == topic_category}
    pr_results = nx.pagerank(graph, topic_weights, alpha=alpha)
    return pr_results


def hypgergeom_test(graph, topic_ids, topic_category):
    """
    """
    pass
