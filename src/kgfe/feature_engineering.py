import networkx as nx
import pandas as pd

def get_feature_pairs(graph, ids, category='Gene'):
    """
    Returns a list of feature pairs within the ids with their specific interactions.
    """
    node_ids = {i for i in ids if i in graph.nodes}
    edges = []
    for i in node_ids:
        for n in graph.neighbors[i]:
            if n in node_ids:
                edges.append((i, n, graph.edges[i, n]))
    return edges


def generate_pairwise_features(df, graph, is_name=True, category='Gene', mode='avg'):
    """
    Generates a new dataframe with
    """

def generate_gene_set_features(df, graph, is_name=True, category='Gene',
        target_category='BiologicalProcess', mode='avg', included_gene_sets=None):
    """
    Generates a new dataframe with new features representing groups of genes (all genes connected to genes in the df)
    """
    # 1. get all BPs
    # 2. for all BPs in included_gene_sets, identify

def get_repressor_features(graph, ids):
    """
    """
