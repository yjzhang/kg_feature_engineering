import networkx as nx
import pandas as pd
from .graph_info import get_names_to_ids

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


def generate_pairwise_features(df, graph, is_name=True, mode='avg'):
    """
    Generates a new dataframe with pairwise features from the graph...

    mode: 'avg', 'product', 'ratio'
    """
    genes = df.columns
    if is_name:
        names_to_ids = get_names_to_ids(graph)
        target_ids = [names_to_ids[n] for n in genes]
    else:
        target_ids = genes
    pairs = set()
    for gene_id in target_ids:
        neighbors = graph.neighbors(gene_id)
        for n in neighbors:
            if n in genes and (n, gene_id) not in pairs:
                pairs.add((gene_id, n))
        pass
    df_new = df.copy()
    for p in pairs:
        if mode == 'product':
            df_new[str(p[0]) + '-x-' + str(p[1])] = df[p[0]]*df[p[1]]
        elif mode == 'avg' or mode == 'average' or mode == 'mean':
            df_new[str(p[0]) + '-x-' + str(p[1])] = (df[p[0]]+df[p[1]])/2
        elif mode == 'harmonic':
            df_new[str(p[0]) + '-x-' + str(p[1])] = 2/(1/df[p[0]] + 1/df[p[1]])
        elif mode == 'ratio':
            df_new[str(p[0]) + '-/-' + str(p[1])] = df[p[0]]/df[p[1]]
            df_new[str(p[1]) + '-/-' + str(p[0])] = df[p[1]]/df[p[0]]
    return pairs, df_new


def generate_gene_set_features(df, graph, is_name=True,
        target_category='BiologicalProcess', mode='avg', included_gene_sets=None):
    """
    Generates a new dataframe with new features representing groups of genes (all genes connected to genes in the df)

    Params:
        df - pandas DataFrame containing the data
        graph -  networkx graph
        is_name - whether the inputs are gene names
        category - usually one of 'Gene',

    Returns a new dataframe
    """
    import numpy as np
    import scipy.stats
    genes = df.columns
    if is_name:
        names_to_ids = get_names_to_ids(graph)
        target_ids = [names_to_ids[n] for n in genes]
    else:
        target_ids = genes
    gene_sets_to_genes = {}
    genes_to_gene_sets = {}
    for gene_id in target_ids:
        neighbors = graph.neighbors(gene_id)
        genes_to_gene_sets[gene_id] = set()
        for n in neighbors:
            if graph.nodes[n]['category'] == target_category:
                if n in gene_sets_to_genes:
                    gene_sets_to_genes[n].add(gene_id)
                else:
                    gene_sets_to_genes[n] = set([gene_id])
                genes_to_gene_sets[gene_id].add(n)
    df_new = df.copy()
    for gs, genes in gene_sets_to_genes.items():
        if mode == 'product':
            df_new[gs] = np.product(df[genes].values, 1)
        elif mode == 'avg' or mode == 'average' or mode == 'mean':
            df_new[gs] = df[genes].mean(1)
        elif mode == 'harmonic':
            df_new[gs] = scipy.stats.hmean(df[genes].values)
    return gene_sets_to_genes, df_new


def get_repressor_features(graph, ids):
    """
    """
