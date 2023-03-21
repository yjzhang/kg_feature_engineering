# TODO: get info on available graphs
import os

import networkx as nx
import pandas as pd

PATH = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_PATH = os.path.join(PATH, 'processed_graphs')

def get_available_graphs():
    files = os.listdir(DATA_PATH)
    return files


def load_graph(filename):
    files = os.listdir(DATA_PATH)
    if filename not in files:
        # if filename not in files, try opening it as a path
        if not os.path.exists(filename):
            raise FileNotFoundError()
        f = open(filename)
    else:
        f = open(os.path.join(DATA_PATH, filename))
    if 'tsv' in filename:
        df = pd.read_csv(f, sep='\t')
    else:
        df = pd.read_csv(f)
    f.close()
    return df


def df_to_networkx(df):
    """
    Converts a panda dataframe to a networkx DiGraph, with node and edge attributes.
    """
    graph = nx.from_pandas_edgelist(df, source='subject_id', target='object_id',
            edge_attr=['predicate',
                       'Primary_Knowledge_Source',
                       'Knowledge_Source',
                       'publications'],
            create_using=nx.DiGraph)
    node_attributes = {}
    for i, row in df.iterrows():
        if row['subject_id'] not in node_attributes:
            node_attributes[row['subject_id']] = {
                    'id_prefix': row['subject_id_prefix'],
                    'name': row['subject_name'],
                    'category': row['subject_category']}
        if row['object_id'] not in node_attributes:
            node_attributes[row['object_id']] = {
                    'id_prefix': row['object_id_prefix'],
                    'name': row['object_name'],
                    'category': row['object_category']}
    nx.set_node_attributes(graph, node_attributes)
    return graph


def get_nodes_table(graph):
    """
    Returns a Pandas DataFrame of the nodes.
    """
    rows = []
    for n, attrs in graph.nodes.items():
        row = {'id': n}
        row.update(attrs)
        rows.append(row)
    return pd.DataFrame(rows)
