# TODO: get info on available graphs
import os
import zipfile

import networkx as nx
import pandas as pd

PATH = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_PATH = os.path.join(PATH, 'processed_graphs')
MSIGDB_PATH = os.path.join(PATH, 'raw_graphs/msigdb_v2023.1.Hs_json_files_to_download_locally.zip')

def get_available_graphs():
    files = os.listdir(DATA_PATH)
    return files

def get_available_msigdb():
    f = zipfile.ZipFile(MSIGDB_PATH)
    files = f.namelist()
    f.close()
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


def _load_msigdb(data, object_category='BiologicalProcess', object_id_prefix='MSigDB',
        predicate='participates_in'):
    """
    Returns a pandas edge table and a graph
    """
    import gene_names
    new_entries = []
    for k, v in data.items():
        gene_set_name = k
        source = v['exactSource']
        ref = v['pmid']
        genes = v['geneSymbols']
        try:
            gene_ids = gene_names.get_ids(genes)
        except:
            gene_ids = []
            for g in genes:
                try:
                    gene_id = gene_names.get_ids([g])
                    gene_ids.append(gene_id)
                except:
                    continue
        for gene_id, symbol in zip(gene_ids, genes):
            entry = {}
            entry['subject_category'] = 'Gene'
            entry['subject_id_prefix'] = 'NCBIGene'
            entry['subject_id'] = gene_id
            entry['subject_name'] = symbol
            entry['predicate'] = predicate
            entry['object_category'] = object_category
            entry['object_id_prefix'] = object_id_prefix
            entry['object_id'] = gene_set_name
            entry['object_name'] = gene_set_name
            entry['Primary_Knowledge_Source'] = source
            entry['Knowledge_Source'] = v['collection']
            entry['publications'] = ref
            new_entries.append(entry)
    return pd.DataFrame(new_entries)


def load_msigdb(name):
    """Load an msigdb graph as a dataframe. Input name is one of the filenames returned by get_available_msigdb."""
    import json
    f = zipfile.ZipFile(MSIGDB_PATH)
    m = f.getinfo(name)
    json_file = f.open(m)
    data = json.load(json_file)
    df = _load_msigdb(data)
    f.close()
    return df


def df_to_networkx(df, directed=False):
    """
    Converts a panda dataframe to a networkx Graph (or DiGraph), with node and edge attributes.
    """
    create_using = nx.Graph
    if directed:
        create_using = nx.DiGraph
    df['subject_id_full'] = df['subject_id_prefix'] + '::' + df['subject_id'].astype(str)
    df['object_id_full'] = df['object_id_prefix'] + '::' + df['object_id'].astype(str)
    graph = nx.from_pandas_edgelist(df, source='subject_id_full', target='object_id_full',
            edge_attr=['predicate',
                       'Primary_Knowledge_Source',
                       'Knowledge_Source',
                       'publications'],
            create_using=create_using)
    node_attributes = {}
    for row in df.rows():
        if row['subject_id'] not in node_attributes:
            node_attributes[row['subject_id_full']] = {
                    'id': row['subject_id'],
                    'id_prefix': row['subject_id_prefix'],
                    'name': row['subject_name'],
                    'category': row['subject_category']}
        if row['object_id'] not in node_attributes:
            node_attributes[row['object_id_full']] = {
                    'id': row['object_id'],
                    'id_prefix': row['object_id_prefix'],
                    'name': row['object_name'],
                    'category': row['object_category']}
    nx.set_node_attributes(graph, node_attributes)
    return graph

# TODO: create a new class for graphs? subclass of networkx graphs?


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

def get_names_to_ids(graph):
    """Returns a dict mapping node names to IDs (ignoring prefixes and categories so on)"""
    names_to_ids = {}
    for n, attrs in graph.nodes.items():
        names_to_ids[attrs['name']] = n
    return names_to_ids

def get_spoke_categories(graph):
    return set(attrs['category'] for attrs in graph.nodes.values())

def get_spoke_sources(graph):
    return set(attrs['source'] for attrs in graph.nodes.values())

def spoke_identifiers_to_ids(graph, category, source=None):
    """
    Returns a mapping from SPOKE identifiers to IDs.

    category: 'Protein', 'Gene', 'Compound', 'Disease', etc
    source: 'KEGG', ...
    """
    identifiers_to_ids = {}
    for n, attrs in graph.nodes.items():
        if 'category' in attrs and attrs['category'] == category:
            if source is None or ('source' in attrs and attrs['source'] == source):
                identifiers_to_ids[attrs['identifier']] = n
    return identifiers_to_ids

def get_category_ids_to_nodes(graph, category):
    """
    Returns a dict that maps from identifiers in the specified category to graph node IDs.
    """
    # TODO
    identifiers_to_ids = {}
    return identifiers_to_ids

def random_nodes_in_category(graph, category, n_nodes):
    """
    Returns a list of random spoke ids in the given category.
    """
    import random
    nodes_in_category = []
    for n, attrs in graph.nodes.items():
        if 'category' in attrs and attrs['category'] == category:
            nodes_in_category.append((n, attrs['identifier']))
    return random.sample(nodes_in_category, n_nodes)

# TODO: random nodes with similar degree distributions? investigative bias - constrain null model to be similar to the problem. We could select random nodes among the nodes that are in the general set...
# subgraph
# randomize the graph... randomly shuffle the nodes/edges?
# shuffle the edges in the network, look at the set again and again. preserve the degree...

# TODO: estimate from whole graph
def randomize_graph(G):
    """
    Randomizes the given graph using edge swapping while preserving the degree distribution.

    Args:
    G (networkx.Graph): Input graph.

    Returns:
    networkx.Graph: Randomized graph.
    """
    import random
    G_rand = G.copy()
    edges = list(G_rand.edges)
    nodes = list(G_rand.nodes)
    n = len(edges)

    for _ in range(n):
        # Pick two edges and swap
        edge1 = random.choice(edges)
        edge2 = random.choice(edges)
