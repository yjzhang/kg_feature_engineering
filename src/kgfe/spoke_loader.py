# Import spoke matrix from a neo4j csv dump.

import csv
import gzip
import json
import os

import numpy as np
from scipy import sparse, io


# TODO: multiple edges between two nodes?
def import_spoke_csv(csv_filename, edges_to_include=None, remove_unused_nodes=False, verbose=True, reindex_edges=True):
    """
    Args:
        csv_filename: name of csv file
        edges_to_include: set of edge types
        remove_unused_nodes: True if nodes with no in- or out-edges are to be removed.
        reindex_edges: whether or not to use indices or original IDs in the edge list.

    Returns:
        nodes: list of (_id, _name, _labels_id) where _labels_id corresponds to a key in node_types
        edges: dict of (node1, node2): _type_id where node1 and node2 index into nodes, and _type_id corresponds to a key in edge_types
        node_types: dict of int: str (_labels)
        edge_types: dict of int: str (_type)
    """
    nodes = []
    n_nodes = 0
    # mapping of _id to index in nodes
    node_index = {}
    # node_types is a map of string (
    node_types = {}
    edges = {}
    # edge_types is a map of string (_type) to node
    edge_types = {}
    # sets of nodes that have in-edges or out-edges (to use when deciding whether to remove nodes)
    node_has_edge = set()
    csv.field_size_limit(99999999)
    if csv_filename.endswith('.gz'):
        # handle gzip
        f = gzip.open(csv_filename, 'rt')
    else:
        f = open(csv_filename)
    dr = csv.DictReader(f, dialect='unix')
    for i, row in enumerate(dr):
        if verbose and i % 100000 == 0:
            print(i, 'nodes: ', len(node_index), 'edges: ', len(edges))
        # if this is a node
        if row['_id']:
            print(row['license'])
            if row['name']:
                row_name = row['name']
                print(row_name)
            else:
                row_name = row['pref_name']
            if row['_labels'] in node_types:
                nodes.append((int(row['_id']), row_name, node_types[row['_labels']]))
            else:
                nodes.append((int(row['_id']), row_name, len(node_types) + 1))
                node_types[row['_labels']] = len(node_types) + 1
            node_index[int(row['_id'])] = n_nodes 
            n_nodes += 1
        # if this row is an edge
        else:
            edge_type = row['_type']
            if edges_to_include is None or edge_type in edges_to_include:
                node1 = int(row['_start'])
                node2 = int(row['_end'])
                node_has_edge.add(node1)
                node_has_edge.add(node2)
                if edge_type in edge_types:
                    edges[(node1, node2)] = edge_types[edge_type]
                else:
                    edges[(node1, node2)] = len(edge_types) + 1
                    edge_types[row['_type']] = len(edge_types) + 1
    if remove_unused_nodes:
        # remove all nodes that don't have edges
        to_remove = set(node_index.keys()).difference(node_has_edge)
        nodes = [n for n in nodes if n[0] not in to_remove]
        # rebuild node_index
        node_index = {n[0]: i for i, n in enumerate(nodes)}
    # convert edge indices
    if reindex_edges:
        new_edges = {}
        for k, e in edges.items():
            node1, node2 = k
            node1 = node_index[node1]
            node2 = node_index[node2]
            new_edges[(node1, node2)] = e
        edges = new_edges
    node_types = {v: k for k, v in node_types.items()}
    edge_types = {v: k for k, v in edge_types.items()}
    return nodes, edges, node_types, edge_types

# TODO: multiple edges between two nodes?
def import_spoke_jsonl(filename, edges_to_include=None, remove_unused_nodes=True, use_edge_types=True, use_node_types=True, verbose=True, reindex_edges=True):
    """
    Imports a jsonl file.
    Args:
        filename: name of jsonl file
        edges_to_include: set of edge types
        remove_unused_nodes: True if nodes with no in- or out-edges are to be removed.
        reindex_edges: whether or not to use indices or original IDs in the edge list.

    Returns:
        nodes: list of (_id, _name, _labels_id) where _labels_id corresponds to a key in node_types
        edges: dict of (node1, node2): _type_id where node1 and node2 index into nodes, and _type_id corresponds to a key in edge_types
        node_types: dict of int: str (_labels)
        edge_types: dict of int: str (_type)
    """
    nodes = []
    n_nodes = 0
    # mapping of _id to index in nodes
    node_index = {}
    # node_types is a map of string (
    node_types = {}
    edges = {}
    # edge_types is a map of string (_type) to node
    edge_types = {}
    # sets of nodes that have in-edges or out-edges (to use when deciding whether to remove nodes)
    node_has_edge = set()
    if filename.endswith('.gz'):
        # handle gzip
        f = gzip.open(filename, 'rt')
    else:
        f = open(filename)
    line = f.readline()
    i = 0
    while line:
        row = json.loads(line)
        if verbose and i % 100000 == 0:
            print(i, 'nodes: ', len(node_index), 'edges: ', len(edges))
        # if this is a node
        if row['type'] == 'node':
            row_name = ''
            row_identifier = ''
            row_source = ''
            if 'name' in row['properties'] and row['properties']['name'] != '':
                row_name = row['properties']['name']
            elif 'pref_name' in row['properties'] and row['properties']['pref_name'] != '':
                row_name = row['properties']['pref_name']
            elif 'identifier' in row['properties'] and row['properties']['identifier'] != '':
                row_name = row['properties']['identifier']
            elif 'id' in row['properties'] and row['properties']['id']:
                row_name = row['properties']['id']
            else:
                row_name = ''
            if 'identifier' in row['properties'] and row['properties']['identifier'] != '':
                row_identifier = row['properties']['identifier']
            if 'source' in row['properties'] and row['properties']['source'] != '':
                row_source = row['properties']['source']
            row_label = row['labels'][0]
            if use_node_types:
                if row_label in node_types:
                    nodes.append((int(row['id']), row_name, node_types[row_label], row_identifier, row_source))
                else:
                    nodes.append((int(row['id']), row_name, len(node_types) + 1, row_identifier, row_source))
                    node_types[row_label] = len(node_types) + 1
            else:
                nodes.append((int(row['id']), row_name, True, row_identifier, row_source))
            node_index[int(row['id'])] = n_nodes 
            n_nodes += 1
        # if this row is an edge
        else:
            edge_type = row['label']
            if edges_to_include is None or edge_type in edges_to_include:
                node1 = int(row['start']['id'])
                node2 = int(row['end']['id'])
                node_has_edge.add(node1)
                node_has_edge.add(node2)
                if use_edge_types:
                    if edge_type in edge_types:
                        edges[(node1, node2)] = edge_types[edge_type]
                    else:
                        edges[(node1, node2)] = len(edge_types) + 1
                        edge_types[edge_type] = len(edge_types) + 1
                else:
                    edges[(node1, node2)] = True
        line = f.readline()
        i += 1
    if remove_unused_nodes:
        # remove all nodes that don't have edges
        to_remove = set(node_index.keys()).difference(node_has_edge)
        nodes = [n for n in nodes if n[0] not in to_remove]
        # rebuild node_index
        node_index = {n[0]: i for i, n in enumerate(nodes)}
    # convert edge indices
    if reindex_edges:
        new_edges = {}
        for k, e in edges.items():
            node1, node2 = k
            node1 = node_index[node1]
            node2 = node_index[node2]
            new_edges[(node1, node2)] = e
        edges = new_edges
    node_types = {v: k for k, v in node_types.items()}
    edge_types = {v: k for k, v in edge_types.items()}
    return nodes, edges, node_types, edge_types


def import_ckg_jsonl(filename, edges_to_include=None, remove_unused_nodes=False, use_edge_types=True, use_node_types=True, n_edges=300000000, n_nodes=20000000, verbose=True):
    """
    Imports a jsonl file.
    This tries to be less memory-intensive than the other import procedure.
    Args:
        filename: name of jsonl file
        edges_to_include: set of edge types
        remove_unused_nodes: True if nodes with no in- or out-edges are to be removed.
        n_edges: An upper bound on the number of edges (does not have to be exact, but should be greater than the actual number of edges)

    Returns:
        nodes: list of (_id, _name, _labels_id) where _labels_id corresponds to a key in node_types
        edges: COO array
        node_types: dict of int: str (_labels)
        edge_types: dict of int: str (_type)
    """
    nodes = []
    n_nodes = 0
    # mapping of _id to index in nodes
    node_index = {}
    # node_types is a map of string (
    node_types = {}
    edges_start = np.zeros(n_edges, dtype=np.int64)
    edges_end = np.zeros(n_edges, dtype=np.int64)
    edges_values = np.zeros(n_edges, dtype=np.uint8)
    # edge_types is a map of string (_type) to node
    edge_types = {}
    # sets of nodes that have in-edges or out-edges (to use when deciding whether to remove nodes)
    if filename.endswith('.gz'):
        # handle gzip
        f = gzip.open(filename, 'rt')
    else:
        f = open(filename)
    line = f.readline()
    i = 0
    # ne is number of current edges
    ne = 0
    while line:
        row = json.loads(line)
        if verbose and i % 100000 == 0:
            print(i, 'nodes: ', len(node_index), 'edges: ', ne)
        # if this is a node
        if row['type'] == 'node':
            if 'name' in row['properties'] and row['properties']['name'] != '':
                row_name = row['properties']['name']
            elif 'pref_name' in row['properties'] and row['properties']['pref_name'] != '':
                row_name = row['properties']['pref_name']
            elif 'identifier' in row['properties'] and row['properties']['identifier'] != '':
                row_name = row['properties']['identifier']
            elif 'id' in row['properties'] and row['properties']['id']:
                row_name = row['properties']['id']
            else:
                row_name = ''
            row_label = row['labels'][0]
            if use_node_types:
                if row_label in node_types:
                    nodes.append((int(row['id']), row_name, node_types[row_label]))
                else:
                    nodes.append((int(row['id']), row_name, len(node_types) + 1))
                    node_types[row_label] = len(node_types) + 1
            else:
                nodes.append((int(row['id']), row_name, True))
            node_index[int(row['id'])] = n_nodes 
            n_nodes += 1
        # if this row is an edge
        # in neo4j exports, edges always come after nodes.
        # assumption: there are less than 255 edge types
        else:
            edge_type = row['label']
            if edges_to_include is None or edge_type in edges_to_include:
                node1 = node_index[int(row['start']['id'])]
                node2 = node_index[int(row['end']['id'])]
                edges_start[ne] = node1
                edges_end[ne] = node2
                if use_edge_types:
                    if edge_type in edge_types:
                        edges_values[ne] = edge_types[edge_type]
                    else:
                        edges_values[ne] = len(edge_types) + 1
                        edge_types[edge_type] = len(edge_types) + 1
                else:
                    edges_values[ne] = 1
                ne += 1
        line = f.readline()
        i += 1
    edges_values = edges_values[:ne]
    edges_start = edges_start[:ne]
    edges_end = edges_end[:ne]
    edges = sparse.coo_array((edges_values, (edges_start, edges_end)), shape=(len(nodes), len(nodes)))
    node_types = {v: k for k, v in node_types.items()}
    edge_types = {v: k for k, v in edge_types.items()}
    return nodes, edges, node_types, edge_types





def to_sparse(nodes, edges):
    """
    Returns a DOK matrix from the edges...
    """
    n_nodes = len(nodes)
    edge_matrix = sparse.dok_array((n_nodes, n_nodes), dtype=int)
    for k, v in sorted(edges.items()):
        n1, n2 = k
        edge_matrix[n1, n2] = v
    return edge_matrix


def load_spoke(filename='spoke.csv', edges_to_include=None, remove_unused_nodes=False, mtx_filename='spoke.mtx', **kwargs):
    if filename.endswith('.csv') or filename.endswith('.csv.gz'):
        nodes, edges, node_types, edge_types = import_spoke_csv(filename, edges_to_include, remove_unused_nodes, **kwargs)
    elif filename.endswith('.json') or filename.endswith('.json.gz') or filename.endswith('.jsonl') or filename.endswith('.jsonl.gz'):
        nodes, edges, node_types, edge_types = import_spoke_jsonl(filename, edges_to_include, remove_unused_nodes, **kwargs)
    if not os.path.exists(mtx_filename):
        edge_matrix = to_sparse(nodes, edges)
        io.mmwrite(mtx_filename, edge_matrix)
    else:
        edge_matrix = io.mmread(mtx_filename)
    return nodes, edges, node_types, edge_types, edge_matrix


def load_spoke_networkx(filename='spoke.csv', edges_to_include=None, remove_unused_nodes=True, directed=False, **kwargs):
    import networkx as nx
    if filename.endswith('.csv') or filename.endswith('.csv.gz'):
        nodes, edges, node_types, edge_types = import_spoke_csv(filename, edges_to_include, remove_unused_nodes, reindex_edges=False, **kwargs)
    elif filename.endswith('.json') or filename.endswith('.json.gz') or filename.endswith('.jsonl') or filename.endswith('.jsonl.gz'):
        nodes, edges, node_types, edge_types = import_spoke_jsonl(filename, edges_to_include, remove_unused_nodes, reindex_edges=False, **kwargs)
    edge_list = edges.keys()
    if directed:
        graph = nx.from_edgelist(edge_list, nx.DiGraph)
    else:
        graph = nx.from_edgelist(edge_list)
    # set node attributes
    node_attributes = {}
    # TODO: get all IDs, not just names
    for n in nodes:
        node_attributes[n[0]] = {
                'name':  n[1],
                'category': node_types[n[2]],
                'identifier': n[3],
                'source': n[4],
        }
    nx.set_node_attributes(graph, node_attributes)
    # set edge attributes
    edge_attributes = {k: {'type': edge_types[v]} for k, v in edges.items()}
    nx.set_edge_attributes(graph, edge_attributes)
    return graph

def load_spoke_igraph(filename='spoke.csv', edges_to_include=None, remove_unused_nodes=True, directed=False, verbose=False, low_memory=False, **kwargs):
    """
    Imports the spoke file as an igraph. The file can be a csv or json/jsonl export from neo4j, and it can be gzipped. The spoke IDs are converted to strings because igraph is very slow if the ids are ints.
    """
    import igraph as ig
    if filename.endswith('.csv') or filename.endswith('.csv.gz'):
        nodes, edges, node_types, edge_types = import_spoke_csv(filename, edges_to_include, remove_unused_nodes, reindex_edges=False, verbose=verbose, **kwargs)
    elif filename.endswith('.json') or filename.endswith('.json.gz') or filename.endswith('.jsonl') or filename.endswith('.jsonl.gz'):
        nodes, edges, node_types, edge_types = import_spoke_jsonl(filename, edges_to_include, remove_unused_nodes, reindex_edges=False, verbose=verbose, **kwargs)
    else:
        raise ValueError('File has to be a csv, csv.gz, json, json.gz file')
    if verbose:
        print('Done loading data, creating edge list')
    # use igraph.graph.DictList
    if low_memory:
        edge_list = ({'source': str(v[0]), 'target': str(v[1])} for v in edges.keys())
    else:
        edge_list = ({'source': str(v[0]), 'target': str(v[1]), 'type': edge_types[e]} for v, e in edges.items())
    if verbose:
        print('creating node list')
    # set node attributes
    # convert the node id to a string, bc
    if low_memory:
        node_list = ({
             'name': str(n[0]),
             'feature_name':  n[1],
             'category': node_types[n[2]],
             'identifier': n[3],
            } for n in nodes)
    else:
        node_list = ({
                'name': str(n[0]),
                'feature_name':  n[1],
                'category': node_types[n[2]],
                'identifier': n[3],
                'source': n[4],
        } for n in nodes)
    if verbose:
        print('calling igraph.Graph.DictList')
    graph = ig.Graph.DictList(node_list, edge_list, directed=directed)
    return graph




def symmetrize_matrix(matrix):
    """
    Symmetrizes an adjacency matrix.

    Warning: this completely destroys any meaning applied to node values. Nonzero = edge exists, zero = edge doesn't exist.
    """
    lower_triangle = sparse.tril(matrix)
    upper_triangle = sparse.triu(matrix)
    return lower_triangle + lower_triangle.T + upper_triangle + upper_triangle.T

