from collections import defaultdict
import math


def all_pairs_shortest_path_lengths_full(G, source_cache=None, dest_cache=None):
    # cache is a dict of start_node : {end_node : length}.
    if source_cache is None:
        source_cache = defaultdict(lambda: {})
    if dest_cache is None:
        dest_cache = defaultdict(lambda: {})
    # lengths is a dict of (n1, n2): length
    lengths = defaultdict(lambda: {})
    node_list = list(G.nodes)
    for i, n1 in enumerate(G.nodes):
        for i2 in range(i+1, len(G.nodes)):
            n2 = node_list[i2]
            lengths[n1][n2] = single_pair_shortest_path_length(G, n1, n2, source_cache, dest_cache)
    return lengths


def all_pairs_shortest_path_lengths(G, nodes, source_cache=None, dest_cache=None):
    """
    Does not do the same thing as the base networkx function - this returns all shortest path lengths for nodes within the set.
    """
    # cache is a dict of start_node : {end_node : length}.
    if source_cache is None:
        source_cache = defaultdict(lambda: {})
    if dest_cache is None:
        dest_cache = defaultdict(lambda: {})
    # lengths is a dict of (n1, n2): length
    lengths = {}
    for n1 in nodes:
        lengths[n1] = single_source_shortest_path_lengths(G, n1, nodes, source_cache, dest_cache)
    return lengths

def single_source_shortest_path_lengths(G, n1, targets, source_cache=None, dest_cache=None):
    # cache is a dict of start_node : {end_node : length}.
    targets = set(targets)
    if source_cache is None:
        source_cache = defaultdict(lambda: {})
    if dest_cache is None:
        dest_cache = defaultdict(lambda: {})
    # lengths is a dict of n2: length
    lengths = {}
    if n1 in targets:
        lengths[n1] = 0
    for n2 in targets:
        if n1 in source_cache and n2 in source_cache[n1]:
            lengths[n2] = source_cache[n1][n2]
        if n2 in source_cache and n1 in source_cache[n2]:
            lengths[n2] = source_cache[n2][n1]
    nextlevel = [n1]
    visited = set([n1])
    level = 0
    while nextlevel:
        level += 1
        thislevel = nextlevel
        nextlevel = []
        for v in thislevel:
            for w in G.neighbors(v):
                if w in visited:
                    continue
                visited.add(w)
                if w in targets:
                    lengths[w] = level
                    source_cache[n1][w] = level
                    if len(lengths) == len(targets):
                        return lengths
                source_cache[n1][w] = level
                nextlevel.append(w)
    for n2 in targets:
        if n2 not in lengths:
            lengths[n2] = math.inf
    return lengths


def single_pair_shortest_path_length(G, n1, n2, source_cache=None, dest_cache=None):
    if n1 == n2:
        return 0
    if source_cache is None:
        source_cache = defaultdict(lambda: {})
    if dest_cache is None:
        dest_cache = defaultdict(lambda: {})
    if n1 in source_cache and n2 in source_cache[n1]:
        return source_cache[n1][n2]
    if n2 in source_cache and n1 in source_cache[n2]:
        return source_cache[n2][n1]
    nextlevel = [n1]
    visited = set(nextlevel)
    level = 0
    while nextlevel:
        level += 1
        thislevel = nextlevel
        nextlevel = []
        for v in thislevel:
            for w in G.neighbors(v):
                if w in visited:
                    continue
                visited.add(w)
                if w == n2:
                    dist = level
                    source_cache[n1][n2] = dist
                    return dist
                source_cache[n1][w] = level
                nextlevel.append(w)
    source_cache[n1][n2] = math.inf
    return math.inf
