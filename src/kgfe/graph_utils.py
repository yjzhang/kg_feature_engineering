import networkx as nx


def all_pairs_shortest_path_lengths(G, nodes, source_cache=None, dest_cache=None):
    """
    Does not do the same thing as the base networkx function - this returns all shortest path lengths for nodes within the set.
    """
    # cache is a dict of start_node : {end_node : length}.
    if source_cache is None:
        source_cache = {}
    if dest_cache is None:
        dest_cache = {}
    # lengths is a dict of (n1, n2): length
    lengths = {}
    for i, n1 in enumerate(nodes):
        for i2 in range(i+1, len(nodes)):
            n2 = nodes[i2]
            lengths[(n1, n2)] = single_pair_shortest_path_length(G, n1, n2, source_cache, dest_cache)
    return lengths

def single_pair_shortest_path_length(G, n1, n2, source_cache=None, dest_cache=None):
    if n1 == n2:
        return 0
    if source_cache is None:
        source_cache = {}
    if dest_cache is None:
        dest_cache = {}
    if n1 in source_cache and n2 in source_cache[n1]:
        return source_cache[n1][n2]
    if n2 in dest_cache and n1 in dest_cache[n2]:
        return dest_cache[n2][n1]
    if n1 not in source_cache:
        source_cache[n1] = {}
    if n2 not in dest_cache:
        dest_cache[n2] = {}
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
                if w == n2:
                    dist = level
                    dest_cache[n2][n1] = dist
                    source_cache[n1][n2] = dist
                    return dist
                if w in dest_cache[n2]:
                    dist = dest_cache[n2][w] + level
                    dest_cache[n2][n1] = dist
                    source_cache[n1][n2] = dist
                    return dist
                if w in source_cache and n2 in source_cache[w]:
                    dist = source_cache[w][n2] + level
                    dest_cache[n2][n1] = dist
                    source_cache[n1][n2] = dist
                    return dist
                # TODO: you could probably memoize this more tbh. Save the path lengths for the in-between nodes...
                source_cache[n1][w] = level
                dest_cache[w][n1] = level
                nextlevel.append(w)
    return None
