# steiner tree implementation in graph, mehlhorn approximation
# based on the networkx implementation at https://networkx.org/documentation/stable/_modules/networkx/algorithms/approximation/steinertree.html#steiner_tree

import igraph as ig


def metric_closure(G, weight="weight"):
    """Return the metric closure of a graph.

    The metric closure of a graph *G* is the complete graph in which each edge
    is weighted by the shortest path distance between the nodes in *G* .

    Parameters
    ----------
    G : NetworkX graph

    Returns
    -------
    NetworkX graph
        Metric closure of the graph `G`.

    """
    M = nx.Graph()

    Gnodes = set(G)

    # check for connected graph while processing first node
    all_paths_iter = nx.all_pairs_dijkstra(G, weight=weight)
    u, (distance, path) = next(all_paths_iter)
    if Gnodes - set(distance):
        msg = "G is not a connected graph. metric_closure is not defined."
        raise nx.NetworkXError(msg)
    Gnodes.remove(u)
    for v in Gnodes:
        M.add_edge(u, v, distance=distance[v], path=path[v])

    # first node done -- now process the rest
    for u, (distance, path) in all_paths_iter:
        Gnodes.remove(u)
        for v in Gnodes:
            M.add_edge(u, v, distance=distance[v], path=path[v])

    return M

# TODO: cached shortest paths?


def multi_source_shortest_paths(G, source_nodes, shortest_paths_cache=None):
    """
    Returns the shortest paths from any of the nodes in source_nodes to every node in G, keyed by the destination (nodes in G).

    Assumes that everything in source_nodes is an index into G.

    shortest_paths_cache is a dict containing the results of G.get_shortest_paths(n) - all-dest shortest paths for node n.
    """
    if shortest_paths_cache is None:
        shortest_paths_cache = {}
    shortest_paths = {}
    for n in source_nodes:
        if n in shortest_paths_cache:
            paths = shortest_paths_cache[n]
        else:
            paths = G.get_shortest_paths(n)
            shortest_paths_cache[n] = paths
        for i, path in enumerate(paths):
            if i in shortest_paths:
                if len(path) < len(shortest_paths[i]):
                    shortest_paths[i] = path
            else:
                shortest_paths[i] = path
    return shortest_paths


def mehlhorn_steiner_tree(G, terminal_nodes):
    """
    Mehlhorn algorithm:

    1. Construct graph G1, where the nodes are terminal nodes and the distances between the nodes are the shortest path in G.
    2. Find a minimum spanning tree G2 of G1.
    3. Construct a subgraph G3 of G by replacing each edge in G2 by the corresponding minimal path.
    4. Find a minimum spanning tree G4 of G3.
    5. Construct a Steiner tree G5 from G4 by deleting edges and nodes from G4, if necessary, so that no leaves in G5 are steiner vertices (that is, remove all leaf nodes that aren't part of the terminal nodes).
    """
    # 1. find all source shortest paths from the terminal nodes
    paths = multi_source_shortest_path(G, terminal_nodes)
    # 2. G1 - construct a complete graph
    shortest_terminals = {}
    distances_1 = {}
    for v in G.vs:
        index = v.index
        shortest_terminals[index] = paths[index][]
        pass

   
def takahashi_matsuyama_steiner_tree(G, terminal_nodes):
    """
    Takahashi and Matsuyama algorithm: I can't find the paper so I'm working off the wikipedia description

    1. start with one arbitrary terminal t
    2. find the terminal s closest to t, add the shortest path to the subgraph G'
    3. Find the closest terminal to any node in G', add that path to G'
    Continue until all nodes in G' have been added to the graph.
    Then return the minimum spanning tree of G', and remove all non-terminals with only one edge.

    According to https://arxiv.org/pdf/1409.8318v1.pdf, it produces smaller steiner trees (better approximation factor) than the Mehlhorn algorithm.
    """


def _mehlhorn_steiner_tree(G, terminal_nodes, weight):
    paths = nx.multi_source_dijkstra_path(G, terminal_nodes)
    d_1 = {}
    s = {}
    for v in G.nodes():
        s[v] = paths[v][0]
        d_1[(v, s[v])] = len(paths[v]) - 1

    # G1-G4 names match those from the Mehlhorn 1988 paper.
    G_1_prime = ig.Graph()
    for u, v, data in G.edges(data=True):
        su, sv = s[u], s[v]
        weight_here = d_1[(u, su)] + data.get(weight, 1) + d_1[(v, sv)]
        if not G_1_prime.has_edge(su, sv):
            G_1_prime.add_edge(su, sv, weight=weight_here)
        else:
            new_weight = min(weight_here, G_1_prime[su][sv][weight])
            G_1_prime.add_edge(su, sv, weight=new_weight)

    G_2 = G_1_prime.spanning_tree(G_1_prime, data=True)

    G_3 = ig.Graph()
    for u, v, d in G_2:
        path = ig.shortest_path(G, u, v, weight)
        for n1, n2 in pairwise(path):
            G_3.add_edge(n1, n2)

    G_3_mst = list(nx.minimum_spanning_edges(G_3, data=False))
    if G.is_multigraph():
        G_3_mst = (
            (u, v, min(G[u][v], key=lambda k: G[u][v][k][weight])) for u, v in G_3_mst
        )
    G_4 = G.edge_subgraph(G_3_mst).copy()
    _remove_nonterminal_leaves(G_4, terminal_nodes)
    return G_4.edges()


def _kou_steiner_tree(G, terminal_nodes, weight):
    # H is the subgraph induced by terminal_nodes in the metric closure M of G.
    M = metric_closure(G, weight=weight)
    H = M.subgraph(terminal_nodes)

    # Use the 'distance' attribute of each edge provided by M.
    mst_edges = nx.minimum_spanning_edges(H, weight="distance", data=True)

    # Create an iterator over each edge in each shortest path; repeats are okay
    mst_all_edges = chain.from_iterable(pairwise(d["path"]) for u, v, d in mst_edges)
    if G.is_multigraph():
        mst_all_edges = (
            (u, v, min(G[u][v], key=lambda k: G[u][v][k][weight]))
            for u, v in mst_all_edges
        )

    # Find the MST again, over this new set of edges
    G_S = G.edge_subgraph(mst_all_edges)
    T_S = nx.minimum_spanning_edges(G_S, weight="weight", data=False)

    # Leaf nodes that are not terminal might still remain; remove them here
    T_H = G.edge_subgraph(T_S).copy()
    _remove_nonterminal_leaves(T_H, terminal_nodes)

    return T_H.edges()


def _remove_nonterminal_leaves(G, terminals):
    terminals_set = set(terminals)
    for n in list(G.nodes):
        if n not in terminals_set and G.degree(n) == 1:
            G.remove_node(n)


def steiner_tree(graph, node_ids):
    pass
