import functools
import multiprocessing as mp

#from numba import jit
import numpy as np

# TODO: rewrite in igraph, parallel
def random_walks_igraph(graph, r, l, verbose=False, parallel=False):
    """
    Random walks starting from all nodes in the input graph.
    r = number of walks per node
    l = length of walk

    Returns a list of lists of node IDs.
    """
    walks = []
    for i, v in enumerate(graph.vs):
        if verbose:
            if i % 10000 == 0:
                print('nodes walked:', i)
        for _ in range(r):
            walk = graph.random_walk(v, l)
            walks.append(walk)
    return walks

def biased_random_walks_igraph(graph, r, l, p=1, q=1, verbose=False):
    """
    Random walk starting from all nodes.
    r = number of walks per node
    l = length of walk
    p and q are probs
    """
    walks = []
    for i, v in enumerate(graph.vs):
        if verbose:
            if i % 10000 == 0:
                print('nodes walked:', i)
        for _ in range(r):
            walk = graph.random_walk(v, l)
            walks.append(walk)
    return walks

def random_walks(adj_list, r, l, p=1, q=1, verbose=False):
    """
    Biased random walk starting from node i.
    r = number of walks per node
    l = length of walk
    p and q are probs
    """
    @functools.lru_cache(10000)
    def tr_probs(prev, node):
        """
        node2vec transition probabilities
        returns array of transition probs indexed to neighbors
        """
        neighbors = adj_list[node]
        transition_probs = np.ones(len(neighbors))
        if p == 1 and q == 1:
            return transition_probs/transition_probs.sum()
        prev_neighbors = adj_list[prev]
        # bias
        for i, t in enumerate(neighbors):
            if t == prev:
                transition_probs[i] = 1./p
            elif t not in prev_neighbors:
                transition_probs[i] = 1./q
        transition_probs /= transition_probs.sum()
        return transition_probs
    walks = []
    # TODO: parallelize this?
    for start in range(len(adj_list)):
        if verbose and start % 1000 == 0:
            print('Node: ', start)
        for _ in range(r):
            n_walk = 0
            prev_node = None
            node = start
            walk = []
            while n_walk < l:
                n_walk += 1
                walk.append(node)
                neighbors = adj_list[node]
                if prev_node is not None:
                    tr = tr_probs(prev_node, node)
                else:
                    tr = np.ones(len(neighbors))/len(neighbors)
                prev_node = node
                node = np.random.choice(neighbors, p=tr)
            walks.append(walk)
    return walks

def run_word2vec(walks, k=10, d=64):
    """
    Uses gensim to run word2vec.
    walks = list of random walks, returned by random_walks.
    k = context size (word2vec parameter)
    d = output dimensionality
    """
    import gensim.models
    model = gensim.models.Word2Vec(sentences=walks, vector_size=d,
            window=k)
    return model
