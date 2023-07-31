import functools

#from numba import jit
import numpy as np
from scipy import sparse

# TODO: rewrite in igraph, parallel

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


if __name__ == '__main__':
    # TODO: try this with a smaller graph?
    import spoke_loader
    import scipy.io
    import umap
    # load graph
    nodes, edges, node_types, edge_types, edge_matrix = spoke_loader.load_spoke('spoke.csv', remove_unused_nodes=True)
    edge_matrix = scipy.io.mmread('spoke.mtx')
    print('matrix loaded')
    edge_matrix = spoke_loader.symmetrize_matrix(edge_matrix)
    edge_matrix = sparse.lil_matrix(edge_matrix)
    print('calculating random walks...')
    walks = random_walks(edge_matrix.rows, r=10, l=50, verbose=True)
    n2v_model = run_word2vec(walks, 8, 50)
    n2v_model.save('spoke_node2vec_gensim_50')
    um = umap.UMAP()
    um.fit_transform(n2v_model.wv.vectors)
    np.savetxt('spoke_umap.txt', um.embedding_)
    # make a 2d plot of spoke nodes, colored by type
    nodes_types = [n[2] for n in nodes]
    # plot using plotly?
    import plotly.express as px
    fig = px.scatter(x=um.embedding_[:, 0], y=um.embedding_[:,1], hover_data=[[(n[1], node_types[n[2]]) for n in nodes]], color=[node_types[n[2]] for n in nodes])
    fig.update_traces(marker=dict(size=1))
    html = fig.to_html()
    with open('spoke.html', 'w') as f:
        f.write(html)
    np.savetxt('spoke_small_nodes.txt', np.array([[str(x) for x in n] for n in nodes]), fmt='%s', delimiter='\t')
