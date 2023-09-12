import time

import kgfe

# node2vec to create node embeddings of the SPOKE-2021 graph

# TODO:
spoke_graph = kgfe.spoke_loader.load_spoke_igraph('../../graph_utils/spoke_2021.jsonl.gz', directed=False, verbose=True, low_memory=False)
spoke_graph = kgfe.graph_info.largest_component(spoke_graph)

# random walks, parallel random walks
t = time.time()
random_walks = kgfe.node2vec.random_walks_igraph(spoke_graph, 10, 50, verbose=True)
print('time for random walks:', time.time() - t)
print(len(random_walks), len(random_walks[0]))
# 6126s

"""
t = time.time()
random_walks_parallel = kgfe.node2vec.random_walks_igraph_parallel(spoke_graph, 10, 50, verbose=True)
print('time for parallel random walks:', time.time() - t)
print(len(random_walks_parallel), len(random_walks_parallel[0]))
# 1.8s
"""

# use word2vec to generate the actual embeddings
t = time.time()
model1 = kgfe.node2vec.run_word2vec(random_walks)
print('time for word2vec:', time.time() - t)
# 1059s

model1.save('spoke_2021_node2vec_64')

import gensim
model1 = gensim.models.Word2Vec.load('spoke_2021_node2vec_64')
# vectors is a np array of shape nodes x 64.
vectors = model1.wv.vectors
import umap
um = umap.UMAP()
vectors_um = um.fit_transform(vectors)
# plot vectors with plotly express
import plotly.express as px
fig = px.scatter(x=vectors_um[:, 0], y=vectors_um[:, 1],
        hover_data=[[(v['name'], v['feature_name'], v['category']) for v in spoke_graph.vs]],
        color=[v['category'] for v in spoke_graph.vs])
fig.update_traces(marker=dict(size=1))
html = fig.to_html()
with open('spoke_node2vec_umap.html', 'w') as f:
    f.write(html)


