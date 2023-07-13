import igraph as ig
import kgfe
from kgfe import steiner_tree

g = ig.Graph(n=4, edges=[[0,1], [1,2], [2,3], [3,0]])

paths = steiner_tree.multi_source_shortest_paths(g, [0, 1])
