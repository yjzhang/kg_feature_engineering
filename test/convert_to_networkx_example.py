import igraph as ig

vertices = [{'name': 'apple'}, {'name': 'pear'}, {'name': 'peach'}]
edges = [{'source': 'apple', 'target': 'pear', 'weight': 1.2}, {'source': 'apple', 'target': 'peach', 'weight': 0.9}]
g = ig.Graph.DictList(vertices, edges)

nxg = g.to_networkx()
