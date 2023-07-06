# TODO: have a graph wrapper class that wraps both networkx and igraph graphs. Use the networkx graph interface and fields and stuff.

class GraphWrapper:

    def __init__(self, graph, is_ig=True):
        """
        is_ig: is the graph an igraph Graph? If not, it's assumed to be a networkx graph.
        """
        self.graph = graph
        self.is_ig = is_ig

    @property
    def vs(self):
        if self.is_ig:
            return self.graph.vs
        else:
            return self.graph.nodes

    @property
    def nodes(self):
        if self.is_ig:
            return self.graph.vs
        else:
            return self.graph.nodes
