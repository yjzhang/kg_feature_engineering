import kgfe
import networkx as nx
import unittest

class GraphUtilsTest(unittest.TestCase):

    def setUp(self):
        pass

    def test_shortest_path_lengths(self):
        random_graph = nx.random_graphs.barabasi_albert_graph(100, 3)
        all_pairs_paths = {}
        for n1 in range(100):
            for n2 in range(n1 + 1, 100):
                path = nx.shortest_path_length(random_graph, n1, n2)
                all_pairs_paths[(n1, n2)] = path
        kgfe.graph_utils.all_pairs_shortest_path_lengths(random_graph, random_graph.nodes)


if __name__ == '__main__':
    unittest.main()
