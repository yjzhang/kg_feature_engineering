import math

import kgfe
import networkx as nx
import unittest

class GraphUtilsTest(unittest.TestCase):

    def setUp(self):
        pass

    def test_shortest_path_lengths_ba_100(self):
        random_graph = nx.random_graphs.barabasi_albert_graph(100, 3)
        paths_nx = dict(nx.shortest_path_length(random_graph))
        paths_kgfe = kgfe.graph_utils.all_pairs_shortest_path_lengths(random_graph, list(random_graph.nodes))
        for k, v in paths_kgfe.items():
            self.assertEqual(v, paths_nx[k])

    def test_shortest_path_lengths_ba_1000(self):
        random_graph = nx.random_graphs.barabasi_albert_graph(1000, 3)
        paths_nx = dict(nx.shortest_path_length(random_graph))
        paths_kgfe = kgfe.graph_utils.all_pairs_shortest_path_lengths(random_graph, list(random_graph.nodes))
        for k, v in paths_kgfe.items():
            self.assertEqual(v, paths_nx[k])


    def test_shortest_path_lengths_er_1000(self):
        random_graph = nx.random_graphs.erdos_renyi_graph(1000, 0.5)
        paths_nx = dict(nx.shortest_path_length(random_graph))
        paths_kgfe = kgfe.graph_utils.all_pairs_shortest_path_lengths(random_graph, list(random_graph.nodes))
        for k, v in paths_kgfe.items():
            for k1, v1 in v.items():
                if not math.isinf(v1):
                    self.assertEqual(v1, paths_nx[k][k1])




if __name__ == '__main__':
    unittest.main()
