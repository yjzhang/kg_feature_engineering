import unittest

import kgfe

class GraphTest(unittest.TestCase):

    def setUp(self):
        pass

    def testLoadGraph(self):
        available_graphs = kgfe.get_available_graphs()
        self.assertTrue(len(available_graphs) > 0)
        self.assertTrue('kegg_pathway_data.csv' in available_graphs)
        df = kgfe.load_graph('kegg_pathway_data.csv')
        graph = kgfe.df_to_networkx(df)
        self.assertIsNotNone(graph)

    def test_hypergeom(self):
        df = kgfe.load_graph('kegg_pathway_data.csv')
        graph = kgfe.df_to_networkx(df)
        hg_results = kgfe.explanations.hypgergeom_test(graph, [7124, 3667, 3630], 'Gene')
        self.assertTrue(len(hg_results) > 0)

    def test_pagerank(self):
        df = kgfe.load_graph('kegg_pathway_data.csv')
        graph = kgfe.df_to_networkx(df)
        pr_results = kgfe.explanations.topic_pagerank(graph, [6504, 5972, 958, 100, 8797, 26762], 'Gene')
        self.assertTrue(len(pr_results) == len(graph.nodes))

if __name__ == '__main__':
    unittest.main()
