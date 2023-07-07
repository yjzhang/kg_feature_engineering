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
        graph = kgfe.df_to_graph(df)
        self.assertIsNotNone(graph)

    def test_pagerank(self):
        df = kgfe.load_graph('reactome_genes_chems.csv.gz')
        graph = kgfe.df_to_graph(df)
        topic_ids = ['NCBIGene::5972',
                'NCBIGene::958',
                'NCBIGene::100', 'NCBIGene::8797', 'NCBIGene::26762']
        pr_results, top_nodes = kgfe.explanations.topic_pagerank(graph, topic_ids)
        self.assertTrue(len(pr_results) > 0 and len(top_nodes) > 0)
        self.assertTrue(len(top_nodes) == len(pr_results))
        self.assertTrue(top_nodes[0]['score'] > 0 and top_nodes[0]['score'] < 1)

    def test_hypergeom(self):
        df = kgfe.load_graph('kegg_pathway_data.csv')
        graph = kgfe.df_to_graph(df)
        topic_ids = ['NCBIGene::5972',
                'NCBIGene::958',
                'NCBIGene::100', 'NCBIGene::8797', 'NCBIGene::26762']
        hypergeom_results = kgfe.explanations.hypgergeom_test(graph, topic_ids, 'Gene')
        self.assertTrue(len(hypergeom_results) > 0)

if __name__ == '__main__':
    unittest.main()
