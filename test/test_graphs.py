import unittest

import kgfe

class GraphTest(unittest.TestCase):

    def setUp(self):
        self.df = kgfe.load_graph('reactome_genes_chems.csv.gz')
        self.graph = kgfe.df_to_graph(self.df)
        self.topic_ids = ['NCBIGene::5972',
                'NCBIGene::958',
                'NCBIGene::100', 'NCBIGene::8797', 'NCBIGene::26762']

    def test_pagerank(self):
        pr_results, top_nodes = kgfe.explanations.topic_pagerank(self.graph, self.topic_ids)
        self.assertTrue(len(pr_results) > 0 and len(top_nodes) > 0)
        self.assertTrue(len(pr_results) == len(self.graph.vs))
        self.assertTrue(len(top_nodes) == len(pr_results) - len(self.topic_ids))
        for node in top_nodes:
            self.assertTrue(node['score'] > 0 and node['score'] < 1)

    def test_hypergeom(self):
        hypergeom_results = kgfe.explanations.hypergeom_test(self.graph, self.topic_ids, 'Gene')
        self.assertTrue(len(hypergeom_results) > 0)
        for k, v in hypergeom_results.items():
            p, overlaps = v
            self.assertTrue(p > 0 and p < 1)
            self.assertTrue(len(overlaps) > 0)

    def test_node_stats(self):
        node_stats = kgfe.explanations.graph_node_stats(self.graph, self.topic_ids)
        pd = node_stats['average_pairwise_distance']
        # has been calculated earlier
        self.assertEqual(pd, 3.8)

    def test_null_model(self):
        null_stats = kgfe.explanations.null_graph_stats(self.graph, 'Gene', 20, 100)
        null_stats = sum(x['average_pairwise_distance'] for x in null_stats)/100.
        # average pairwise distance within groups of 20 randomly selected genes
        self.assertTrue(null_stats > 3 and null_stats < 4)

    def test_steiner_tree(self):
        st = kgfe.explanations.steiner_tree(self.graph, self.topic_ids)
        self.assertTrue(st.is_connected())
        self.assertTrue(st.is_tree())
        for topic_id in self.topic_ids:
            self.assertIsNotNone(st.vs.find(name=topic_id))
            node = st.vs.find(name=topic_id)
            self.assertEqual(len(node.neighbors()), 1)


if __name__ == '__main__':
    unittest.main()
