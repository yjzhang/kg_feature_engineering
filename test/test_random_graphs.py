import unittest

import igraph as ig

import kgfe

class GraphTest(unittest.TestCase):

    def setUp(self):
        # TODO: generate a random graph
        self.graph = ig.Graph.Barabasi(1000, 4)
        for v in self.graph.vs:
            v['name'] = str(v.index)
            v['category'] = 'Gene'
        self.topic_ids = ['1', '2', '3', '4']

    def test_pagerank(self):
        pr_results, top_nodes = kgfe.explanations.topic_pagerank(self.graph, self.topic_ids)
        self.assertTrue(len(pr_results) > 0 and len(top_nodes) > 0)
        self.assertTrue(len(pr_results) == len(self.graph.vs))
        self.assertTrue(len(top_nodes) == len(pr_results) - len(self.topic_ids))
        for node in top_nodes:
            self.assertTrue(node['score'] > 0 and node['score'] < 1)


    def test_node_stats(self):
        node_stats = kgfe.explanations.graph_node_stats(self.graph, self.topic_ids)
        pd = node_stats['average_pairwise_distance']
        print(pd)
        # has been calculated earlier

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

    def test_steiner_tree_2(self):
        for i in range(50):
            topic_ids = kgfe.graph_info.random_nodes(self.graph, 10)
            st = kgfe.explanations.steiner_tree(self.graph, topic_ids, method='takahashi')
            self.assertTrue(st.is_connected())
            self.assertTrue(st.is_tree())
            for topic_id in topic_ids:
                self.assertIsNotNone(st.vs.find(name=topic_id))
        for i in range(50):
            topic_ids = kgfe.graph_info.random_nodes(self.graph, 10)
            st = kgfe.explanations.steiner_tree(self.graph, topic_ids, method='shortest_paths')
            self.assertTrue(st.is_connected())
            self.assertTrue(st.is_tree())
            for topic_id in topic_ids:
                self.assertIsNotNone(st.vs.find(name=topic_id))
        for i in range(50):
            topic_ids = kgfe.graph_info.random_nodes(self.graph, 10)
            st = kgfe.explanations.steiner_tree(self.graph, topic_ids, method='mehlhorn')
            self.assertTrue(st.is_connected())
            self.assertTrue(st.is_tree())
            for topic_id in topic_ids:
                self.assertIsNotNone(st.vs.find(name=topic_id))




if __name__ == '__main__':
    unittest.main()
