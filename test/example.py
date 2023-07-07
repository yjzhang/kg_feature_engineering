# basic usage example for kgfe

import kgfe

df = kgfe.load_graph('reactome_genes_chems.csv.gz')
graph = kgfe.df_to_graph(df)
topic_ids = ['NCBIGene::5972',
        'NCBIGene::958',
        'NCBIGene::100', 'NCBIGene::8797', 'NCBIGene::26762']
pr_results, top_nodes = kgfe.explanations.topic_pagerank(graph, topic_ids)
# graph size: 15957 nodes, 176595 edges
# timing for pr: 53 ms (also includes postprocessing)
shortest_path_lengths = graph.distances()
# timing for shortest paths: 19s

import networkx as nx
nx_graph = kgfe.graph_info.df_to_networkx(df)
nx_pr_results = nx.pagerank(nx_graph, personalization={i: 1 for i in topic_ids})
# this graph has 158959 edges (possibly removing duplicate edges between nodes?)
# note: networkx by default does not have duplicate edges between nodes (requres nx.MultiGraph class to have that)
# timing for pr: 137 ms (does not include postprocessing)
shortest_path_lengths = list(nx.shortest_path_length(nx_graph))
# timing for shortest paths: gets killed on 16gb of memory (takes at least many minutes)

# also try doing a hypergeometric test
