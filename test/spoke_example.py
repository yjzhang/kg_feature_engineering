import time
import kgfe

# desktop cpu: Intel(R) Core(TM) i7-9700 CPU @ 3.00GHz with 64gb memory
# laptop cpu: 12th Gen Intel(R) Core(TM) i7-1265U, running through WSL2, 16gb memory

t = time.time()
spoke_graph = kgfe.spoke_loader.load_spoke_igraph('../../graph_utils/spoke_2021.jsonl.gz', directed=False, verbose=True, low_memory=True)
print('time for loading spoke from jsonl with igraph:', time.time() - t)
# time for loading spoke from jsonl with igraph: 138s - 147s on desktop
# laptop: 161.52715277671814

proteins = kgfe.graph_info.nodes_in_category(spoke_graph, 'Protein')

# time null model
t = time.time()
kgfe.explanations.null_graph_stats(spoke_graph, 'Protein', 20, 100, method='distances')
print('time for null model with 20 random proteins (distances method):', time.time() - t)
# desktop: 170.5931851863861
# laptop: 171.632248878479

t = time.time()
kgfe.explanations.null_graph_stats(spoke_graph, 'Protein', 30, 100, method='distances')
print('time for null model with 30 random proteins (distances method):', time.time() - t)
# desktop:250 
# laptop: 237.89269495010376

t = time.time()
kgfe.explanations.null_graph_stats(spoke_graph, 'Protein', 20, 100, method='get_shortest_paths')
print('time for null model with 20 random proteins (get_shortest_paths method):', time.time() - t)
# desktop: 
# laptop: 544.4069104194641

t = time.time()
kgfe.explanations.null_graph_stats(spoke_graph, 'Protein', 30, 100, method='get_shortest_paths')
print('time for null model with 30 random proteins (get_shortest_paths method):', time.time() - t)
# desktop: 
# laptop: 898.4321784973145
