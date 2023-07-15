import time
import kgfe

t = time.time()
spoke_graph = kgfe.spoke_utils.load_spoke_igraph('../../graph_utils/spoke_2021.jsonl.gz')
print('time for loading spoke from jsonl with igraph:', time.time() - t)

