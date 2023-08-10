# basic usage example for kgfe
import time

import pandas as pd

import kgfe

# timing ran on intel i7-9700 CPU @ 3.00GHz with 64gb memory

df = kgfe.load_graph('reactome_genes_chems.csv.gz')
df2 = kgfe.load_graph('uniprot_genes.csv')
combined_df = pd.concat([df, df2], axis=0)
t = time.time()
graph = kgfe.df_to_graph(combined_df)
graph.simplify()
print('igraph graph from df time:', time.time() - t)
# timing: 12s
topic_ids = ['NCBIGene::5972',
        'NCBIGene::958',
        'NCBIGene::100', 'NCBIGene::8797', 'NCBIGene::26762']

for i in topic_ids:
    neighbors = graph.neighbors(i)
    has_protein = False
    for n in neighbors:
        if graph.vs[n]['category'] == 'Protein':
            has_protein = True
    if not has_protein:
        print('WARNING: gene {0} did not have a corresponding protein'.format(i))
