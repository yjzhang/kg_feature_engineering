A framework for feature engineering using knowledge graphs

Data ingestion:
1. Download raw knowledge graphs as csvs from various sources - see raw_graphs/
2. Process these graphs into the biolink model format, with defined columns. Store these processed graphs as CSVs in processed_graphs/.

Data usage:
1. Load graphs from csvs. Have some standard graph representation.
2. Map features in data to graphs. Create a tool that lets us convert between different feature formats - NCBI gene IDs, uniprot IDs, ensembl IDs, for metabolites/chems - HMDB ids, chembl ids, etc.

Integration:
- use Bioservices for ID mapping? Or store the ID maps locally...
- use NetworkX? Or a bespoke graph implementation? Probably just use NetworkX. Use the from_pandas_edgelist function to import the graphs.


BioLink format:
"subject_id":[],
"object_id":[],
"subject_id_prefix":[],
"object_id_prefix":[],
"subject_name":[],
"object_name":[],
"predicate":[],
"Primary_Knowledge_Source":[],
"Knowledge_Source":[],
"publications":[],
"subject_category":[],
"object_category":[]

## Installation

Run `pip install -e .` in this directory.

To run tests, run `python -m unittest discover test`

## Usage

```
import kgfe

# get available knowledge graphs
available_graphs = kgfe.get_available_graphs()

pathway_graph_df = kgfe.load_graph('kegg_pathway_data.csv')
pathway_graph = kgfe.df_to_networkx(graph_df)
nodes_table = kgfe.get_nodes_table(pathway_graph)

# using graphs for explanation
# pagerank results
pr_results = kgfe.explanations.topic_pagerank(pathway_graph, query_genes, 'Genes')
# hypergeom results
```

TODO: include more graphs - MeshGeneGraph, CellMesh, MSigDB 
