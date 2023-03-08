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
