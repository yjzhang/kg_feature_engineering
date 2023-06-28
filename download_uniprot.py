import pandas as pd
from kgfe import gene_names

# TODO: create a graph from uniprot node-edge relationships.
path = '../uniprot/HUMAN_9606_idmapping.dat.gz'

data = pd.read_csv(path, sep='\t', names=['uniprot', 'ID_type', 'ID'])

# TODO: get the graph
uniprotkb_ids = {}
graph_entries = []
for i, row in data.iterrows():
    if row.ID_type == 'UniProtKB-ID':
        uniprotkb_ids[row.uniprot] = row.ID
    if row.ID_type != 'GeneID':
        continue
    new_entry = {}
    new_entry['subject_category'] = 'Gene'
    new_entry['subject_id_prefix'] = 'NCBIGene'
    new_entry['subject_id'] = row.ID
    new_entry['subject_name'] = gene_names.get_symbol(int(row.ID))
    new_entry['predicate'] = 'translates_to'
    new_entry['object_category'] = 'Protein'
    new_entry['object_id_prefix'] = 'UNIPROT'
    new_entry['object_id'] = row.uniprot
    new_entry['object_name'] = uniprotkb_ids[row.uniprot]
    new_entry['Primary_Knowledge_Source'] = 'UNIPROT'
    new_entry['Knowledge_Source'] = 'UNIPROT'
    new_entry['publications'] = 'PMID:36408920'
    graph_entries.append(new_entry)

columns = 'subject_id  object_id   subject_id_prefix   object_id_prefix    subject_name    object_name predicate   Primary_Knowledge_Source    Knowledge_Source    publications    subject_category    object_category'.split()

uniprot_table = pd.DataFrame(graph_entries, columns=columns)
uniprot_table.to_csv('uniprot_genes.csv', index=False)

