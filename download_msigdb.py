# TODO: process msigdb into a single formatted csv (or a collection of csvs)

# data source: http://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp
import json
import os

import gene_names

BASE_PATH = '/media/yjzhang/easystore-5gb-1/research_big_data/msigdb/msigdb_v2023.1.Hs_json_files_to_download_locally'
files = [x for x in os.listdir(BASE_PATH) if x.endswith('.json')]
# TODO: object categories for every msigdb category, also need predicates and object_id_prefix
files_to_object_category = {}
for f in files:
    pass

def load_file(json_filename, object_category='BiologicalProcess', object_id_prefix='MSigDB',
        predicate='participates_in'):
    with open(json_filename) as f:
        data = json.load(f)
    new_entries = []
    for k, v in data.items():
        gene_set_name = k
        source = v['exactSource']
        ref = v['pmid']
        genes = v['geneSymbols']
        gene_ids = gene_names.get_ids(genes)
        for gene_id, symbol in zip(gene_ids, genes):
            entry = {}
            entry['subject_category'] = 'Gene'
            entry['subject_id_prefix'] = 'NCBIGene'
            entry['subject_id'] = gene_id
            entry['subject_name'] = symbol
            entry['predicate'] = predicate
            entry['object_category'] = object_category
            entry['object_id_prefix'] = object_id_prefix
            entry['object_id'] = gene_set_name
            entry['object_name'] = gene_set_name
            entry['Primary_Knowledge_Source'] = source
            entry['Knowledge_Source'] = 'KEGG_Pathways'
            entry['publications'] = ref
            new_entries.append(entry)
    return new_entries

