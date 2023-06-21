import gene_names
import pandas as pd
# TODO: load reactome data
# download links: https://reactome.org/download-data
# files: NCBI2Reactome.txt - NCBI gene IDs to reactome lowest-level pathways (which are usually but not always reactions)
# files: ChEBI2Reactome.txt - ChEBI IDs to reactome lowest-level pathways
# NCBI2ReactomeReactions.txt - NCBI gene IDs to reactions only

#ncbi_reactions_file = open('../reactome/NCBI2ReactomeReactions.txt')
#chebi_reactions_file = open('../reactome/ChEBI2ReactomeReactions.txt')
ncbi_all_levels_file = open('../reactome/NCBI2Reactome_All_Levels.txt')
chebi_all_levels_file = open('../reactome/ChEBI2Reactome_All_Levels.txt')

ncbi_entries = []
all_reactome_ids = {}

# mapping of reactome reactions and pathway ids to the name
reactome_reactions = {}
reactome_pathways = {}

def load_reactome_file(data_file, subject_category='Gene',
        subject_id_prefix='NCBIGene',
        predicate = 'participates_in',
        object_category='MolecularActivity',
        object_id_prefix='REACT',
        knowledge_source='Reactome',
        publications='PMID:34788843',
        all_reactome_ids=None,
        use_gene_name=True):
    ncbi_entries = []
    if all_reactome_ids is None:
        all_reactome_ids = {}
    for line in data_file.readlines():
        line_data = line.strip().split('\t')
        gene_id = line_data[0]
        reactome_id = line_data[1]
        reactome_name = line_data[3]
        species = line_data[5]
        if species != 'Homo sapiens' or 'HSA' not in reactome_id:
            continue
        if use_gene_name:
            try:
                gene_name = gene_names.get_symbols([int(gene_id)])[0]
            except:
                continue
        else:
            gene_name = gene_id
        all_reactome_ids[reactome_id] = reactome_name
        new_entry = {}
        new_entry['subject_category'] = subject_category
        new_entry['subject_id_prefix'] = subject_id_prefix
        new_entry['subject_id'] = gene_id
        new_entry['subject_name'] = gene_name
        new_entry['predicate'] = predicate
        new_entry['object_category'] = object_category
        new_entry['object_id_prefix'] = object_id_prefix
        new_entry['object_id'] = reactome_id
        new_entry['object_name'] = reactome_name
        new_entry['Primary_Knowledge_Source'] = reactome_id
        new_entry['Knowledge_Source'] = knowledge_source
        new_entry['publications'] = publications
        ncbi_entries.append(new_entry)
    return ncbi_entries

#ncbi_entries_reactions = load_reactome_file(ncbi_reactions_file, all_reactome_ids=all_reactome_ids, object_category='Pathway')
ncbi_entries_pathways = load_reactome_file(ncbi_all_levels_file, object_category='Pathway', all_reactome_ids=all_reactome_ids)
all_reactome_ids_chebi = {}
#chebi_entries_reactions = load_reactome_file(chebi_reactions_file, subject_category='SmallMolecule', subject_id_prefix='CHEBI',
#        all_reactome_ids=all_reactome_ids_chebi, object_category='Pathway')
chebi_entries_pathways = load_reactome_file(chebi_all_levels_file, subject_category='SmallMolecule', subject_id_prefix='CHEBI',
        object_category='Pathway',
        all_reactome_ids=all_reactome_ids_chebi, use_gene_name=False)

# get the pathway hierarchy
reactome_relations = open('../reactome/ReactomePathwaysRelation.txt')
reactome_relations_entries = []
for line in reactome_relations.readlines():
    line_data = line.split()
    parent_id = line_data[0]
    child_id = line_data[1]
    if parent_id not in all_reactome_ids or child_id not in all_reactome_ids:
        continue
    parent_name = all_reactome_ids[parent_id]
    child_name = all_reactome_ids[child_id]
    new_entry = {}
    new_entry['subject_category'] = 'Pathway'
    new_entry['subject_id_prefix'] = 'REACT'
    new_entry['subject_id'] = child_id
    new_entry['subject_name'] = child_name
    new_entry['predicate'] = 'subclass_of'
    new_entry['object_category'] = 'Pathway'
    new_entry['object_id_prefix'] = 'REACT'
    new_entry['object_id'] = parent_id
    new_entry['object_name'] = parent_name
    new_entry['Primary_Knowledge_Source'] = 'Reactome'
    new_entry['Knowledge_Source'] = 'Reactome'
    new_entry['publications'] = 'PMID:34788843'
    reactome_relations_entries.append(new_entry)

columns = 'subject_id  object_id   subject_id_prefix   object_id_prefix    subject_name    object_name predicate   Primary_Knowledge_Source    Knowledge_Source    publications    subject_category    object_category'.split()

# TODO: add PPIs?
ppi_filename = '../reactome/reactome.homo_sapiens.interactions.psi-mitab.txt'
data_ppi = pd.read_csv(ppi_filename, sep='\t')
ppi_entries = []
id_prefix_map = {'uniprotkb': 'UNIPROT', 'ChEBI': 'CHEBI', 'reactome': 'REACT'}
id_prefix_to_category = {'uniprotkb': 'Protein', 'ChEBI': 'SmallMolecule', 'reactome': 'Pathway'}
for i, row in data_ppi.iterrows():
    new_entry = {}
    new_entry['subject_category'] = id_prefix_to_category[row['#ID(s) interactor A'].split(':')[0]]
    new_entry['subject_id_prefix'] = id_prefix_map[row['#ID(s) interactor A'].split(':')[0]]
    if new_entry['subject_id_prefix'] != 'UNIPROT':
        print(new_entry['subject_id_prefix'])
    new_entry['subject_id'] = row['#ID(s) interactor A'].split(':')[1]
    new_entry['subject_name'] = row['#ID(s) interactor A'].split(':')[1]
    new_entry['predicate'] = row['Interaction type(s)']
    new_entry['object_category'] = id_prefix_to_category[row['ID(s) interactor B'].split(':')[0]]
    new_entry['object_id_prefix'] = id_prefix_map[row['ID(s) interactor B'].split(':')[0]]
    new_entry['object_id'] = row['ID(s) interactor B'].split(':')[1]
    new_entry['object_name'] = row['ID(s) interactor B'].split(':')[1]
    new_entry['Primary_Knowledge_Source'] = row['Interaction annotation(s)']
    new_entry['Knowledge_Source'] = 'Reactome'
    new_entry['publications'] = row['Publication Identifier(s)']
    ppi_entries.append(new_entry)



ncbi_table = pd.DataFrame(ncbi_entries_pathways, columns=columns)
chebi_table = pd.DataFrame(chebi_entries_pathways, columns=columns)
combined_table = pd.DataFrame(ncbi_entries_pathways + chebi_entries_pathways + reactome_relations_entries,
        columns=columns)

combined_table.to_csv('reactome_genes_chems.csv', index=False)

ppi_table = pd.DataFrame(ppi_entries, columns=columns)
ppi_table.to_csv('reactome_ppis.csv', index=False)
