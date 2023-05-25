import gene_names
# TODO: load reactome data
# download links: https://reactome.org/download-data
# files: NCBI2Reactome.txt - NCBI gene IDs to reactome lowest-level pathways (which are usually but not always reactions)
# files: ChEBI2Reactome.txt - ChEBI IDs to reactome lowest-level pathways
# NCBI2ReactomeReactions.txt - NCBI gene IDs to reactions only

ncbi_reactions_file = open('../reactome/NCBI2ReactomeReactions.txt')
chebi_reactions_file = open('../reactome/ChEBI2ReactomeReactions.txt')
ncbi_all_levels_file = open('../reactome/NCBI2Reactome_All_Levels.txt')
chebi_all_levels_file = open('../reactome/ChEBI2Reactome_All_Levels.txt')

ncbi_entries = []
all_reactome_ids = {}

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
        all_reactome_ids[reactome_id] = reactome_name
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

ncbi_entries_reactions = load_reactome_file(ncbi_reactions_file, all_reactome_ids=all_reactome_ids)
ncbi_entries_pathways = load_reactome_file(ncbi_all_levels_file, object_category='Pathway', all_reactome_ids=all_reactome_ids)
all_reactome_ids_chebi = {}
chebi_entries_reactions = load_reactome_file(chebi_reactions_file, subject_category='SmallMolecule', subject_id_prefix='CHEBI',
        all_reactome_ids=all_reactome_ids_chebi)
chebi_entries_pathways = load_reactome_file(chebi_reactions_file, subject_category='SmallMolecule', subject_id_prefix='CHEBI',
        object_category='Pathway',
        all_reactome_ids=all_reactome_ids_chebi)
