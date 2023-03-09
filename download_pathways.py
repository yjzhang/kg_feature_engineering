# crawls through all kegg pathways, downloads stuff

# source: https://www.genome.jp/kegg/pathway.html

from bioservices.kegg import KEGG

s = KEGG()
s.organism = "hsa"

mapping_kegg_uniprot = s.conv("hsa", "uniprot")
mapping_kegg_ncbi = s.conv("ncbi-geneid", 'hsa')

mapping_kegg_pubchem = s.conv('pubchem', 'cpd')

mapping_kegg_drug = s.conv('pubchem', 'drug')

print('Finished mappings')
# cpd:C00162


def postprocess_data(data, kegg_id='KEGG'):
    data_id = data['ENTRY'].split()[0]
    new_entries = []
    if kegg_id == 'KEGG':
        kegg_id = data_id
    if 'GENE' in data:
        for gene_id, gene_value in data['GENE'].items():
            new_entry = {}
            new_entry['subject_category'] = 'Gene'
            new_entry['subject_id_prefix'] = 'NCBIGene'
            new_entry['subject_id'] = gene_id
            new_entry['subject_name'] = gene_value.split(';')[0]
            new_entry['predicate'] = 'participates_in'
            new_entry['object_category'] = 'Pathway'
            new_entry['object_id_prefix'] = 'KEGG'
            new_entry['object_id'] = data_id
            new_entry['object_name'] = data['NAME'][0]
            new_entry['Primary_Knowledge_Source'] = kegg_id
            new_entry['Knowledge_Source'] = 'KEGG_Pathways'
            new_entry['publications'] = 'PMID:27899662'
            new_entries.append(new_entry)
    if 'COMPOUND' in data:
        for c_id, c_value in data['COMPOUND'].items():
            new_entry = {}
            new_entry['subject_category'] = 'SmallMolecule'
            new_entry['subject_id_prefix'] = 'PubChem'
            new_entry['subject_id'] = mapping_kegg_pubchem['cpd:'+c_id].split(':')[1]
            new_entry['subject_name'] = c_value
            new_entry['predicate'] = 'participates_in'
            new_entry['object_category'] = 'Pathway'
            new_entry['object_id_prefix'] = 'KEGG'
            new_entry['object_id'] = data_id
            new_entry['object_name'] = data['NAME'][0]
            new_entry['Primary_Knowledge_Source'] = kegg_id
            new_entry['Knowledge_Source'] = 'KEGG_Pathways'
            new_entry['publications'] = 'PMID:27899662'
            new_entries.append(new_entry)
    if 'DRUG' in data:
        for drug_id, drug_value in data['DRUG'].items():
            new_entry = {}
            new_entry['subject_category'] = 'Drug'
            new_entry['subject_id_prefix'] = 'PubChem'
            new_entry['subject_id'] = mapping_kegg_drug['dr:'+drug_id].split(':')[1]
            new_entry['subject_name'] = drug_value
            new_entry['predicate'] = 'participates_in'
            new_entry['object_category'] = 'Pathway'
            new_entry['object_id_prefix'] = 'KEGG'
            new_entry['object_id'] = data_id
            new_entry['object_name'] = data['NAME'][0]
            new_entry['Primary_Knowledge_Source'] = kegg_id
            new_entry['Knowledge_Source'] = 'KEGG_Pathways'
            new_entry['publications'] = 'PMID:27899662'
            new_entries.append(new_entry)
    return new_entries


def postprocess_kgml(data, kegg_id='KEGG'):
    """
    Processes a parsed KGML to return the relations... convert to biolink format?
    """
    data['relations']
    entries_map = {}
    data['entries']
    for e in data['entries']:
        entries_map[e['id']] = e
    # create new relation for each relation.
    # list of subject_id, object_id, subject_id_prefix, object_id_prefix,
    # subject_name, object_name, predicate, Primary_Knowledge_Source,
    # Knowledge_Source, publications, subject_category, object_category
    # the link could
    new_entries = []
    for r in data['relations']:
        new_entry = {}
        e1 = entries_map[r['entry1']]
        e2 = entries_map[r['entry2']]
        if e1['type'] == 'gene':
            new_entry['subject_category'] = 'Gene'
            new_entry['subject_id_prefix'] = 'NCBIGene'
            name = e1['name'].split()[0]
            new_entry['subject_id'] = mapping_kegg_ncbi[name].split(':')[1]
            new_entry['subject_name'] = e1['gene_names'].split(', ')[0]
        elif e1['type'] == 'compound':
            new_entry['subject_category'] = 'SmallMolecule'
            new_entry['subject_id_prefix'] = 'PubChem'
            name = e1['name'].split()[0]
            new_entry['subject_id'] = mapping_kegg_pubchem[name].split(':')[1]
            new_entry['subject_name'] = e1['gene_names'].split(', ')[0]
        else:
            print('skipping non-gene or compound entry...', e1)
            continue
        if e2['type'] == 'gene':
            new_entry['object_category'] = 'Gene'
            new_entry['object_id_prefix'] = 'NCBIGene'
            name = e2['name'].split()[0]
            new_entry['object_id'] = mapping_kegg_ncbi[name].split(':')[1]
            new_entry['object_name'] = e2['gene_names'].split(', ')[0]
        elif e2['type'] == 'compound':
            new_entry['object_category'] = 'SmallMolecule'
            new_entry['object_id_prefix'] = 'PubChem'
            name = e2['name'].split()[0]
            new_entry['object_id'] = mapping_kegg_pubchem[name].split(':')[1]
            new_entry['object_name'] = e2['gene_names'].split(', ')[0]
        else:
            print('skipping non-gene or compound entry...', e2)
            continue
        new_entry['predicate'] = r['name']
        new_entry['Primary_Knowledge_Source'] = kegg_id
        new_entry['Knowledge_Source'] = 'KEGG_Pathways'
        new_entry['publications'] = 'PMID:27899662'
        new_entries.append(new_entry)
    return new_entries

    # subject_id: hsa:...


parsed_kgml = s.parse_kgml_pathway('hsa04930')
relations = postprocess_kgml(parsed_kgml, kegg_id='hsa04930')

import pandas as pd
columns = 'subject_id  object_id   subject_id_prefix   object_id_prefix    subject_name    object_name predicate   Primary_Knowledge_Source    Knowledge_Source    publications    subject_category    object_category'.split()
df = pd.DataFrame(relations, columns=columns)
df.to_csv('hsa04930.csv')

parsed_data = s.parse(s.get('hsa04930'))


pathway_strings = [s.get(x) for x in s.pathwayIds]
pathway_parsed_kgmls = [s.parse_kgml_pathway(x) for x in s.pathwayIds]
pathway_data = [s.parse(x) for x in pathway_strings]
pathway_ids = {d['ENTRY'].split()[0]: d for d in pathway_data}

# pathway relations of the form [x] PARTICIPATES_IN [pathway]
pathway_relations = []
for data in pathway_data:
    try:
        pathway_relations += postprocess_data(data)
    except:
        continue

pathway_data_table = pd.DataFrame(pathway_relations, columns=columns)
pathway_data_table.to_csv('kegg_pathway_data.csv', index=False)

# relationships within pathways - [x] (something) [y]
interaction_data = []
for x, kgml in zip(s.pathwayIds, pathway_parsed_kgmls):
    interaction_data += postprocess_kgml(parsed_kgml, kegg_id=x)

interaction_data_table = pd.DataFrame(interaction_data, columns=columns)
interaction_data_table.to_csv('kegg_interaction_data.csv', index=False)
# t2d pathway: hsa04930


