# API for https://academic.oup.com/bioinformatics/article/38/17/4194/6633929, https://togoid.dbcls.jp/apidoc/
import functools
import json
import urllib.parse
import urllib.request

# TODO: make this cached?
@functools.cache
def convert_ids(ids, source, dest):
    """
    List of databases for source/dest: https://academic.oup.com/view-large/401921930

    Example/commonly used databases: chebi, chembl_compound, chembl_target, clinvar, drugbank, ensembl_gene, ensembl_protein, ensembl_transcript, go (gene ontology), hgnc_symbol (gene names), hmdb (human metabolome db), inchi_key, mesh, mondo (disease), ncbigene, pubchem_compound, puchem_substance, pubmed, reactome (reactions), taxonomy, uniprot, wikipathways

    Examples:
    > convert_ids(['Q13427', 'P02790'], 'uniprot', 'ncbigene')

    Args:
        ids (list of strings)
        source (database)
        dest (database)
    """
    # URL: https://api.togoid.dbcls.jp
    # /convert
    url = 'https://api.togoid.dbcls.jp/convert'
    values = {'ids': ','.join(str(i) for i in ids),
              'route': source + ',' + dest,
              'report': 'pair'}
    url_values = urllib.parse.urlencode(values)
    full_url = url + '?' + url_values
    response = urllib.request.urlopen(full_url)
    data = json.load(response)
    return data

def convert_prots_to_genes(prot_ids):
    """
    convert uniprot IDs to NCBI gene IDs
    """
    data = convert_ids(tuple(prot_ids), 'uniprot', 'ncbigene')
    pairs = data['results']
    return pairs

def convert_prots_to_gene_symbols(prot_ids):
    """
    convert uniprot IDs to gene symbols from HGNC
    """
    data = convert_ids(tuple(prot_ids), 'uniprot', 'hgnc,hgnc_symbol')
    pairs = data['results']
    return pairs
