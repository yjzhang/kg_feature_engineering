# API for https://academic.oup.com/bioinformatics/article/38/17/4194/6633929, https://togoid.dbcls.jp/apidoc/
import functools
import urllib.request

# TODO: make this cached?
@functools.cache
def convert_ids(ids, source, dest):
    """
    List of databases for source/dest: https://academic.oup.com/view-large/401921930

    Example/commonly used databases: chebi, chembl_compound, chembl_target, clinvar, drugbank, ensembl_gene, ensembl_protein, ensembl_transcript, go (gene ontology), hmdb (human metabolome db), inchi_key, mesh, mondo (disease), ncbigene, pubchem_compound, puchem_substance, pubmed, reactome (reactions), taxonomy, uniprot, wikipathways

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
