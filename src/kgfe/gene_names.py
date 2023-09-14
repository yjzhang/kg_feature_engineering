# extremely simple way of converting gene names
import gzip
import os

base_dir = os.path.dirname(os.path.abspath(__file__))
ID_TO_SYMBOL = {}
SYMBOL_TO_ID = {}
ID_TO_UNIPROT = {}
UNIPROT_TO_ID = {}
ID_TO_ENSEMBL = {}
ENSEMBL_TO_ID = {}


def _load_gene_info():
    with gzip.open(os.path.join(base_dir, 'Homo_sapiens.gene_info.gz'), 'rt') as f:
        for row in f.readlines():
            row = row.split('\t')
            if row[0] == '#tax_id':
                continue
            gene_id = int(row[1])
            ID_TO_SYMBOL[gene_id] = row[2]
            if row[2] not in SYMBOL_TO_ID:
                SYMBOL_TO_ID[row[2]] = gene_id

def _load_uniprot_info():
    with open(os.path.join(base_dir, 'gene_uniprot.txt')) as f:
        for row in f.readlines():
            row = row.split()
            gene_id = int(row[0])
            ID_TO_UNIPROT[gene_id] = row[1]
            UNIPROT_TO_ID[row[1]] = gene_id

def _load_ensembl_info():
    with open(os.path.join(base_dir, 'geneid_ensembl.txt')) as f:
        for row in f.readlines():
            row = row.split()
            gene_id = int(row[0])
            ID_TO_ENSEMBL[gene_id] = row[1]
            ENSEMBL_TO_ID[row[1]] = gene_id

_load_gene_info()
    
def convert(source, dest, ids):
    """
    General conversion function.

    source/dest can be 'ncbi'/'geneid', 'symbol', 'uniprot', 'ensembl'
    """
    gene_ids = []
    if source == 'ncbi' or source == 'geneid':
        gene_ids = ids
    elif source == 'symbol':
        gene_ids = get_ids(ids)
    elif source == 'uniprot':
        gene_ids = uniprot_to_gene_ids(ids)
    elif source == 'ensembl':
        gene_ids = ensembl_to_gene_ids(ids)
    if dest == 'ncbi' or dest == 'geneid':
        return gene_ids
    elif dest == 'symbol':
        return get_symbols(gene_ids)
    elif dest == 'uniprot':
        return gene_ids_to_uniprot(gene_ids)
    elif dest == 'ensembl':
        return gene_ids_to_ensembl(gene_ids)

def get_symbols(gene_ids):
    """
    Get symbols given a list of gene ids.
    """
    return [ID_TO_SYMBOL[int(x)] for x in gene_ids]

def get_symbol(gene_id):
    return ID_TO_SYMBOL[int(gene_id)]

def get_ids(gene_symbols):
    """
    Get ids given a list of gene symbols.
    """
    return [SYMBOL_TO_ID[x] for x in gene_symbols]

def get_id(gene_symbol):
    return SYMBOL_TO_ID[gene_symbol]


def gene_ids_to_uniprot(gene_ids):
    if len(ID_TO_UNIPROT) == 0:
        _load_uniprot_info()
    return [ID_TO_UNIPROT[int(x)] for x in gene_ids]


def uniprot_to_gene_ids(uniprot_ids):
    if len(ID_TO_UNIPROT) == 0:
        _load_uniprot_info()
    return [UNIPROT_TO_ID[x] for x in uniprot_ids]


def gene_ids_to_ensembl(gene_ids):
    if len(ID_TO_ENSEMBL) == 0:
        _load_ensembl_info()
    return [ID_TO_ENSEMBL[int(x)] for x in gene_ids]


def ensembl_to_gene_ids(ensembl_ids):
    if len(ID_TO_ENSEMBL) == 0:
        _load_ensembl_info()
    return [ENSEMBL_TO_ID[x] for x in ensembl_ids]
