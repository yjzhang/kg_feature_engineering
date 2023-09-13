# extremely simple way of converting gene names
import gzip
import os

base_dir = os.path.dirname(os.path.abspath(__file__))
ID_TO_SYMBOL = {}
SYMBOL_TO_ID = {}
SYMBOL_TO_UNIPROT = {}


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


_load_gene_info()
    

def get_symbols(gene_ids):
    """
    Get symbols given a list of gene ids.
    """
    return [ID_TO_SYMBOL[x] for x in gene_ids]

def get_symbol(gene_id):
    return ID_TO_SYMBOL[gene_id]

def get_ids(gene_symbols):
    """
    Get ids given a list of gene symbols.
    """
    return [SYMBOL_TO_ID[x] for x in gene_symbols]

def get_id(gene_symbol):
    return SYMBOL_TO_ID[gene_symbol]

# TODO: convert genes to uniprot.
def gene_symbols_to_uniprot(gene_symbols):
    pass
