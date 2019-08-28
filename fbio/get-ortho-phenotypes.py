import argparse
import logging
import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


parser = argparse.ArgumentParser(
    description='Fetches ortholog phenotypes from list of gene symbols'
)
parser.add_argument('--input', '-i', type=str, required=True,
                    help='input list of gene symbols')
parser.add_argument('--output', '-o', type=str, required=True,
                    help='Location of output file')

args = parser.parse_args()

monarch_api = 'https://api-dev.monarchinitiative.org/api'
genes = set()
gene_ids = []
#ortholog_phenotypes = monarch_api + '/bioentity/gene/{id}/ortholog/phenotypes'

# i/o
output_file = open(args.output, 'w')

with open(args.input) as gene_list:
    for gene in gene_list:
        genes.add(gene.rstrip('\n'))

for gene in genes:
    if gene == 'ACA64':
        gene = 'SNORA78'
    # Get the gene id from the symbol
    monarch_search =  monarch_api + '/search/entity/' + gene
    params = {
        'category': 'gene',
        'taxon': 'NCBITaxon:9606',
        'rows': 1
    }
    req = requests.get(monarch_search, params=params)
    resp = req.json()
    gene_ids.append((resp['docs'][0]['id'], gene))

for gene_id, gene_symbol in gene_ids:
    # Get the gene id from the symbol
    # https://api-dev.monarchinitiative.org/api/bioentity/gene/HGNC%3A795/ortholog/phenotypes?rows=5000&facet=false&unselect_evidence=false&exclude_automatic_assertions=false&fetch_objects=false&use_compact_associations=false
    ortholog_phenotypes = monarch_api + '/bioentity/gene/' + gene_id + '/ortholog/phenotypes'
    params = {
        'rows': 1000,
        'start': 0,
        'unselect_evidence': True,
    }

    numFound = params['rows']
    while params['start'] < numFound:
        req = requests.get(ortholog_phenotypes, params=params)
        resp = req.json()
        numFound = resp['numFound']
        params['start'] += params['rows']

        # Process current set of results
        for doc in resp['associations']:
            ortholog_id = doc['subject']['id']
            ortholog_label = doc['subject']['label']
            ortholog_tax = doc['subject']['taxon']['id']
            ortholog_tax_label =  doc['subject']['taxon']['label']
            phenotype_id = doc['object']['id']
            phenotype_label = doc['object']['label']
            if doc['publications'] is not None:
                pmids = [pmid['id'] for pmid in doc['publications']]
            else:
                pmids = []

            output_file.write('\t'.join([
                gene_symbol,
                gene_id,
                ortholog_id,
                ortholog_label,
                ortholog_tax,
                ortholog_tax_label,
                phenotype_id,
                phenotype_label,
                '|'.join(pmids)
            ]) + '\n')
