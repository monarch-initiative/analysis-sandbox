import csv
import argparse
import logging
import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Custom session for services that require more than 3 retries
session = requests.Session()
adapter = requests.adapters.HTTPAdapter(max_retries=10)
session.mount('https://', adapter)

"""
Search for lethal phenotypes from set of genotypes
"""
SOLR_URL = 'https://solr.monarchinitiative.org/solr/golr/select'

def main():
    parser = argparse.ArgumentParser(
        description='Search for lethal phenotypes from set of genotypes')
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='Location of input file')
    parser.add_argument('--output', '-o', type=str, required=True,
                        help='Location of output file')

    args = parser.parse_args()
    output_file = open(args.output, 'w')
    output_file.write("\t".join(['genotype', 'genotype_label', 'phenotype',
                                 'phenotype_label', 'gene', 'gene_label',
                                 'sources', 'ortholog_id',
                                 'ortholog_label']) + "\n")

    genotype_set = set()

    with open(args.input, 'r') as tsvfile:
        reader = csv.reader(tsvfile, delimiter="\t")
        for row in reader:
            (genotype, zygosity, allele, gene) = row
            genotype = genotype.replace(
                "http://www.informatics.jax.org/accession/", '')
            genotype_set.add(genotype)

    for genotype in genotype_set:
        params = {
            'wt': 'json',
            'rows': 100,
            'start': 0,
            'q': '*:*',
            'fl': 'subject, subject_label, subject_gene, subject_gene_label, '
                  'object, object_label, source',
            'fq': ['subject_closure:"{0}"'.format(genotype),
                   'object_category:"phenotype"',
                   'object_closure: "{0}" OR object_closure: "{1}"'
                   ' OR object_closure: "{2}" OR object_closure: "{3}"'
                   .format("MP:0011400", "MP:0010831", "MP:0008569",
                           "MP:0010770")
            ]
        }

        resultCount = params['rows']

        while params['start'] < resultCount:
            solr_request = session.get(SOLR_URL, params=params)
            response = solr_request.json()
            resultCount = response['response']['numFound']
            for doc in response['response']['docs']:
                source = "|".join(doc['source'])

                if 'subject_gene' not in doc:
                    subject_gene = ''
                    subject_gene_label = ''
                else:
                    subject_gene = doc['subject_gene']
                    subject_gene_label = doc['subject_gene_label']

                ortholog_id, ortholog_label = get_human_ortholog(subject_gene)

                output_file.write("\t".join([doc['subject'], doc['subject_label'],
                                             doc['object'], doc['object_label'],
                                             subject_gene,
                                             subject_gene_label,
                                             source, ortholog_id,
                                             ortholog_label]) + "\n")

            params['start'] += params['rows']

def get_human_ortholog(gene):
    params = {
            'wt': 'json',
            'rows': 100,
            'start': 0,
            'q': '*:*',
            'fl': 'subject, subject_label,'
                  'object, object_label',
            'fq': ['subject_closure: "{0}"'.format(gene),
                   'relation_closure: "RO:HOM0000017"',
                   'object_taxon: "NCBITaxon:9606"'
            ]
    }
    solr_request = session.get(SOLR_URL, params=params)
    response = solr_request.json()
    resultCount = response['response']['numFound']
    if resultCount > 1:
        logger.info("More than one ortholog found for {}".format(gene))
    gene_ids = []
    gene_labels = []
    for doc in response['response']['docs']:
        gene_ids.append(doc['object'])
        gene_labels.append(doc['object_label'])

    return "|".join(gene_ids), "|".join(gene_labels)


if __name__ == "__main__":
    main()






