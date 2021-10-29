import argparse
import logging
import requests
from json import JSONDecodeError

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(
    description='Given gene labels/descriptions fetch the symbol'
)
parser.add_argument('--input', '-i', type=str, required=True,
                    help='input list of gene symbols')
parser.add_argument('--output', '-o', type=str, default='genes.tsv',
                    help='Location of output file')

args = parser.parse_args()

monarch_api = 'https://api-dev.monarchinitiative.org/api'

# i/o
output_file = open(args.output, 'w')

columns = '\t'.join(['Target.based.Actions', 'total', 'match label', 'match id', 'matched text'])

output_file.write(f"#{columns}\n")

biolink_search = monarch_api + "/search/entity/"

parameters = {
    'rows': 1,
    'category': 'gene',
    'taxon': 'NCBITaxon:9606',
    'min_match': '67%'
}

with open(args.input, 'r') as gene_file:
    for row in gene_file:
        target, total = row.strip().split('\t')

        if total == 'total': continue  # Header

        targets = [tar.strip() for tar in target.split(';')]

        ids = []
        matched_labels = []
        matched_text = []

        for term in targets:
            # strip last word in each term?
            tokens = term.split(' ')
            trimmed_term = ' '.join(tokens[:-1])
            url = biolink_search + trimmed_term
            biolink_request = requests.get(url, params=parameters)
            try:
                biolink_response = biolink_request.json()
            except JSONDecodeError:
                # TODO URL encode terms for slashes
                print(term)
                continue
            if not biolink_response['docs']:
                continue  # TODO
            matched_labels.append(biolink_response['docs'][0]['label'][0])
            match_id = biolink_response['docs'][0]['id']
            ids.append(match_id)
            matched_text.append(biolink_response['highlighting'][match_id]['match'])

        output_file.write(
            '\t'.join(
                [
                    target,
                    total,
                    ';'.join(ids),
                    ';'.join(matched_labels),
                    ';'.join(matched_text)
                ]
            ) + '\n'
        )
