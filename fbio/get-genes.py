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
parser.add_argument('--num-match', '-n', type=int, default=1,
                    help='Number of matches to return')

args = parser.parse_args()

monarch_api = 'https://api-dev.monarchinitiative.org/api'

# i/o
output_file = open(args.output, 'w')

cols = ['Target.based.Actions', 'total']

for num_match in range(1, args.num_match + 1):
    cols = cols + [f'{num_match} match label', f'{num_match} match id', f'{num_match} match text']

columns = '\t'.join(cols)


output_file.write(f"#{columns}\n")

biolink_search = monarch_api + "/search/entity/"

parameters = {
    'rows': args.num_match,
    'category': 'gene',
    'taxon': 'NCBITaxon:9606',
    'min_match': '67%'
}
parameters_no_min = {
    'rows': args.num_match,
    'category': 'gene',
    'taxon': 'NCBITaxon:9606',
}

with open(args.input, 'r') as gene_file:
    for row in gene_file:
        target, total = row.strip().split('\t')

        if total == 'total': continue  # Header

        targets = [tar.strip() for tar in target.split(';')]

        for term in targets:

            output_file.write(
                '\t'.join(
                    [
                        term,
                        total
                    ]
                ) + '\t'
            )
            # strip last word in each term?
            tokens = term.split(' ')
            # todo urllib.parse.urlencode
            trimmed_term = ' '.join(tokens[:-1]).replace("/", "%2F")
            url = biolink_search + trimmed_term
            biolink_request = requests.get(url, params=parameters)
            try:
                biolink_response = biolink_request.json()
            except JSONDecodeError:
                # TODO URL encode terms for slashes
                print(trimmed_term)
                output_file.write(''.join(['\t' for _ in range(args.num_match*3)]))
                continue
            if not biolink_response['docs']:
                biolink_request = requests.get(url, params=parameters_no_min)
                biolink_response = biolink_request.json()
                if not biolink_response['docs']:
                    output_file.write(''.join(['\t' for _ in range(args.num_match*3)]))
                    continue

            for index, doc in enumerate(biolink_response['docs']):
                matched_label = doc['label'][0]
                match_id = doc['id']
                matched_text = biolink_response['highlighting'][match_id]['match']

                output_file.write(
                    '\t'.join(
                        [
                            matched_label,
                            match_id,
                            matched_text
                        ]
                    )
                )
                if index != args.num_match - 1:
                    output_file.write('\t')

            num_found = biolink_response['numFound']
            if num_found < args.num_match:
                output_file.write(''.join(['\t' for _ in range((args.num_match - num_found)*3)]))
        output_file.write('\n')
