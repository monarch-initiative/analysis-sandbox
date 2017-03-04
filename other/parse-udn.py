import json
import argparse
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

"""
Parse UDN patient dump file
"""

parser = argparse.ArgumentParser(
    description='Converts UDN patient dump file into table'
)
parser.add_argument('--input', '-i', type=str, required=True,
                    help='Location of input file')
parser.add_argument('--output', '-o', type=str, required=True,
                    help='Location of output file')

args = parser.parse_args()

input_file = open(args.input, 'r')
output_file = open(args.output, 'w')

json_data = json.load(input_file)

count = 0

for entry in json_data:
    simple_id = entry['fields']['simpleid']
    clinical_sites = "|".join(entry['fields']['seenatclinicalsites'])

    if 'phenotips' in entry['fields']\
            and entry['fields']['phenotips'] is not None\
            and 'features' in entry['fields']['phenotips']:
        phenotypes = entry['fields']['phenotips']['features']
        phenotype_list = "|".join([obj['id']for obj in phenotypes if obj['observed'] == 'yes'])
        solved = entry['fields']['phenotips']['solved']['status']

        sufficiency = entry['fields']['phenotips']['specificity']['score']
        pmid = ""
        if solved == 'solved' and 'pubmed_id' in entry['fields']['phenotips']['solved']:
            pmid = entry['fields']['phenotips']['solved']['pubmed_id']
        if solved == 'solved' and 'pubmed_id' in entry['fields']['phenotips']['solved']:
            pmid = entry['fields']['phenotips']['solved']['pubmed_id']
        count += 1
        output_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n"
                          .format(simple_id, phenotype_list,
                                  sufficiency, solved, count, pmid, clinical_sites))
