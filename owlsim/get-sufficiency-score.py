from monarch import monarch
import argparse
import logging
import re

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


parser = argparse.ArgumentParser(description='Generates annotation'
                                 'sufficiency scores for set of phenotype ids')
parser.add_argument('--input', '-i', type=str, required=True,
                    help='Location of input file'
                         'Expected fields:'
                         '1. Disease id as curie'
                         '2. Disease label'
                         '3. Pipe separated phenotype IDs as curies'
                         '4. Pipe separated phenotype labels')
parser.add_argument('--output', '-o', type=str, required=True,
                    help='Location of output file')

args = parser.parse_args()

input_file = open(args.input, 'r')
output_file = open(args.output, 'w')

for index, line in enumerate(input_file):
    fields = re.split(r'\t', line)
    phenotypes = fields[2].split("|")
    score = monarch.get_annotation_sufficiency_score(phenotypes)
    output_file.write("{0}\t{1}\t{2}\t{3}\t{4}".format(fields[0], fields[1], score, fields[2], fields[3]))
    if index % 100 == 0:
        logger.info("processed {0} profiles".format(index))

output_file.close()
input_file.close()