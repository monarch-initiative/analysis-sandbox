from monarch import monarch
import argparse
import logging
import re

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


parser = argparse.ArgumentParser(description='Generates annotation'
                                 'sufficiency scores for set of phenotype ids')
parser.add_argument('--input', '-i', type=str, required=True,
                    help='Location of input file')
parser.add_argument('--output', '-o', type=str, required=True,
                    help='Location of output file')

args = parser.parse_args()

input_file = open(args.input, 'r')
output_file = open(args.output, 'w')

for index, line in enumerate(input_file):
    fields = re.split(r'\t', line)
    phenotypes = fields[2].split("|")
    score = monarch.get_annotation_sufficiency_score(phenotypes)
    output_file.write("{0}\t{1}\t{2}\t{3}".format(fields[0], fields[1], score, fields[2]))
    if index % 100 == 0:
        logger.info("processed {0} profiles".format(index))

output_file.close()
input_file.close()