import re
import argparse
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

"""
Compare too list files with the format
1   foo|bar|baz
2   foo|baz|qux

And determine if identical sets (ignoring dupes) exist
"""

parser = argparse.ArgumentParser(
    description='Checks for equality of lists in two files'
)
parser.add_argument('--reference', '-r', type=str, required=True,
                    help='Location of input file')
parser.add_argument('--query', '-q', type=str, required=True,
                    help='Location of output file')

args = parser.parse_args()

reference_file = open(args.reference, 'r')
query_file = open(args.query, 'r')
query_list = query_file.readlines()
phenotypes_list = []

for line in reference_file:
    line = line.rstrip('\n')
    fields = re.split(r'\t', line)
    phenotypes = fields[1].split("|")
    phenotypes_list.append(phenotypes)

for i, ref_phenotypes in enumerate(phenotypes_list):
    for line_num, line in enumerate(query_list):
        line = line.rstrip('\n')
        fields = re.split(r'\t', line)
        query_phenotypes = fields[1].split("|")
        set1 = set(ref_phenotypes)
        set2 = set(query_phenotypes)
        if set1 == set2:
            print("{0} line in query file same as ref line {1}".format(line_num, i))