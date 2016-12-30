import json
import argparse
import logging

import re
import monarch


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

"""
Generates similarity matrix based on OwlSim2 (PhenoDigm)
phenotypic similarity scores.

Input is a tab delimited file in the format:
foobar:123\tHP:123|HP:345|HP:567

Output is a json formatted asymmetric 2d array

Supports batch requests using grequests; however, when using
the monarch owlsim endpoint --chunk should not be set above 3

--temp /path/to/file generates a temporary file that can be used
as a cache for subsequent runs

TODOs
1. Get rid of distance matrix as this can be generated
   from the sim matrix in one line
2. Add owlsim3 support and support for symmetric algorithms
   (code for symmetric matrix in old commits)
"""

parser = argparse.ArgumentParser(description='Generates similarity'
                                 ' and distance matrices for'
                                 'pairwise phenotype comparison')
parser.add_argument('--input', '-i', type=str, required=True,
                    help='Location of input file')
parser.add_argument('--similarity', '-sim', type=str, required=True,
                    help='Location of output file')
parser.add_argument('--distance', '-dist', type=str, required=True,
                    help='Location of output file')
parser.add_argument('--temp', '-t', type=str, required=False,
                    help='Location of temp file')
parser.add_argument('--cache', '-c', type=str, required=False,
                    help='Location of cache file')
parser.add_argument('--chunk', '-ch', type=int, required=False,
                    default=20, help='Number of ids to send to server in one call')
args = parser.parse_args()

input_file = open(args.input, 'r')
similarity_file = open(args.similarity, 'w')
distance_file = open(args.distance, 'w')
if args.temp:
    temp_file = open(args.temp, 'w')
    temp_file.close()

sample_ids = []

for line in input_file:
    line = line.rstrip('\n')
    fields = re.split(r'\t', line)
    phenotypes = fields[1].split("|")
    sample_ids.append(phenotypes)

input_file.close()

if args.cache:
    cached_matrix = open(args.cache, 'r')
    similarity_matrix = json.load(cached_matrix)
    distance_matrix = [
        [
            100-value if value != 0 and index != i else 0
            for index, value in enumerate(similarity_matrix[i])
        ]
        for i in range(len(similarity_matrix))
    ]
    cached_matrix.close()
else:
    # This could be simplified if we use numpy
    similarity_matrix = [[0 for k in range(len(sample_ids))] for i in range(len(sample_ids))]
    distance_matrix = [[0 for k in range(len(sample_ids))] for i in range(len(sample_ids))]

for index, value in enumerate(sample_ids):
    for chunk_num, query_list in \
            enumerate([sample_ids[i:i + args.chunk] for i in range(0, len(sample_ids), args.chunk)]):
        # Determine where we are starting on the x axis of the matrix
        x_axis_index = chunk_num * args.chunk

        end_index = (x_axis_index + args.chunk) - 1
        if end_index > len(sample_ids):
            end_index = len(sample_ids) - 1

        if args.cache:
            is_matrix_filled = True
            if sum(1 for i in range(x_axis_index, (end_index + 1)) if (similarity_matrix[index][i] == 0
                   and distance_matrix[index][i] != 100)
                   or (similarity_matrix[index][i] == 100 and index != i)) > 0:
                is_matrix_filled = False
        else:
            is_matrix_filled = False

        if not is_matrix_filled:
            try:
                scores = monarch.compare_attribute_sets(value, query_list)
                for score_index, score in enumerate(scores):
                    if index != x_axis_index and score == 100:
                        # Similarity of 100 is unexpected
                        logger.warn("Errant score at {0} {1}".format(index, x_axis_index))
                    similarity_matrix[index][x_axis_index] = score
                    distance_matrix[index][x_axis_index] = 100 - score
                    x_axis_index += 1
            except ConnectionError:
                for index_query in range(x_axis_index, (end_index + 1)):
                    similarity_matrix[index][index_query] = 0
                    distance_matrix[index][index_query] = 0

    # Dump matrix to temp file every now and then
    if index % 50 == 0:
        logger.info("Processed {0} samples".format(index))
        if args.temp:
            temp_file = open(args.temp, 'w')
            temp_file.write(json.dumps(similarity_matrix))
            temp_file.close()

similarity_file.write(json.dumps(similarity_matrix))
distance_file.write(json.dumps(distance_matrix))
