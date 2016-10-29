import json
import argparse
import copy
import logging
from monarch import monarch

logger = logging.getLogger(__name__)


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
    sample_ids.append(line.rstrip('\n'))

sample_tmp = copy.deepcopy(sample_ids)

if args.cache:
    cached_matrix = open(args.cache, 'r')
    similarity_matrix = json.load(cached_matrix)
    distance_matrix = [[100-k if k != 0 else 0 for k in similarity_matrix[i]] for i in range(len(similarity_matrix))]
else:
    similarity_matrix = [[0 for k in range(len(sample_ids))] for i in range(len(sample_ids))]
    distance_matrix = [[0 for k in range(len(sample_ids))] for i in range(len(sample_ids))]


for index, value in enumerate(sample_ids):
    offset = len(sample_ids) - len(sample_tmp)
    for chunk_num, query_list in \
            enumerate([sample_tmp[i:i + args.chunk] for i in range(0, len(sample_tmp), args.chunk)]):
        # Determine where we are starting on the x axis of the matrix
        x_axis_index = chunk_num * args.chunk + offset

        end_index = (x_axis_index + args.chunk) - 1 + offset
        if end_index > len(sample_ids):
            end_index = len(sample_ids)

        if args.cache:
            is_matrix_filled = True
            for query_index in range(x_axis_index, end_index):
                if similarity_matrix[index][query_index] is 0 \
                        and distance_matrix[index][query_index] is not 100:
                    is_matrix_filled = False
                    break
        else:
            is_matrix_filled = False

        if not is_matrix_filled:
            try:
                scores = monarch.get_score_from_compare(value, query_list)
            except ConnectionError:
                similarity_score = [0 for i in range(chunk_num)]
                distance_score = [0 for i in range(chunk_num)]

            for score_index, score in enumerate(scores):
                index_query = x_axis_index
                similarity_matrix[index][index_query] = score
                similarity_matrix[index_query][index] = score
                distance_matrix[index][index_query] = 100 - score
                distance_matrix[index_query][index] = 100 - score
                x_axis_index += 1

    sample_tmp.pop(0)

    # Dump matrix to temp file every now and then
    if index % 50 == 0:
        logger.info("Processed {0} samples".format(index))
        if args.temp:
            temp_file = open(args.temp, 'w')
            temp_file.write(json.dumps(similarity_matrix))
            temp_file.close()


similarity_file.write(json.dumps(similarity_matrix))
distance_file.write(json.dumps(distance_matrix))

