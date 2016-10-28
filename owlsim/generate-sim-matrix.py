import json
import argparse
from monarch import monarch

parser = argparse.ArgumentParser(description='description')
parser.add_argument('--input', '-i', type=str, required=True,
                    help='Location of input file')
parser.add_argument('--similarity', '-sim', type=str, required=True,
                    help='Location of output file')
parser.add_argument('--distance', '-dist', type=str, required=True,
                    help='Location of output file')
args = parser.parse_args()

input_file = open(args.input, 'r')
similarity_file = open(args.similarity, 'w')
distance_file = open(args.distance, 'w')

sample_ids = []

for line in input_file:
    sample_ids.append(line.rstrip('\n'))

similarity_matrix = [[0 for i in range(len(sample_ids))] for i in range(len(sample_ids))]
distance_matrix = [[0 for i in range(len(sample_ids))] for i in range(len(sample_ids))]


for index, value in enumerate(sample_ids):
    for index_query, value_query in enumerate(sample_ids):
        if similarity_matrix[index][index_query] is 0:
            score = monarch.get_score_from_compare(value, value_query)
            try:
                similarity_score = int(score)
            except ValueError:
                similarity_score = 0
            distance_score = 100 - similarity_score
            similarity_matrix[index][index_query] = similarity_score
            similarity_matrix[index_query][index] = similarity_score
            distance_matrix[index][index_query] = distance_score
            distance_matrix[index_query][index] = distance_score


similarity_file.write(json.dumps(similarity_matrix))
distance_file.write(json.dumps(distance_matrix))

