import itertools
import requests
from collections import Counter

solr = "https://solr-dev.monarchinitiative.org/solr/golr/select/"

closure_list = []


def get_solr_docs(solr, params):
    resultCount = params['rows']
    while params['start'] < resultCount:
        solr_request = requests.get(solr, params=params)
        response = solr_request.json()
        resultCount = response['response']['numFound']
        params['start'] += params['rows']
        for doc in response['response']['docs']:
            yield doc

params = {
    'wt': 'json',
    'rows': 500,
    'start': 0,
    'q': '*:*',
    'fl': 'object_closure_label',
    'fq': [
        'subject_category:"disease"',
        'object_category:"phenotype"',
        'relation_closure:"RO:0002200"'
     ]
}

for doc in get_solr_docs(solr, params):
    closure_list.append(doc['object_closure_label'])

# Credit https://stackoverflow.com/a/49079618

# generate combinations for each sub list seperately
lists_of_pairs = [list(itertools.combinations(sub_list, 2)) for sub_list in closure_list]
# flatten the lists of pairs to 1 large list of pairs
all_pairs = [pair for pairs_list in lists_of_pairs for pair in pairs_list]
# let the Counter do the rest for you
co_occurrence_counts = Counter(all_pairs)

for k,v in co_occurrence_counts.items():
    print("{}\t{}\t{}".format(k[0], k[1], v))
