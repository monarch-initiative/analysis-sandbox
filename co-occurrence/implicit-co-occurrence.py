import itertools
import requests
from collections import Counter


solr = "https://solr-dev.monarchinitiative.org/solr/golr/select/"

closure_map = {}

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
    'fl': 'subject,object_closure',
    'fq': [
        'subject_category:"disease"',
        'object_category:"phenotype"',
        'relation_closure:"RO:0002200"'
     ]
}

for doc in get_solr_docs(solr, params):
    try:
        closures = [pheno[3:] for pheno in doc['object_closure'] if pheno.startswith("HP:")]
        closure_map[doc['subject']].extend(closures)
    except KeyError:
        closures = [pheno[3:] for pheno in doc['object_closure'] if pheno.startswith("HP:")]
        closure_map[doc['subject']] = closures

closure_list = [closures for closures in closure_map.values()]
closure_map = None

# Get co-occurring terms in a list of lists, credit
# https://stackoverflow.com/a/42272583

# generate combinations for each sub list seperately
lists_of_pairs = [list(itertools.combinations(sub_list, 2)) for sub_list in closure_list]
closure_list = None
# flatten the lists of pairs to 1 large list of pairs
all_pairs = [pair if pair[0] >= pair[1] else pair[::-1] for pairs_list in lists_of_pairs for pair in pairs_list]
lists_of_pairs = None

# let the Counter do the rest for you
co_occurrence_counts = Counter(all_pairs)

# Get a id-label map
solr_search = "https://solr-dev.monarchinitiative.org/solr/search/select/"

params = {
    'wt': 'json',
    'rows': 500,
    'start': 0,
    'q': '*:*',
    'fl': 'id,label',
    'fq': [
        'prefix:HP',
        'category:Phenotype',
     ]
}

label_map = dict()

for doc in get_solr_docs(solr_search, params):
   label_map[doc['id']] = doc['label'][0]
   try:
       for curie in doc['equivalent_curie']:
           if curie.startswith('HP'):
               label_map[curie] = doc['label'][0]
   except KeyError:
       pass

label_map['HP:0000001'] = 'All'

for k,v in co_occurrence_counts.items():
    curie1 = "HP:{}".format(k[0])
    curie2 = "HP:{}".format(k[1])
    try:
        label1 = label_map[curie1]
    except KeyError:
        label1 = ''
    try:
        label2 = label_map[curie2]
    except KeyError:
        label2 = ''
    print("{}\t{}\t{}\t{}\t{}".format(curie1, label1, curie2, label2, v))

