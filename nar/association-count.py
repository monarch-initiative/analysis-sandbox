import requests

SOLR_URL = 'https://solr.monarchinitiative.org/solr/golr/select'

params = {
    'wt': 'json',
    'rows': 10000,
    'start': 0,
    'fq': '!subject_category:"publication"',
    'q': '*:*',
    'fl': 'subject,object,is_defined_by,relation'
}

resultCount = params['rows']
associations = set()
source_count = 0

while params['start'] < resultCount:
    solr_request = requests.get(SOLR_URL, params)
    response = solr_request.json()
    resultCount = response['response']['numFound']
    for doc in response['response']['docs']:
        if 'relation' in doc:
            subj_obj = "{0}-{1}-{2}".format(doc['subject'], doc['relation'], doc['object'])
        else:
            subj_obj = "{0}-{1}".format(doc['subject'], doc['object'])
        if subj_obj not in associations:
            associations.add(subj_obj)
            if 'is_defined_by' in doc and len(doc['is_defined_by']) > 1:
                source_count += 1
    if params['start'] % 100000 == 0:
        print("processed {0} out of {1}".format(params['start'], resultCount))
    params['start'] += params['rows']

print("Number of associations: {0}".format(len(associations)))
print("Number of associations with >1 source: {0}".format(source_count))
