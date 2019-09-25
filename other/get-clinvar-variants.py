from ontobio.obograph_util import convert_json_object, obograph_to_assoc_results
import requests
import json

url = "https://solr-dev.monarchinitiative.org/solr/golr/select/?q=*%3A*&fq=association_type:%22marker_disease%22&rows=1000&wt=json&fq=is_defined_by:%22https://data.monarchinitiative.org/ttl/clinvar.nt%22&fq=-is_defined_by:%22https://data.monarchinitiative.org/ttl/omim.ttl%22&fl=subject,subject_label,relation_label,object,object_label,evidence_graph"

req = requests.get(url)
resp = req.json()

output = open('/home/kshefchek/clinvar-variants.tsv', 'w')
columns = ['gene', 'symbol', 'variant', 'hgvs', 'relation', 'mondo', 'mondo_label', 'disease_id', 'pubs']
output.write("\t".join(columns) + "\n")

for doc in resp['response']['docs']:
    ev = json.loads(doc['evidence_graph'])
    eg = {'graphs': [ev]}
    digraph = convert_json_object(eg, reverse_edges=False)['graph']
    association_results = obograph_to_assoc_results(digraph)
    for assoc in association_results:
        output.write("\t".join([
            doc['subject'],
            doc['subject_label'],
            assoc['subject']['id'],
            assoc['subject']['label'],
            doc['relation_label'],
            doc['object'],
            doc['object_label'],
            "|".join([foo for foo in assoc['object_eq'] if 'omim' in foo or 'ORDO' in  foo]),
            "|".join([foo['id'] for foo in assoc['publications']])
        ]) + "\n")
