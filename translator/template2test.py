"""
Automatically fill out a testing template using the neo4j HTTP endpoint
"""

import requests
import json

sri_graph = 'http://scigraph.ncats.io/db/data/transaction'

template_file = './onehop/templates/KP/Standards_Reference_Implementation_Team/' \
                'SRI_Reference_Knowledge_Graph_API.json'

with open(template_file, 'r') as template_fh:
    template = json.load(template_fh)

output_json = {
    'url': template['url'],
    'TRAPI': template['TRAPI'],
    'edges': []
}

for triple in template['edges']:
    cypher_query = \
        f"MATCH (subject:`{triple['subject_category']}`)" \
        f"-[:`{triple['predicate']}`]" \
        f"->(object:`{triple['object_category']}`) " \
        f"RETURN subject, object " \
        f"LIMIT 1"
    statements = {
        'statements': [
            {
                'statement': cypher_query
            }
        ]
    }
    req = requests.post(sri_graph, data=json.dumps(statements))

    try:
        row = req.json()['results'][0]['data'][0]['row']
        triple['subject'] = row[0]['id']
        triple['object'] = row[1]['id']
        output_json['edges'].append(triple)
    except Exception:
        continue

with open('./output.json', 'w') as output:
    output.write(json.dumps(output_json))
