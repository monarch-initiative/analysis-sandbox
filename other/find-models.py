#!/usr/bin/env python3
import requests

# requires python 3.5
SOLR_BASE = "https://solr.monarchinitiative.org/solr/golr/select"

default_params = {
    'wt': 'json',
    'rows': 100,
    'start': 0,
    'q': '*:*'
}

# Get all diseases with phenotype closure CL:0000540PHENOTYPE
pheno_disease_params = {
    'fq': ['object_closure:"CL:0000540PHENOTYPE"', 'subject_category:"disease"']
}

params = {k: v for d in [default_params, pheno_disease_params] for k, v in d.items()}
resultCount = params['rows']
disease_set = set()

while params['start'] < resultCount:
    solr_request = requests.get(SOLR_BASE, params)
    response = solr_request.json()
    resultCount = response['response']['numFound']
    temp_set = {doc['subject'] for doc in response['response']['docs']}
    disease_set = (disease_set | temp_set)
    params['start'] += params['rows']


no_gene_set = set()
for disease in disease_set:
    disease_gene_params = {
        'fq': ['object_closure:"{0}"'.format(disease), 'subject_category:"gene"']
    }
    params = {k: v for d in [default_params, disease_gene_params] for k, v in d.items()}
    params['rows'] = 0
    solr_request = requests.get(SOLR_BASE, params)
    response = solr_request.json()
    if response['response']['numFound'] == 0:
        no_gene_set.add(disease)

print("\n".join(disease for disease in no_gene_set))


