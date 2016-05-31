#!/usr/bin/env python3
import requests
import argparse
import re

"""
Prints gene id if it has a gene-disease association in monarch
Could theoretically do this with our solr index
"""

SCIGRAPH_BASE = "https://scigraph-data-dev.monarchinitiative.org/scigraph/cypher/execute.json"

parser = argparse.ArgumentParser(description='description')
parser.add_argument('--input', '-i', type=str, required=True,
                    help='Location of input file')
parser.add_argument('--output', '-o', type=str, required=True,
                    help='Location of output file')
args = parser.parse_args()

input_file = open(args.input, 'r')
output_file = open(args.output, 'w')
output_file.write("gene\tdisease\n")
disease_count = 0

for gene in input_file:
    curie = gene
    gene = gene.rstrip()
    gene_iri = re.sub('NCBIGene:', 'http://www.ncbi.nlm.nih.gov/gene/', gene)

    query = """
        START gene = node:node_auto_index(iri='{0}')
        MATCH (gene)<-[:GENO:0000408|GENO:0000639!]-(variant)-[:RO:0002200|RO:0002610|RO:0002326!]->(object:disease)
        RETURN DISTINCT
        gene.label as gene, object.label as disease
        UNION
        START gene = node:node_auto_index(iri='{0}')
        MATCH (gene)<-[:GENO:0000408|GENO:0000639!]-(variant)<-[:BFO:0000051!*]-(genotype:genotype)-[:RO:0002200|RO:0002610|RO:0002326!]->(object:disease)
        RETURN DISTINCT
        gene.label as gene, object.label as disease
        UNION
        START gene = node:node_auto_index(iri='{0}')
        MATCH (gene)<-[:GENO:0000408|GENO:0000639!]-(variant)<-[:BFO:0000051!*]-(genotype:genotype)<-[:GENO:0000222|RO:0001000*1..2]-(person)-[:RO:0002200!]->(object:disease)
        RETURN DISTINCT
        gene.label as gene, object.label as disease
    """.format(gene_iri)

    params = {'cypherQuery': query,
              'limit': 300}
    scigraph_request = requests.get(SCIGRAPH_BASE, params=params)
    response = scigraph_request.json()

    try:
        if len(list(response)) > 0:
            disease_count += 1
            print(gene)
            for result in list(response):
                output_file.write(
                    "{0}\t{1}\n".format(
                        gene, result['disease']))
    except TypeError:
        pass

print("Total gene-disease associations: {0}".format(disease_count))
