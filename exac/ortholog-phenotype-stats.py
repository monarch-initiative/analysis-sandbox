#!/usr/bin/env python3
import requests
import argparse
import re

SCIGRAPH_BASE = "https://scigraph-data-dev.monarchinitiative.org/scigraph/cypher/execute.json"

parser = argparse.ArgumentParser(description='description')
parser.add_argument('--input', '-i', type=str, required=True,
                    help='Location of input file')
parser.add_argument('--output', '-o', type=str, required=True,
                    help='Location of output file')
args = parser.parse_args()

input_file = open(args.input, 'r')
output_file = open(args.output, 'w')
output_file.write("ortholog\tgene\tphenotype\ttaxon\n")
phenotype_count = 0

taxon_stats = {}
taxon_stats[1] = 0
taxon_stats[2] = 0
taxon_stats[3] = 0
taxon_stats[4] = 0

for gene in input_file:
    gene = gene.rstrip()
    gene = re.sub('NCBIGene:', 'http://www.ncbi.nlm.nih.gov/gene/', gene)

    query = """
        START gene = node:node_auto_index(iri='{0}')
        MATCH (gene)-[rel:RO:HOM0000017|RO:HOM0000020]-(ortholog)<-[:GENO:0000408!*0..1]-(feature)-[p:RO:0002200|RO:0002326!*]->(object:Phenotype)
        MATCH (ortholog)-[tax:RO:0002162]->(taxon)
        RETURN DISTINCT
        ortholog.label as ortholog, gene.label as gene, object.label as phenotype , taxon.label as taxon,
        'pattern1' as pattern
        UNION
        START gene = node:node_auto_index(iri='{0}')
        MATCH (gene)-[rel:RO:HOM0000017|RO:HOM0000020]-(ortholog)<-[:GENO:0000408|GENO:0000639!*0..1]-(feature)<-[:BFO:0000051!*]-(genotype:genotype)-[p:RO:0002200|RO:0002326!*]->(object:Phenotype)
        MATCH (ortholog)-[tax:RO:0002162]->(taxon)
        RETURN DISTINCT
        ortholog.label as ortholog, gene.label as gene, object.label as phenotype , taxon.label as taxon,
        'pattern2' as pattern
        UNION
        START gene = node:node_auto_index(iri='{0}')
        MATCH (gene)-[rel:RO:HOM0000017|RO:HOM0000020]-(ortholog)<-[:GENO:0000408|GENO:0000639!*0..1]-(feature)<-[:BFO:0000051!*]-(genotype:genotype)<-[:GENO:0000222|RO:0001000*1..2]-(person)-[p:RO:0002200|RO:0002326!*]->(object:Phenotype)
        MATCH (ortholog)-[tax:RO:0002162]->(taxon)
        RETURN DISTINCT
        ortholog.label as ortholog, gene.label as gene, object.label as phenotype , taxon.label as taxon,
        'pattern3' as pattern
    """.format(gene)

    params = {'cypherQuery': query,
              'limit': 300}
    scigraph_request = requests.get(SCIGRAPH_BASE, params)
    response = scigraph_request.json()

    taxes = set()
    try:
        if len(list(response)) > 0:
            phenotype_count += 1
            for result in list(response):
                output_file.write(
                    "{0}\t{1}\t{2}\t{3}\n".format(
                        result['ortholog'], result['gene'],
                        result['phenotype'], result['taxon']))
                taxes.add(result['taxon'])
            tax_len = len(taxes)
            taxon_stats[tax_len] += 1

    except TypeError:
        pass

print("Total orthologs with phenotypes: {0}".format(phenotype_count))
print("Number of genes with othorlog-phenotypes from 1 taxon: {0}".format(taxon_stats[1]))
print("Number of genes with othorlog-phenotypes from 2 taxa: {0}".format(taxon_stats[2]))
print("Number of genes with othorlog-phenotypes from 3 taxa: {0}".format(taxon_stats[3]))
print("Number of genes with othorlog-phenotypes from 4 taxa: {0}".format(taxon_stats[4]))
