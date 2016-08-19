#!/usr/bin/env python3
import requests
import argparse
import re

SCIGRAPH_BASE = "https://scigraph-data.monarchinitiative.org/scigraph/cypher/execute.json"

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
taxon_stats[5] = 0


for gene in input_file:
    gene = gene.rstrip()
    gene_iri = re.sub('NCBIGene:', 'http://www.ncbi.nlm.nih.gov/gene/', gene)

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
        UNION
        START gene = node:node_auto_index(iri='{0}')
        MATCH (gene)-[rel:RO:HOM0000017|RO:HOM0000020]-(ortholog)<-[:GENO:0000408]-(allele)-[:BFO:0000051!*]->(feature)-[:RO:0002610!]->(object:Phenotype)
        MATCH (ortholog)-[tax:RO:0002162]->(taxon)
        RETURN DISTINCT
        ortholog.label as ortholog, gene.label as gene, object.label as phenotype , taxon.label as taxon,
        'pattern4' as pattern
    """.format(gene_iri)

    # Catch case where gene label is null
    if re.match(r'NCBIGene:(7145|3678|23326|158880)', gene):
        query = re.sub(r'ortholog\.label', 'ortholog.iri', query)

    params = {'cypherQuery': query,
              'limit': 300}
    scigraph_request = requests.get(SCIGRAPH_BASE, params)
    response = scigraph_request.json()

    taxes = set()
    try:
        if len(list(response)) > 0:
            if 'code' in response:
                if response['code'] == 500:
                    print("Code 500 - Likely error in query")
                    print(query)
                    print(scigraph_request.url)
                    continue
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
print("Number of genes with ortholog-phenotypes from 1 taxon: {0}".format(taxon_stats[1]))
print("Number of genes with ortholog-phenotypes from 2 taxa: {0}".format(taxon_stats[2]))
print("Number of genes with ortholog-phenotypes from 3 taxa: {0}".format(taxon_stats[3]))
print("Number of genes with ortholog-phenotypes from 4 taxa: {0}".format(taxon_stats[4]))
print("Number of genes with ortholog-phenotypes from 5 taxa: {0}".format(taxon_stats[5]))

