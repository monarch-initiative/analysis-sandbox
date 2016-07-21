#!/usr/bin/env python3

import requests
import argparse
import csv
import re
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)
logging.getLogger('urllib3').setLevel(logging.ERROR)

# Globals and Constants
SCIGRAPH_URL = 'https://scigraph-data.monarchinitiative.org/scigraph'
OWLSIM_URL = 'https://monarchinitiative.org/simsearch/phenotype'
SOLR_URL = 'https://solr.monarchinitiative.org/solr/golr/select'

CURIE_MAP = {
    "http://purl.obolibrary.org/obo/OMIM_": "OMIM",
    "http://purl.obolibrary.org/obo/DOID_": "DOID",
    "http://www.informatics.jax.org/accession/MGI:": "MGI",
    "http://www.ncbi.nlm.nih.gov/gene/": "NCBIGene",
    "http://zfin.org/": "ZFIN",
    "http://purl.obolibrary.org/obo/NCBITaxon_": "NCBITaxon"
}

TAXON_MAP = {
    "NCBITaxon:9606": "Human",
    "NCBITaxon:10090": "Mouse",
    "NCBITaxon:7955": "Zebrafish",
    "NCBITaxon:10116": "Rat",
}


def main():

    parser = argparse.ArgumentParser(description='description')
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='Location of input file')
    args = parser.parse_args()

    logger.info("Processing input file")
    disease_dictionary = process_input_file(args.input)
    logger.info("Replacing NCBIGene with ZFIN")
    disease_dictionary = get_zfin_ids(disease_dictionary)
    logger.info("Getting disease leaf/group info")
    disease_dictionary = get_disease_info(disease_dictionary)
    logger.info("Getting disease-phenotype and gene-phenotype counts from solr")
    disease_dictionary = get_solr_counts(disease_dictionary)
    logger.info("Getting owlsim score/rank")

    print(disease_dictionary)


def process_input_file(input_file):
    """
    Convert input tsv with fields disease, model, model gene, human gene
    to
    {
        'OMIM:123-modGene:123': {
            disease:
            model_gene:
            marker_gene:
            models: 3

    }
    :param input_file:
    :return:
    """
    disease_dictionary = dict()
    with open(input_file, 'rt') as tsvfile:
        reader = csv.reader(tsvfile, delimiter="\t")
        next(reader, None)
        for row in reader:
            (disease, model, model_gene, human_gene) = [map_iri_to_curie(val) for val in row]

            key = "{0}-{1}".format(disease, model_gene)
            gene_taxon = get_taxon(model_gene)

            if gene_taxon == "NCBITaxon:10116" and model.startswith("MGI"):
                # Catch and skip over rat transgenes
                continue

            if key not in disease_dictionary:
                disease_dictionary[key] = {
                    "disease": disease,
                    "marker_gene": human_gene,
                    "model_gene": model_gene,
                    "models": 1,
                    "model_gene_taxon": gene_taxon
                }
            else:
                disease_dictionary[key]['models'] += 1

    return disease_dictionary


def get_solr_counts(disease_dictionary):
    for disease_id, disease in disease_dictionary.items():
        disease_pheno_filter = \
            ['subject_closure:"{0}"'.format(disease["disease"]),
             'object_category:"phenotype"']
        gene_pheno_filter = \
            ['subject_closure:"{0}"'.format(disease["model_gene"]),
             'object_category:"phenotype"']
        disease["disease_pheno"] = get_solr_result_count(disease_pheno_filter)
        disease["model_gene_pheno"] = get_solr_result_count(gene_pheno_filter)

    return disease_dictionary


def get_solr_result_count(filters):
    """
    Get result Count from solr
    :param filters: list
    :return:
    """
    params = {
        'wt': 'json',
        'rows': 0,
        'start': 0,
        'q': '*:*',
        'fq': filters
    }
    solr_request = requests.get(SOLR_URL, params=params)
    response = solr_request.json()
    result_count = response['response']['numFound']
    return result_count


def get_disease_info(disease_dictionary):
    scigraph_service = SCIGRAPH_URL + "/graph/neighbors"
    for disease_id, disease in disease_dictionary.items():
            params = {
               "id": disease["disease"],
               "relationshipType": "subClassOf",
               "direction": "INCOMING",
               "depth": 2
            }

            request = requests.get(scigraph_service, params=params)
            results = request.json()

            if len(results["edges"]) == 0:
                disease["isLeafNode"] = "false"
            else:
                disease["isLeafNode"] = "true"
    return disease_dictionary


def get_zfin_ids(disease_dictionary):
    """
    Owlsim can be limited by ID prefix, but some ZFIN genes show
    up with NCBIGene prefixes.  This function converts them to
    their zfin equivalents

    :param disease_dictionary:
    :return:
    """
    eq_predicate = "equivalentClass"
    scigraph_service = SCIGRAPH_URL + "/graph/neighbors"
    for disease_id, disease in disease_dictionary.items():
        if disease["model_gene_taxon"] == "NCBIGene:7955" \
                and not disease["model_gene"].startswith("ZFIN"):

            gene = disease["model_gene"]

            params = {
               "id": gene,
               "relationshipType": eq_predicate,
               "direction": "BOTH",
               "depth": 3
            }

            request = requests.get(scigraph_service, params=params)
            results = request.json()
            node_set = set()
            for node in results["nodes"]:
                if node["id"].startswith("ZFIN"):
                    node_set.add(node["id"])

            if len(node_set) > 1:
                logger.error("More than one equivalent zfin id found for {0}".format(gene))
            else:
                disease["model_gene"] = node_set[0]

    return disease_dictionary


def map_iri_to_curie(iri):
    curie = iri
    for prefix in CURIE_MAP:
        if iri.startswith(prefix):
            curie = re.sub(r"{0}".format(prefix),
                           "{0}:".format(CURIE_MAP[prefix]), iri)
            break
    return curie


def get_taxon(curie):
    taxon = None
    taxon_predicate = "http://purl.obolibrary.org/obo/RO_0002162"
    scigraph_service = SCIGRAPH_URL + "/graph/neighbors"
    params = {
        "id": curie,
        "relationshipType": taxon_predicate,
        "direction": "OUTGOING"
    }
    request = requests.get(scigraph_service, params=params)
    results = request.json()

    if len(results["edges"]) == 1:
        taxon = results["edges"][0]["obj"]
    elif len(results["edges"]) < 1:
        logger.warn("No taxon found for {0}".format(curie))
    elif len(results["edges"]) > 1:
        taxon = results["edges"][0]["obj"]
        taxon_set = set()
        for edge in results["edges"]:
            taxon_set.add(edge["obj"])
        if len(taxon_set) > 1:
            logger.warn("More than one taxon found for {0}".format(curie))

    return taxon


if __name__=="__main__":
    main()
