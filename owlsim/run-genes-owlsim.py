#!/usr/bin/env python3

import requests
import argparse
import csv
import re
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)
logging.getLogger('urllib3').setLevel(logging.ERROR)

# Custom session for services that require more than 3 retries
session = requests.Session()
adapter = requests.adapters.HTTPAdapter(max_retries=10)
session.mount('https://', adapter)

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
    parser.add_argument('--output', '-o', type=str, required=True,
                        help='Location of output file')
    args = parser.parse_args()

    output_file = open(args.output, 'w')

    logger.info("Processing input file")
    disease_dictionary = process_input_file(args.input)
    logger.info("Replacing NCBIGene with ZFIN")
    disease_dictionary = get_zfin_ids(disease_dictionary)
    logger.info("Getting disease leaf/group info")
    disease_dictionary = get_disease_info(disease_dictionary)
    logger.info("Getting disease-phenotype, gene-phenotype, and disease-gene counts from solr")
    disease_dictionary = get_solr_counts(disease_dictionary)
    logger.info("Getting owlsim score/rank")
    disease_dictionary = get_owlsim_scores(disease_dictionary)

    output_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\n".format(
                      "disease", "disease_label", "model_gene", "gene_label", "human_gene", "model_gene_taxon",
                      "owlsim_score", "owlsim_rank", "models", "disease_pheno_count",
                      "disease_gene_count", "modgene_pheno_count", "isLeafNode", "absolute_rank", "human_gene_label"
                      ))

    for disease_id, disease in disease_dictionary.items():
        output_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\n".format(
                          disease["disease"], disease["disease_label"],
                          disease["model_gene"], disease["model_gene_label"],
                          disease["human_gene"], TAXON_MAP[disease["model_gene_taxon"]],
                          disease["owlsim_score"], disease["owlsim_rank"],
                          disease["models"], disease["disease_pheno_count"],
                          disease["disease_gene_count"], disease["modgene_pheno_count"],
                          disease["isLeafNode"], disease["abs_rank"], disease["human_gene_label"]
        ))


def get_owlsim_scores(disease_dictionary):

    # Get all phenotypes
    for disease_id, disease in disease_dictionary.items():
        params = {
            'wt': 'json',
            'rows': 100,
            'start': 0,
            'q': '*:*',
            'fl': 'object',
            'fq': ['subject_closure:"{0}"'.format(disease["disease"]),
                   'object_category:"phenotype"']
        }

        resultCount = params['rows']
        phenotype_list = []

        while params['start'] < resultCount:
            solr_request = requests.get(SOLR_URL, params=params)
            response = solr_request.json()
            temp_list = [doc['object'] for doc in response['response']['docs']]
            phenotype_list = phenotype_list + temp_list
            params['start'] += params['rows']

        phenotypes = "+".join(phenotype_list)
        taxon = re.sub(r"NCBITaxon:", "", disease["model_gene_taxon"])

        url = "{0}?input_items={1}&target_species={2}".format(OWLSIM_URL, phenotypes, taxon)
        # Run through owlsim through monarch simsearch endpoint
        owlsim_request = requests.post(url)
        owlsim_results = owlsim_request.json()
        rank = 0
        last_score = -1
        if "b" in owlsim_results:
            is_found = False
            for result in owlsim_results["b"]:
                if result["score"]["score"] != last_score:
                    rank += 1
                last_score = result["score"]["score"]
                if result["id"] == disease["model_gene"]:
                    disease["owlsim_score"] = result["score"]["score"]
                    disease["owlsim_rank"] = rank
                    disease["abs_rank"] = result["score"]["rank"] + 1
                    is_found = True
            if not is_found:
                logger.warn("No owlsim results found for {0}"
                            " in disease {1}".format(disease["model_gene"],
                                                     disease["disease"]))
                disease["owlsim_score"] = ""
                disease["owlsim_rank"] = ""
                disease["abs_rank"] = ""
        else:
            logger.warn("No owlsim results found for {0}".format(disease["disease"]))

    return disease_dictionary


def process_input_file(input_file):
    """
    Convert input tsv with fields disease, model, model gene, human gene
    to
    {
        'OMIM:123-modGene:123': {
            disease:
            model_gene:
            human_gene:
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

            # Get label, next time get this from the cypher query
            disease_label = get_label_from_scigraph(disease)
            human_gene_label = get_label_from_scigraph(human_gene)
            model_gene_label = get_label_from_scigraph(model_gene)

            if gene_taxon == "NCBITaxon:10116" and model.startswith("MGI"):
                # Catch and skip over rat transgenes
                continue

            if key not in disease_dictionary:
                disease_dictionary[key] = {
                    "disease": disease,
                    "disease_label": disease_label,
                    "human_gene": human_gene,
                    "human_gene_label": human_gene_label,
                    "model_gene": model_gene,
                    "model_gene_label": model_gene_label,
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
        disease_gene_filter = \
            ['subject_closure:"{0}"'.format(disease["model_gene"]),
             'object_category:"phenotype"']

        gene_pheno_filter = \
            ['object_closure:"{0}"'.format(disease["disease"]),
             'subject_category:"gene"']
        disease["disease_pheno_count"] = get_solr_result_count(disease_pheno_filter)
        disease["disease_gene_count"] = get_solr_result_count(disease_gene_filter)
        disease["modgene_pheno_count"] = get_solr_result_count(gene_pheno_filter)

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
                disease["isLeafNode"] = "true"
            else:
                disease["isLeafNode"] = "false"
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
        if disease["model_gene_taxon"] == "NCBITaxon:7955" \
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
                    disease["model_gene"] = node["id"]
                    node_set.add(node["id"])

            if len(node_set) > 1:
                logger.error("More than one equivalent zfin id found for {0}".format(gene))

    return disease_dictionary


def get_label_from_scigraph(curie):
    scigraph_service = SCIGRAPH_URL + "/graph/" + curie + ".json"
    request = requests.get(scigraph_service)
    results = request.json()
    label = results["nodes"][0]["lbl"]

    return label


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
    request = session.get(scigraph_service, params=params)
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
