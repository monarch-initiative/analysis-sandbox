import requests
import re
import logging
import json

logger = logging.getLogger(__name__)

requests_log = logging.getLogger("requests.packages.urllib3")
requests_log.setLevel(logging.ERROR)

# Custom session for services that require more than 3 retries
session = requests.Session()
adapter = requests.adapters.HTTPAdapter(max_retries=10)
session.mount('https://', adapter)

# Globals and Constants
SCIGRAPH_URL = 'https://scigraph-data.monarchinitiative.org/scigraph'
OWLSIM_URL = 'https://monarchinitiative.org/simsearch/phenotype'
OWLSIM_COMPARE = 'https://monarchinitiative.org/compare'
MONARCH_SCORE = 'https://monarchinitiative.org/score'
SOLR_URL = 'https://solr.monarchinitiative.org/solr/golr/select'

CURIE_MAP = {
    "http://purl.obolibrary.org/obo/OMIM_": "OMIM",
    "http://purl.obolibrary.org/obo/DOID_": "DOID",
    "http://www.informatics.jax.org/accession/MGI:": "MGI",
    "http://www.ncbi.nlm.nih.gov/gene/": "NCBIGene",
    "http://zfin.org/": "ZFIN",
    "http://purl.obolibrary.org/obo/NCBITaxon_": "NCBITaxon",
    "https://catalog.coriell.org/0/Sections/Search/Sample_Detail.aspx?Ref=": "Coriell"
}


TAXON_MAP = {
    "NCBITaxon:9606": "Human",
    "NCBITaxon:10090": "Mouse",
    "NCBITaxon:7955": "Zebrafish",
    "NCBITaxon:10116": "Rat",
}


def get_owlsim_scores(disease_dictionary):

    # Get all phenotypes
    for disease_id, disease in disease_dictionary.items():
        disease["owlsim_rank"] = ""
        disease["owlsim_score"] = ""

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
            solr_request = session.get(SOLR_URL, params=params)
            response = solr_request.json()
            resultCount = response['response']['numFound']
            temp_list = [doc['object'] for doc in response['response']['docs']]
            phenotype_list.extend(temp_list)
            params['start'] += params['rows']

        phenotypes = "+".join(phenotype_list)
        if disease["model_gene_taxon"] is not None:
            taxon = re.sub(r"NCBITaxon:", "", disease["model_gene_taxon"])
            params = "input_items={0}&target_species={1}".format(phenotypes, taxon)
        else:
            params = "input_items={0}".format(phenotypes)
        # Run through owlsim through monarch simsearch endpoint
        try:
            owlsim_request = session.post(OWLSIM_URL, params=params)
        except requests.exceptions.ConnectionError:
            logger.warn("Error hitting owlsim for {0}".format(disease["disease"]))
            disease["owlsim_score"] = get_score_from_compare(disease["disease"],
                                                             disease["model_gene"])
            continue

        try:
            owlsim_results = owlsim_request.json()
        except ValueError:
            logger.warn("Error parsing json for {0}".format(disease["disease"]))
            disease["owlsim_rank"] = ""
            disease["owlsim_score"] = get_score_from_compare(disease["disease"],
                                                             disease["model_gene"])
            continue
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
                    is_found = True
            if not is_found:
                logger.warn("No owlsim results found for {0}"
                            " in disease {1}".format(disease["model_gene"],
                                                     disease["disease"]))
                disease["owlsim_score"] = get_score_from_compare(disease["disease"],
                                                                 disease["model_gene"])
        else:
            logger.warn("No owlsim results found for {0}".format(disease["disease"]))

    return disease_dictionary


def get_score_from_compare(reference, query):

    query_ids = ",".join(query)
    compare_url = OWLSIM_COMPARE + "/{0}/{1}.json".format(reference, query_ids)
    results = []
    try:
        owlsim_request = session.get(compare_url)
    except requests.exceptions.ConnectionError:
        logger.warn("Connection error fetch owlsim compare for gene {0}"
                    " in disease {1}".format(query, reference))
        raise ConnectionError

    try:
        owlsim_results = owlsim_request.json()
        if "b" not in owlsim_results:
            logger.warn("No owlsim compare results found for {0}"
                        " in disease {1}".format(query, reference))
            results = [0 for i in range(len(query))]
        else:
            for query_id in query:
                result = [result["score"]["score"] for result in owlsim_results["b"] if result["id"] == query_id]
                if len(result) == 1:
                    results.append(result[0])
                else:
                    logger.warn("No owlsim compare results found for {0}"
                                " in disease {1}".format(query_id, reference))
                    score = 0
                    results.append(score)

    except ValueError:
        logger.warn("Error parsing json for {0} and {1} for request {2}".format(reference, query, owlsim_request))

    return results


def get_annotation_sufficiency_score(id_list):
    phenotype_dictionary = dict()
    score = dict()

    phenotype_dictionary["features"] = list()
    for hp_id in id_list:
        phenotype_dictionary["features"].append({
            "id": hp_id,
            "label": "",
            "observed": "positive",
            "isPresent": "true"
        })
    phenotypes = json.dumps(phenotype_dictionary)

    params = {
        'annotation_profile': phenotypes
    }

    score_request = session.post(MONARCH_SCORE, data=params)

    response = score_request.json()
    score['simple_score'] = response['simple_score']
    score['scaled_score'] = response['scaled_score']

    return score


def get_solr_counts(disease_dictionary):
    for disease_id, disease in disease_dictionary.items():
        disease_pheno_filter = \
            ['subject_closure:"{0}"'.format(disease["disease"]),
             'object_category:"phenotype"']
        gene_pheno_filter = \
            ['subject_closure:"{0}"'.format(disease["model_gene"]),
             'object_category:"phenotype"']

        disease_gene_filter = \
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
    solr_request = session.get(SOLR_URL, params=params)
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

            request = session.get(scigraph_service, params=params)
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
               "depth": 1
            }

            request = session.get(scigraph_service, params=params)
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
    request = session.get(scigraph_service)
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