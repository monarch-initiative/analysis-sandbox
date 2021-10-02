#! /usr/bin/env python3

"""

TODO move this to https://github.com/monarch-initiative/release-utils/
The purpose of this script to determine what percentage of human protein coding genes
have associated phenotype data in either humans or in orthologous sequences in other
species

Some caveats:
1. There are a 1 to many mapping between a human gene and its orthologs
2. The stats are somewhat misleading in that we are showing coverage for predominately
    loss of function variants of a gene
3. This does not include evidence of functional conservation across species

Where there are cases of many orthologs in a model organism to one human gene, these are
collapsed into 1 gene in the output count
"""
import requests
import re
from json import JSONDecodeError
import urllib.request


# Globals and Constants
SCIGRAPH_URL = "https://scigraph-data.monarchinitiative.org/scigraph"
SOLR_URL = "https://solr.monarchinitiative.org/solr/golr/select"


CURIE_MAP = {
    "http://identifiers.org/hgnc/HGNC:": "HGNC",
    "http://purl.obolibrary.org/obo/NCBITaxon_": "NCBITaxon",
}


TAXON_MAP = {
    "Mouse": "http://purl.obolibrary.org/obo/NCBITaxon_10090",
    "ZebraFish": "http://purl.obolibrary.org/obo/NCBITaxon_7955",
    "Worm": "http://purl.obolibrary.org/obo/NCBITaxon_6239",
    "Fly": "http://purl.obolibrary.org/obo/NCBITaxon_7227",
    "Yeast": "http://purl.obolibrary.org/obo/NCBITaxon_559292",
    "Rat": "http://purl.obolibrary.org/obo/NCBITaxon_10116",
}


def main():

    hgnc = "ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_types/gene_with_protein_product.txt"

    protein_coding_genes = set()

    hgnc_resp = urllib.request.urlopen(hgnc).read()
    for line in hgnc_resp.decode("utf-8").splitlines():
        protein_coding_genes.add(re.split(r"\t", line)[0])

    print("Number of hgnc protein coding genes: {0}".format(len(protein_coding_genes)))

    human_causal = get_causal_gene_phenotype_assocs(protein_coding_genes)
    print(
        "Number of human genes with a causal g2p association: {0}".format(
            len(human_causal)
        )
    )

    # human_genes = get_human_genes()
    # print("Number of human gene cliques: {0}".format(len(human_genes)))

    human_genes_pheno = get_gene_phenotype_list("NCBITaxon:9606")

    for taxon, taxon_iri in TAXON_MAP.items():
        taxon_curie = map_iri_to_curie(taxon_iri)
        gene_list = get_gene_phenotype_list(taxon_curie)
        print("{0}: {1} gene counts".format(taxon, len(gene_list)))

    print(
        "Total human gene-phenotype/disease counts: {0}".format(len(human_genes_pheno))
    )

    model_only = set()
    all_models = dict()
    model_human_set = set()
    multi_model_set = set()

    for taxon, taxon_iri in TAXON_MAP.items():
        gene_counts = get_orthology_stats(taxon_iri)
        print(
            "Human-{0} orthology count, human: {1}, {2}: {3}".format(
                taxon, gene_counts[0]["human"], taxon, gene_counts[0]["ortholog"]
            )
        )

    for taxon, taxon_iri in TAXON_MAP.items():
        taxon_curie = map_iri_to_curie(taxon_iri)
        results = get_model_gene_stats(
            taxon_curie,
            protein_coding_genes,
            human_causal,
            model_human_set,
            model_only,
            multi_model_set,
        )

        model_only = results["model_only"]
        all_models[taxon] = results["model_set"]
        model_human_set = results["model_human_set"]
        print("{0}: {1} ortholog counts".format(taxon, len(results["model_set"])))
        # print("{0}: {1} unmatched count".format(taxon, len(results['unmatched_set'])))

    print("Models only: {0}".format(len(model_only)))
    print("Models AND human: {0}".format(len(model_human_set)))
    print("Human only: {0}".format(len(human_causal) - len(model_human_set)))
    print(
        "Models OR human: {0}".format(
            len(model_human_set)
            + len(model_only)
            + (len(human_causal) - len(model_human_set))
        )
    )
    print("Human total: {}".format(len(human_causal)))
    print("Models total: {0}".format(len(model_human_set) + len(model_only)))

    print("##########################")

    print(
        "Orthologs with >2 species and no human data: {0}".format(len(multi_model_set))
    )
    for taxon, taxon_iri in TAXON_MAP.items():
        one_species_count = (
            len(all_models[taxon])
            - len(all_models[taxon].intersection(multi_model_set))
            - len(all_models[taxon].intersection(model_human_set))
        )

        print("{0} only: {1}".format(taxon, one_species_count))


def get_model_gene_stats(
    taxon_curie,
    human_genes,
    human_genes_pheno,
    model_human_set,
    model_only,
    multi_model_set,
):

    results = {
        "model_set": set(),
        "model_only": model_only,
        "model_human_set": model_human_set,
        "multi_model_set": multi_model_set,
        "unmatched_set": set(),
    }

    filters = [
        "association_type:gene_phenotype",
        'subject_taxon: "{0}"'.format(taxon_curie),
    ]
    params = {
        "wt": "json",
        "rows": 1000,
        "start": 0,
        "q": "*:*",
        "fq": filters,
        "fl": "subject, subject_ortholog_closure",
    }
    resultCount = params["rows"]
    while params["start"] < resultCount:
        solr_request = requests.get(SOLR_URL, params=params)
        response = solr_request.json()
        resultCount = response["response"]["numFound"]

        for doc in response["response"]["docs"]:
            foundOrtholog = False
            if "subject_ortholog_closure" in doc:
                for ortholog in doc["subject_ortholog_closure"]:
                    if ortholog in human_genes_pheno:
                        results["model_human_set"].add(ortholog)
                        results["model_set"].add(ortholog)
                        foundOrtholog = True
                    elif ortholog in human_genes:
                        results["model_set"].add(ortholog)
                        if ortholog in results["model_only"]:
                            results["multi_model_set"].add(ortholog)
                        else:
                            results["model_only"].add(ortholog)
                        foundOrtholog = True
            if not foundOrtholog:
                results["unmatched_set"].add(doc["subject"])

        params["start"] += params["rows"]

    return results


def get_orthology_stats(taxon_iri):

    stats = dict()

    query = (
        "MATCH (:Node{{iri:'NCBITaxon:9606'}})<-[:RO:0002162]-(gene:gene)"
        "-[rel:RO:HOM0000017|RO:HOM0000020]-(ortholog:gene)-[:RO:0002162]->"
        "({{iri:'{0}'}}) "
        "RETURN COUNT(DISTINCT(ortholog)) as ortholog, COUNT(DISTINCT(gene)) as human".format(
            taxon_iri
        )
    )

    scigraph_service = SCIGRAPH_URL + "/cypher/execute.json"
    params = {"cypherQuery": query, "limit": 10}
    request = requests.get(scigraph_service, params=params)
    results = request.json()

    stats["human_gene_count"] = results[0]["human"]
    stats["model_gene_count"] = results[0]["ortholog"]
    return results


def get_gene_phenotype_list(taxon_curie):
    """
    Get a list of genes with phenotype or disease info indexed in solr
    :param taxon_curie:
    :return:
    """
    filters = [
        "association_type:({} OR {})".format("gene_phenotype", "gene_disease"),
        'subject_taxon: "{0}"'.format(taxon_curie),
    ]
    params = {
        "wt": "json",
        "rows": 0,
        "start": 0,
        "q": "*:*",
        "fq": filters,
        "facet": "true",
        "facet.mincount": 1,
        "facet.sort": "count",
        "json.nl": "arrarr",
        "facet.limit": -1,
        "facet.field": "subject",
    }
    solr_request = requests.get(SOLR_URL, params=params)
    response = solr_request.json()
    genes = {val[0] for val in response["facet_counts"]["facet_fields"]["subject"]}

    return genes


def get_human_genes():

    scigraph_service = SCIGRAPH_URL + "/cypher/execute.json"
    query = (
        "MATCH (gene:gene)-[tax:RO:0002162]->(taxon:Node{iri:'NCBITaxon:9606'}) "
        "RETURN gene.iri"
    )
    params = {"cypherQuery": query, "limit": 100000}

    request = requests.get(scigraph_service, params=params)
    try:
        results = request.json()
        genes = [key["gene.iri"] for key in results]
        gene_set = {
            map_iri_to_curie(val)
            for val in genes
            if not val.startswith("http://flybase.org")
        }
    except JSONDecodeError as e:
        print(request.text)
        gene_set = set()

    return gene_set


def get_causal_gene_phenotype_assocs(protein_coding_genes):
    print("Fetching causal human gene phenotype and disease associations")
    result_set = set()
    filters = [
        "association_type:({} OR {})".format("gene_phenotype", "gene_disease"),
        'subject_taxon: "{0}"'.format("NCBITaxon:9606"),
    ]
    params = {
        "wt": "json",
        "rows": 1000,
        "start": 0,
        "q": "*:*",
        "fq": filters,
        "fl": "subject, relation, is_defined_by, qualifier, association_type",
    }

    resultCount = params["rows"]
    while params["start"] < resultCount:
        solr_request = requests.get(SOLR_URL, params=params)
        response = solr_request.json()
        resultCount = response["response"]["numFound"]

        for doc in response["response"]["docs"]:

            if (
                doc["association_type"] == "gene_phenotype"
                and doc["qualifier"] != "direct"
            ):
                continue

            if doc["subject"] not in protein_coding_genes:
                continue

            result_set.add(doc["subject"])

        params["start"] += params["rows"]

    return result_set


def map_iri_to_curie(iri):
    curie = iri
    for prefix in CURIE_MAP:
        if iri.startswith(prefix):
            curie = re.sub(r"{0}".format(prefix), "{0}:".format(CURIE_MAP[prefix]), iri)
            break
    return curie


if __name__ == "__main__":
    main()
