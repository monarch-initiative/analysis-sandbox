import argparse
import logging
import requests
import monarch
import mygene


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

"""
Filter genes as potential candidates for co-crispr
experiments in mice
"""

MONARCH_SEARCH = 'https://solr.monarchinitiative.org/solr/search/select'
SIMSEARCH_URL = 'https://monarchinitiative.org/simsearch/phenotype'
GOLR_URL = 'https://solr-dev.monarchinitiative.org/solr/golr/select'
SCIGRAPH_URL = 'https://scigraph-data-dev.monarchinitiative.org/scigraph/cypher/execute.json'

CURIE_MAP = {
    "OMIM": "http://purl.obolibrary.org/obo/OMIM_",
    "DOID": "http://purl.obolibrary.org/obo/DOID_",
    "MGI": "http://www.informatics.jax.org/accession/MGI:",
    "NCBIGene": "http://www.ncbi.nlm.nih.gov/gene/",
    "ZFIN": "http://zfin.org/"
}

"""
# Convert kinase, ion channel, gpcr gene symbols to ids
# Get mouse orthologs
# Find all phenotypically similar genes (>85)
# co expression (BGEE)
# Check if pairwise genes directly interact (biogrid), add qualifier
# Check if pairwise genes are in same pathway (reactome, kegg, ctd), add pipe delim pathways
"""


def main():

    parser = argparse.ArgumentParser(
        description='Converts UDN patient dump file into table'
    )
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='Location of input file with gene symbols')
    parser.add_argument('--output', '-o', type=str, required=True,
                        help='Location of output file')
    parser.add_argument('--taxon', '-t', type=str,
                        help='NCBI taxon ID', default="NCBITaxon:10090")

    args = parser.parse_args()
    mg = mygene.MyGeneInfo()

    input_file = open(args.input, 'r')
    output_file = open(args.output, 'w')
    output_file.write("gene label\tortholog id\tortholog label\t"
                      "similar gene label\tsimilarity score\tdirect interactors\t"
                      "pathway_labels\tco-expressed in tissue\t"
                      "similar gene id\tpathway id\ttissue ids\tphenotype_count\n")

    symbol_list = []
    gene_ids = []
    for symbol in input_file:
        symbol = symbol.rstrip('\n')
        symbol_list.append(symbol)

    logger.info("Querying mygene")
    mygene_results = mg.querymany(symbol_list, scopes='symbol', fields='entrezgene', species='human', returnall=True)

    mygene_output = mygene_results['out']

    for symbol in symbol_list:
        isFound = False
        for result in mygene_output:
            if symbol == result['query'] and 'entrezgene' in result:
                gene_ids.append("NCBIGene:{0}".format(result['entrezgene']))
                isFound = True
                break
        if not isFound:
            gene_ids.append(None)

    ortholog_ids = []
    ortholog_labels = []
    for gene in gene_ids:
        ortholog_id, ortholog_label = get_orthologs(gene, args.taxon)
        ortholog_ids.append(ortholog_id)
        ortholog_labels.append(ortholog_label)

    for gene_group in zip(symbol_list, gene_ids, ortholog_ids, ortholog_labels):
        for ortholog, ortholog_label\
                in zip(gene_group[2].split("|"), gene_group[3].split("|")):
            phenotype_profile = monarch.get_phenotype_profile(ortholog)
            matches = get_top_similar_genes(phenotype_profile, args.taxon)
            for match in matches:
                pathway_ids, pathway_labels \
                    = get_pathways(ortholog, match['id'])
                print(pathway_ids, pathway_labels)
                biogrid_interactors = get_direct_interactors(ortholog,
                                                             match['id'])
                print(biogrid_interactors)
                uberon_ids, uberon_labels \
                    = get_tissue_coexpression(ortholog, match['id'])
                print(uberon_ids, uberon_labels)
                output_file.write("{}\t{}\t{}\t{}\t{}\t{}"
                                  "\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                    gene_group[0], ortholog, ortholog_label,
                                    match['label'], match['score'],
                                    biogrid_interactors, pathway_labels,
                                    uberon_labels, match['id'],
                                    pathway_ids, uberon_ids, len(phenotype_profile)))


def get_pathways(gene1, gene2):
    prefix, id = gene1.split(':')
    gene1_iri = CURIE_MAP[prefix] + id
    prefix, id = gene2.split(':')
    gene2_iri = CURIE_MAP[prefix] + id

    query = "MATCH (gene1:Node{{iri:'{}'}})-[:RO:0002205*0..1]->()-[:RO:0002331]->(pathway:pathway)<-" \
            "[:RO:0002331]-()<-[:RO:0002205*0..1]-(gene2:Node{{iri:'{}'}}) " \
            "WHERE gene1 <> gene2 " \
            "RETURN pathway".format(gene1_iri, gene2_iri)

    params = {
        'cypherQuery': query,
        'limit': 1000
    }

    request = requests.get(SCIGRAPH_URL, params=params)
    results = request.json()
    pathway_ids = []
    pathway_labels = []
    if len(results) > 0:
        pathway_ids = [result['pathway']['iri'] for result in results]
        pathway_labels = [result['pathway']['label'] for result in results]

    return "|".join(pathway_ids), "|".join(pathway_labels)


def get_orthologs(gene, taxon):
    params = {
        'wt': 'json',
        'rows': 100,
        'start': 0,
        'q': '*:*',
        'fl': 'subject, subject_label,'
              'object, object_label',
        'fq': ['subject_closure: "{0}"'.format(gene),
               'relation_closure: "RO:HOM0000017"',
               'object_taxon: "{0}"'.format(taxon)
        ]
    }
    solr_request = requests.get(GOLR_URL, params=params)
    response = solr_request.json()
    resultCount = response['response']['numFound']
    if resultCount > 1:
        logger.info("More than one ortholog found for {}".format(gene))
    gene_ids = []
    gene_labels = []
    for doc in response['response']['docs']:
        gene_ids.append(doc['object'])
        gene_labels.append(doc['object_label'])

    return "|".join(gene_ids), "|".join(gene_labels)


def get_top_similar_genes(pheno_profile, taxon=None, cutoff=75):
    phenotypes = "+".join(pheno_profile)
    params = "input_items={}".format(phenotypes)
    if taxon is not None:
        params += "&target_species={}".format(taxon.replace("NCBITaxon:", ""))

    owlsim_request = requests.post(SIMSEARCH_URL, params=params)
    owlsim_results = owlsim_request.json()
    matches = []
    if 'b' in owlsim_results:
        for index, result in enumerate(owlsim_results['b']):
            if index == 0:
                continue
            if result['score']['score'] >= cutoff:
                result_slim = {
                    'id': result['id'],
                    'label': result['label'],
                    'score': result['score']['score']
                }
                matches.append(result_slim)
    return matches


def get_tissue_coexpression(gene1, gene2):
    prefix, id = gene1.split(':')
    gene1_iri = CURIE_MAP[prefix] + id
    prefix, id = gene2.split(':')
    gene2_iri = CURIE_MAP[prefix] + id

    query = "MATCH (gene1:Node{{iri:'{}'}})-[:RO:0002206]->(tissue:Node)<-[:RO:0002206]-" \
            "(gene2:Node{{iri:'{}'}}) " \
            "WHERE gene1 <> gene2 " \
            "RETURN tissue".format(gene1_iri, gene2_iri)

    params = {
        'cypherQuery': query,
        'limit': 1000
    }
    request = requests.get(SCIGRAPH_URL, params=params)
    results = request.json()
    uberon_ids = []
    uberon_labels = []
    if len(results) > 0:
        uberon_ids = [result['tissue']['iri']
                        .replace('http://purl.obolibrary.org/obo/', '')
                      for result in results]
        uberon_labels = [result['tissue']['label'] for result in results if 'label' in result['tissue']]

    return "|".join(uberon_ids), "|".join(uberon_labels)


def get_ncbi_id_from_symbol(gene_symbol, taxon="NCBITaxon:9606"):
    """
    Get ncbi gene id from symbol using monarch services
    :param gene_symbol: str, gene symbol
    :param taxon: str, ncbi taxon id formatted
                  as a curie with prefix NCBITaxon, ie NCBITaxon:1234
    :return: str, NCBI gene id as curie, ie NCBIGene:1234
    """
    params = _get_solr_weight_settings()
    params["q"] = "{0} \"{0}\"".format(gene_symbol)
    params["fq"] = ["taxon:\"{0}\"".format(taxon), "category:\"gene\""]
    gene_id = None
    try:
        monarch_request = requests.get(MONARCH_SEARCH, params=params)
        response = monarch_request.json()
        count = response['response']['numFound']
        if count > 0:
            gene_id = response['response']['docs'][0]['id']
    except requests.ConnectionError:
        logger.error("error fetching {} params: {}".format(MONARCH_SEARCH, params))

    return gene_id


def get_direct_interactors(gene1, gene2):

    areInteractors = False
    params = {
        'wt': 'json',
        'rows': 0,
        'start': 0,
        'q': '*:*',
        'fl': 'object',
        'fq': ['subject_closure:"{0}"'.format(gene1),
               'object_closure:"{0}"'.format(gene2),
               'relation_closure:"RO:0002434"']
    }

    solr_request = requests.get(GOLR_URL, params=params)
    response = solr_request.json()
    resultCount = response['response']['numFound']
    if resultCount > 0:
        areInteractors = True

    params['subject_closure'] = '"{0}"'.format(gene2)
    params['object_closure'] = '"{0}"'.format(gene1)

    solr_request = requests.get(GOLR_URL, params=params)
    response = solr_request.json()
    resultCount = response['response']['numFound']
    if resultCount > 0:
        areInteractors = True

    return 1 if areInteractors is True else 0


def _get_solr_weight_settings():
    return {
        "qt": "standard",
        "json.nl": "arrarr",
        "fl": "*,score",
        "start": "0",
        "rows": "5",
        "defType": "edismax",
        "personality": "monarch_search",
        "qf": [
            "label_searchable^1",
            "definition_searchable^1",
            "synonym_searchable^1",
            "label_std^2",
            "synonym_std^1"
        ],
        "wt": "json"
    }


if __name__ == "__main__":
    main()
