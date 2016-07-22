#!/usr/bin/env python3

import argparse
import csv
import logging
from monarch import monarch

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


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
    disease_dictionary = monarch.get_zfin_ids(disease_dictionary)
    logger.info("Getting disease leaf/group info")
    disease_dictionary = monarch.get_disease_info(disease_dictionary)
    logger.info("Getting disease-phenotype, gene-phenotype, and disease-gene counts from solr")
    disease_dictionary = monarch.get_solr_counts(disease_dictionary)
    logger.info("Getting owlsim score/rank")
    disease_dictionary = monarch.get_owlsim_scores(disease_dictionary)

    output_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\n".format(
                      "disease", "disease_label", "model_gene", "gene_label", "model_gene_taxon",
                      "owlsim_score", "owlsim_rank", "models", "disease_pheno_count",
                      "disease_gene_count", "modgene_pheno_count", "isLeafNode"
                      ))

    for disease_id, disease in disease_dictionary.items():
        output_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\n".format(
                          disease["disease"], disease["disease_label"],
                          disease["model_gene"], disease["model_gene_label"],
                          monarch.TAXON_MAP[disease["model_gene_taxon"]],
                          disease["owlsim_score"], disease["owlsim_rank"],
                          disease["models"], disease["disease_pheno_count"],
                          disease["disease_gene_count"], disease["modgene_pheno_count"],
                          disease["isLeafNode"]
        ))


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
            (model, model_label, model_gene, model_gene_label, disease, disease_label)\
                = [monarch.map_iri_to_curie(val) for val in row]

            key = "{0}-{1}".format(disease, model_gene)
            gene_taxon = monarch.get_taxon(model_gene)


            if gene_taxon == "NCBITaxon:10116" and model.startswith("MGI"):
                # Catch and skip over rat transgenes
                continue

            if key not in disease_dictionary:
                disease_dictionary[key] = {
                    "disease": disease,
                    "disease_label": disease_label,
                    "model_gene": model_gene,
                    "model_gene_label": model_gene_label,
                    "models": 1,
                    "model_gene_taxon": gene_taxon
                }
            else:
                disease_dictionary[key]['models'] += 1

    return disease_dictionary


if __name__=="__main__":
    main()
