#!/usr/bin/env python3

import xml.etree.ElementTree as etree
import requests
import argparse
import sys

parser = argparse.ArgumentParser(description='description')
parser.add_argument('--input', '-i', type=str, required=True,
                    help='Location of input file')
args = parser.parse_args()

BIOMART = 'http://www.ensembl.org/biomart/martservice'

query_attributes = {
    "virtualSchemaName": "default",
    "formatter": "TSV",
    "header": "0",
    "uniqueRows": "1",
    "count": "0",
    "datasetConfigVersion": "0.6"
}

object_attributes = {
    "name": 'hsapiens_gene_ensembl',
    "interface": "default"
}

# prepend the query
prepend = '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'

file = open(args.input, 'r')
for transcript in file:
    transcript = transcript.rstrip()

    query = etree.Element("Query", query_attributes)
    dataset = etree.SubElement(query, "Dataset", object_attributes)
    attributes = etree.SubElement(dataset, "Attribute", {"name": "entrezgene"})
    attributes = etree.SubElement(dataset, "Attribute", {"name": "external_gene_name"})

    filters = etree.SubElement(dataset, "Filter", {"name": 'ensembl_transcript_id',
                                                   "value": transcript})

    query_string = etree.tostring(query, encoding="unicode")
    query_string = prepend + query_string

    params = {'query':  query_string}

    response = requests.get(BIOMART, params)
    res_text = response.text.strip('\n')
    rows = res_text.split('\n')
    if len(rows) > 1:
        print("ambiguous mapping: {0}".format(transcript), file=sys.stderr)
    for row in rows:
        if len(row) == 0:
            print(transcript)
            continue
        records = row.split('\t')
        entrez_id = "NCBIGene:{0}".format(records[0])
        print("{0}\t{1}\t{2}".format(transcript, entrez_id, records[1]))


