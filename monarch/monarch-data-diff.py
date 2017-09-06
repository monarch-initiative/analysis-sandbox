#!/usr/bin/python3
""" Monarch Data Diff

Generate a diff of solr facet counts
and scigraph queries as markdown and html
"""
import requests
import argparse
from pathlib import Path
import logging
import json
from typing import List, Tuple, Iterable
import copy
import itertools
import re
import markdown

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--config', '-c', required=True, help='JSON configuration file')
    parser.add_argument('--out', '-o', required=False, help='output directory', default="./")
    args = parser.parse_args()
    dir_path = Path(args.out)

    md_path = dir_path / 'monarch-diff.md'
    md_file = md_path.open("w")

    html_path = dir_path / 'monarch-diff.html'
    html_file = html_path.open("w")

    golr_facet_params = {
        'q': '*:*',
        'wt': 'json',
        'rows': '0',
        'facet': 'true',
        'facet.limit': '3000',
        'facet.method': 'enum',
        'facet.mincount': '1',
        'facet.sort': 'count'
    }

    conf_fh = open(args.config, 'r')
    config = json.load(conf_fh)

    solr_dev = config['solr-dev']
    solr_prod = config['solr-prod']
    scigraph_data_dev = config['scigraph-data-dev']
    scigraph_data_prod = config['scigraph-data-prod']
    scigraph_ontolgoy_dev = config['scigraph-ontology-dev']
    scigraph_ontology_prod = config['scigraph-ontology-prod']

    md_file.write("{}\n".format(add_md_header("Solr Queries", 3)))

    # Process solr queries
    for q_name, query in config['solr_facet_queries'].items():
        golr_facet_params['fq'] = query['filters']
        golr_facet_params['facet.field'] = query['facet.field']
        prod_results = get_facets(solr_prod, golr_facet_params)
        dev_results = get_facets(solr_dev, golr_facet_params)
        diff = diff_facets(dev_results, prod_results)

        md_file.write("{}\n".format(add_md_header(q_name, 4)))
        diff_list = [(k, v) for k, v in diff.items()]
        diff_list.sort(key=lambda tup: int(re.search(r'\d+', tup[1]).group(0)), reverse=True)
        md_file.write(add_md_table(diff_list, query['headers']))
        md_file.write("\n\n")

    md_file.write("{}\n".format(add_md_header("SciGraph Queries", 3)))

    for q_name, query in config['scigraph_data_queries'].items():
        prod_results = get_scigraph_results(scigraph_data_prod, query['query'])
        dev_results = get_scigraph_results(scigraph_data_dev, query['query'])
        diff = diff_facets(dev_results, prod_results)

        md_file.write("{}\n".format(add_md_header(q_name, 4)))
        diff_list = [(k, v) for k, v in diff.items()]
        diff_list.sort(key=lambda tup: int(re.search(r'\d+', tup[1]).group(0)), reverse=True)
        md_file.write(add_md_table(diff_list, query['headers']))
        md_file.write("\n\n")

    for q_name, query in config['scigraph_ontology_queries'].items():
        prod_results = get_scigraph_results(scigraph_data_prod, query['query'])
        dev_results = get_scigraph_results(scigraph_data_dev, query['query'])
        diff = diff_facets(dev_results, prod_results)

        md_file.write("{}\n".format(add_md_header(q_name, 4)))
        diff_list = [(k, v) for k, v in diff.items()]
        diff_list.sort(key=lambda tup: int(re.search(r'\d+', tup[1]).group(0)), reverse=True)
        md_file.write(add_md_table(diff_list, query['headers']))
        md_file.write("\n\n")

    md_file.close()
    md_file = md_path.open("r")
    html = markdown.markdown(md_file.read(), output_format='html5', extensions=['markdown.extensions.tables'])
    html_file.write(html)
    html_file.close()
    md_file.close()



def add_md_table(data: Iterable[Tuple], headers: List[str]=None) -> str:
    """
    Convert tuple to markdown table
    :param data: Tuple[Tuple]
    :param headers: List[str]
    :return: str
    """
    table = '| {} |\n'.format(' | '.join(str(header) for header in headers))
    table += '|-'*(len(headers)) + '|\n'
    for row in data:
        table += '| {} |\n'.format(' | '.join(str(cell) for cell in row))

    return table


def add_md_header(header: str, size: int=1) -> str:
    """
    :param header: str
    :param size: int
    :return: str
    """
    if size not in range(7)[1:7]:
        logger.error("Invalid header size: {}".format(size))
    return "{} {}".format(('#'*size), header)


def add_italics(text: str) -> str:
    return "_{}_".format(text)


def add_bold(text: str) -> str:
    return "__{}__".format(text)


def get_facets(solr_server: str, params:dict) -> dict:
    solr_request = requests.get(solr_server, params=params)
    response = solr_request.json()
    facet = params['facet.field']
    result_count = response['response']['numFound']
    facet_result = response['facet_counts']['facet_fields'][facet]
    facet_obj = dict(itertools.zip_longest(*[iter(facet_result)] * 2, fillvalue=""))
    res_sum = sum([v for k,v in facet_obj.items()])
    other_count = result_count - res_sum
    facet_obj['Total'] = result_count
    facet_obj['Other'] = other_count
    return facet_obj

def get_scigraph_results(scigraph_server: str, query:dict) -> dict:
    params = {
        'cypherQuery': query
    }
    scigraph_request = requests.get(scigraph_server, params=params)
    response = scigraph_request.json()
    results = response[0]
    return results


def diff_facets(query: dict, reference: dict) -> dict:
    diff_obj = copy.deepcopy(query)
    for key in reference:
        if key not in query:
            query[key] = 0
    for key in query:
        if key not in reference:
            reference[key] = 0
        diff = query[key] - reference[key]
        if diff == 0:
            diff_obj[key] = "{} ({})".format(query[key], diff)
        elif diff > 0:
            diff_obj[key] = "{} (+{})".format(query[key], add_bold(diff))
        else:
            diff_obj[key] = "{} (-{})".format(query[key], add_italics(abs(diff)))

    return diff_obj


if __name__ == "__main__":
    main()