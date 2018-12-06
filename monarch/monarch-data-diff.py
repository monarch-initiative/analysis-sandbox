#!/usr/bin/python3
""" Monarch Data Diff

Generate a diff of solr facet counts
and scigraph queries as markdown and html

TODO support lists of facets
"""
import requests
from requests import Request, Session
import argparse
from pathlib import Path
import logging
import json
from typing import List, Tuple, Iterable, Dict, Union, Any, Set
import copy
import itertools
import re
import markdown
from json.decoder import JSONDecodeError


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

JSONType = Union[str, int, float, bool, None, Dict[str, Any], List[Any]]


def main():

    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--config', '-c', required=True, help='JSON configuration file')
    parser.add_argument('--out', '-o', required=False, help='output directory', default="./")
    parser.add_argument('--threshold', '-t', required=False, help='diff threshold', default=10, type=int)
    parser.add_argument('--quick', '-q', required=False, help='diff threshold',
                        default=False, action='store_true')
    args = parser.parse_args()
    dir_path = Path(args.out)
    threshold = float(args.threshold/100)

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
        'facet.sort': 'count',
        'indent': 'on'
    }

    golr_default_params = {
        'q': '*:*',
        'wt': 'json',
        'fl': 'subject,object'
    }

    conf_fh = open(args.config, 'r')
    config = json.load(conf_fh)

    solr_dev = config['solr-dev']
    solr_prod = config['solr-prod']
    scigraph_data_dev = config['scigraph-data-dev']
    scigraph_data_prod = config['scigraph-data-prod']
    scigraph_ontology_dev = config['scigraph-ontology-dev']
    scigraph_ontology_prod = config['scigraph-ontology-prod']

    md_file.write("{}\n".format(add_md_header("Solr Queries", 3)))

    # Process solr queries
    for q_name, query in config['solr_facet_queries'].items():
        golr_facet_params['fq'] = query['filters']
        golr_facet_params['facet.field'] = query['facet.field']
        prod_results = get_facets(solr_prod, golr_facet_params)
        dev_results = get_facets(solr_dev, golr_facet_params)
        diff = diff_facets(dev_results, prod_results)
        formatted_diff = convert_diff_to_md(diff)

        # Create a filtered object based on the threshold passed in (or 10%)
        # If there is a x% decrease in data pass this on to generate a data
        # diff
        if args.quick:
            filtered_diff = {}
        else:
            filtered_diff = {k:v for k,v in diff.items()
                             if (v[0] == 0 and v[1] < 0) or
                             (v[1] < 0 and 1 - v[0] / (v[0] + abs(v[1])) >= threshold)}

        for dropped_data in filtered_diff:
            if dropped_data == 'Total' or \
                   dropped_data.startswith(":.well-known"):
                continue

            diff_dir = dir_path / q_name.replace(' ', '_').lower()
            diff_dir.mkdir(parents=True, exist_ok=True)

            diff_path = diff_dir / "{}.tsv".format(dropped_data.replace(' ', '_').lower())
            diff_file = diff_path.open("w")

            params = copy.deepcopy(golr_default_params)
            params['fq'] = copy.deepcopy(query['filters'])
            if dropped_data == 'Other':
                params['fq'].append('-{}:[* TO *]'.format(query['facet.field'], dropped_data))
            else:
                params['fq'].append('{}:"{}"'.format(query['facet.field'], dropped_data))
            logger.info("Processing data diff for params {}".format(params['fq']))

            query_pairs = get_solr_so_pairs(solr_dev, params)
            reference_pairs = get_solr_so_pairs(solr_prod, params)

            solr_diff = diff_solr_so_data(query_pairs, reference_pairs)
            for sub, objects in solr_diff.items():
                for obj in objects:
                    diff_file.write("{}\t{}\n".format(sub, obj))
            diff_file.close()

        md_file.write("{}\n".format(add_md_header(q_name, 4)))
        sesh = Session()
        prod_req = sesh.prepare_request(Request('GET', solr_prod, params=golr_facet_params))
        dev_req = sesh.prepare_request(Request('GET', solr_dev, params=golr_facet_params))

        md_file.write(add_href(prod_req.url, "Production Query"))
        md_file.write('\n\n')
        md_file.write(add_href(dev_req.url, "Dev Query"))
        md_file.write('\n\n')

        diff_list = [(k, v) for k, v in formatted_diff.items()]
        diff_list.sort(key=lambda tup: int(re.search(r'\d+', tup[1]).group(0)), reverse=True)
        md_file.write(add_md_table(diff_list, query['headers']))
        md_file.write("\n\n")

    md_file.write("{}\n".format(add_md_header("SciGraph Queries", 3)))

    for q_name, query in config['scigraph_data_queries'].items():
        md_file.write(get_scigraph_diff(
            scigraph_data_prod, scigraph_data_dev, query, q_name))

    for q_name, query in config['scigraph_ontology_queries'].items():
        md_file.write(get_scigraph_diff(
            scigraph_ontology_prod, scigraph_ontology_dev, query, q_name))

    md_file.close()
    md_file = md_path.open("r")
    html = markdown.markdown(md_file.read(), output_format='html5', extensions=['markdown.extensions.tables'])
    html_file.write(html)
    html_file.close()
    md_file.close()


def get_scigraph_diff(scigraph_prod: str, scigraph_dev: str,
                      conf: dict, query_name: str) -> str:
    output_md = str()
    prod_results = get_scigraph_results(scigraph_prod, conf['query'])
    dev_results = get_scigraph_results(scigraph_dev, conf['query'])
    diff = diff_facets(dev_results, prod_results)
    formatted_diff = convert_diff_to_md(diff)

    output_md += "{}\n".format(add_md_header(query_name, 4))

    params = {
        'cypherQuery': conf['query']
    }

    sesh = Session()
    prod_req = sesh.prepare_request(Request('GET', scigraph_prod, params=params))
    dev_req = sesh.prepare_request(Request('GET', scigraph_dev, params=params))

    output_md += add_href(prod_req.url, "Production Query")
    output_md += '\n\n'
    output_md += add_href(dev_req.url, "Dev Query")
    output_md += '\n\n'

    diff_list = [(k, v) for k, v in formatted_diff.items()]
    diff_list.sort(key=lambda tup: int(re.search(r'\d+', tup[1]).group(0)), reverse=True)
    output_md += add_md_table(diff_list, conf['headers'])
    output_md += "\n\n"

    return output_md


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


def add_href(link: str, display_name: str) -> str:
    return "[{}]({})".format(display_name, link)


def get_facets(solr_server: str,
               params: Dict[str, Union[str, List[str]]]) -> Dict[str, int]:
    solr_request = requests.get(solr_server, params=params)
    response = solr_request.json()
    facet = params['facet.field']
    result_count = response['response']['numFound']
    facet_result = response['facet_counts']['facet_fields'][facet]
    facet_obj = dict(itertools.zip_longest(*[iter(facet_result)] * 2, fillvalue=""))
    res_sum = sum([v for k, v in facet_obj.items()])
    other_count = result_count - res_sum
    facet_obj['Total'] = result_count
    facet_obj['Other'] = other_count
    return facet_obj


def get_scigraph_results(scigraph_server: str, query: str) -> JSONType:
    params = {
        'cypherQuery': query
    }
    scigraph_request = requests.get(scigraph_server, params=params, timeout=120)
    try:
        response = scigraph_request.json()
    except JSONDecodeError as json_exc:
         raise JSONDecodeError(
                "Cannot parse scigraph response for {}: {}".format(scigraph_request.url, json_exc.msg),
                json_exc.doc,
                json_exc.pos
         )

    results = response[0]
    return results


def diff_facets(query: Dict, reference: Dict) -> Dict[str, List[int]]:
    diff_obj = copy.deepcopy(query)
    for key in reference:
        if key not in query:
            query[key] = 0
    for key in query:
        if key not in reference:
            reference[key] = 0
        diff = query[key] - reference[key]
        diff_obj[key] = [query[key], diff]

    return diff_obj


def convert_diff_to_md(diff: Dict[str, List[int]]) -> Dict[str, str]:
    diff_obj = copy.deepcopy(diff)
    for k, v in diff.items():
        count, diff = v
        if diff == 0:
            diff_obj[k] = "{} ({})".format(count, diff)
        elif diff > 0:
            diff_obj[k] = "{} (+{})".format(count, add_bold(diff))
        else:
            diff_obj[k] = "{} (-{})".format(count, add_italics(abs(diff)))
    return diff_obj


def diff_solr_so_data(query: Dict[str, Set[str]],
                      reference: Dict[str, Set[str]]) -> Dict[str, Set[str]]:
    results = {}
    for sub in reference:
        results[sub] = set()
        try:
            results[sub] = reference[sub] - query[sub]
        except KeyError:
            results[sub] = reference[sub]

    return results


def get_solr_so_pairs(
        solr_server: str,
        params: Dict[str, Union[str, List[str]]]) -> Dict[str, Set[str]]:
    results = {}
    params['start'] = 0
    params['rows'] = 1000

    result_count = params['rows']

    while params['start'] < result_count:
        solr_request = requests.get(solr_server, params=params)
        response = solr_request.json()
        result_count = response['response']['numFound']
        if result_count == 0:
            break
        for doc in response['response']['docs']:
            try:
                results[doc['subject']].add(doc['object'])
            except KeyError:
                results[doc['subject']] = {doc['object']}

        params['start'] += params['rows']
    return results


if __name__ == "__main__":
    main()
