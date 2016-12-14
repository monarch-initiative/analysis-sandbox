#!/usr/bin/env python3
import argparse
import requests


def main():

    SCIGRAPH_BASE = "http://localhost:9000/scigraph/vocabulary/search/"

    parser = argparse.ArgumentParser(description='description',
                                     formatter_class=
                                     argparse.RawTextHelpFormatter)
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='Location of input file')
    parser.add_argument('--output', '-o', type=str,
                        help='Location of output file')

    args = parser.parse_args()

    # Open file
    input_file = open(args.input, 'r')
    output_file = open(args.output, 'w')
    output_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n"
                      .format("Trait", "Match", "Match Label",
                              "Second Match", "Second Label",
                              "Third Match", "Third Label"))

    scigraph_params = {
        'limit': 20,
        'searchSynonyms': 'true',
        'searchAbbreviations': 'true',
        'searchAcronyms': 'false',
        'prefix': 'CHEBI'
    }

    for trait in input_file:
        trait = trait.rstrip()
        scigraph_url = SCIGRAPH_BASE + trait
        scigraph_request = requests.get(scigraph_url, params=scigraph_params)
        try:
            response = scigraph_request.json()
            match_curie = response[0]['curie']
            match_label = response[0]['labels'][0]
            output_file.write("{0}\t{1}\t{2}".format(trait,
                                                       match_curie,
                                                       match_label))
            if len(response) > 1:
                match_curie = response[1]['curie']
                match_label = response[1]['labels'][0]
                output_file.write("\t{0}\t{1}".format(match_curie,
                                                      match_label))

            if len(response) > 2:
                match_curie = response[2]['curie']
                if len(response[2]['labels']) > 0:
                    match_label = response[2]['labels'][0]
                else:
                    match_label = ""
                output_file.write("\t{0}\t{1}".format(match_curie,
                                                      match_label))
            output_file.write("\n")

        except ValueError:
            output_file.write(trait+"\n")

    input_file.close()
    output_file.close()


if __name__ == "__main__":
    main()
