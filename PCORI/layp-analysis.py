from monarch import monarch
import argparse
import logging
import requests
from json.decoder import JSONDecodeError

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

requests_log = logging.getLogger("requests.packages.urllib3")
requests_log.setLevel(logging.ERROR)

# Custom session for services that require more than 3 retries
session = requests.Session()
adapter = requests.adapters.HTTPAdapter(max_retries=10)
session.mount('https://', adapter)

parser = argparse.ArgumentParser(description='Generates annotation'
                                 'sufficiency scores for set of phenotype ids')
parser.add_argument('--phenotypes', '-p', type=str, required=True)
parser.add_argument('--diseases', '-d', type=str, required=True)
parser.add_argument('--output', '-o', type=str, required=False,
                    help='Location of output file', default="./output.tsv")

args = parser.parse_args()

OWLSIM_URL = 'https://monarchinitiative.org/owlsim/'

pheno_fh = open(args.phenotypes, 'r')
disease_fh = open(args.diseases, 'r')
output = open(args.output, 'w')

pheno_list = pheno_fh.read().splitlines()
disease_list = disease_fh.read().splitlines()

for disease in disease_list:
    # Get clique leader
    clique_leader = monarch.get_clique_leader(disease)
    mondo = clique_leader['id']
    mondo_label = clique_leader['label']

    # Get phenotypes
    pheno_profile = monarch.get_direct_phenotypes(mondo)

    lay_profile = set(pheno_profile).intersection(pheno_list)

    # Get annot sufficiency of whole profile
    scores = monarch.get_annotation_sufficiency_score(pheno_profile)
    disease_score = scores['scaled_score']

    if (len(lay_profile) == 0):
        output.write("{}\t{}\t{:0.3f}\t{}\t{}\t{}\t{}\n".format(
            mondo,
            mondo_label,
            disease_score,
            "",
            "",
            "",
            ""
        ))
        continue

    # Get annot sufficiency of lay profile
    scores = monarch.get_annotation_sufficiency_score(lay_profile)
    lay_score = scores['scaled_score']

    # Get sim score
    compare_url = OWLSIM_URL + "compareAttributeSets"
    params = {
        'a': pheno_profile,
        'b': lay_profile
    }
    sim_req = session.post(compare_url, data=params)
    try:
        owlsim_results = sim_req.json()
        sim_score = owlsim_results['results'][0]['combinedScore']
    except JSONDecodeError:
        sim_score = 0
    except IndexError:
        sim_score = 0
    annot_suff_ratio = lay_score / disease_score
    owlsim_ratio = ((sim_score/100) + annot_suff_ratio)/2

    output.write("{}\t{}\t{:0.3f}\t{:0.3f}\t{:0.3f}\t{}\t{:0.3f}\n".format(
        mondo,
        mondo_label,
        disease_score,
        lay_score,
        annot_suff_ratio,
        sim_score,
        owlsim_ratio
    ))

