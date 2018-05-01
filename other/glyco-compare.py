import numpy as np
from monarch import monarch
import requests

fucosidosis = 'MONDO:0009254'

output = open('./fuco-analysis.tsv', 'w')

output.write("{}\t{}\t{}\t{}\t{}\n".format(
        "sample_score",
        "sample_rank",
        "glyco_score",
        "glyco_rank",
        "sample"
))

def get_sim_and_rank(pheno_list, disease):
    OWLSIM_URL = 'https://beta.monarchinitiative.org/owlsim/'
    search_url = OWLSIM_URL + "searchByAttributeSet"
    params = {
        'a': pheno_list,
        'target': 'MONDO'
    }
    sim_req = requests.post(search_url, data=params)
    sim_results = sim_req.json()
    sim_list = sim_results['results']
    rank = 1
    last_score = 0
    sample_rank = 'NA'
    sample_score = 'NA'
    for res in sim_list:
        if res["j"]["id"] == fucosidosis:
            sample_rank = rank
            sample_score = res["combinedScore"]
            break
        elif res["combinedScore"] != last_score:
            rank += 1
            last_score = res["combinedScore"]

    return sample_score, sample_rank


phenotypes = [
    'HP:0010864',
    'HP:0001263',
    'HP:0000975',
    'HP:0008430',
    'HP:0100578',
    'HP:0000280',
    'HP:0002808',
    'HP:0000248',
    'HP:0000943',
    'HP:0000821',
    'HP:0011220',
    'HP:0002240',
    'HP:0011276',
    'HP:0000365',
    'HP:0005595',
    'HP:0001508'
]

glyco_phenotypes = [
    'HP:0010471',
    'HP:0003541'
]

pheno_profile = set(monarch.get_direct_phenotypes(fucosidosis))

seen = []
i = 0
while(i < 100):
    sample = np.random.choice(phenotypes, 10, False)
    if set(sample) in seen:
        continue
    seen.append(set(sample))

    sample_sim, sample_rank = get_sim_and_rank(sample, fucosidosis)
    w_glyco = np.append(sample, glyco_phenotypes)
    sim, rank = get_sim_and_rank(w_glyco, fucosidosis)

    output.write("{}\t{}\t{}\t{}\t{}\n".format(
        sample_sim,
        sample_rank,
        sim,
        rank,
        "|".join(sample)
    ))
    i += 1

