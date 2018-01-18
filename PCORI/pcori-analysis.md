### PCORI Lay Person Analysis

#### Background

As of the December 2017 HPO release, 4555/12868 (35.4%) of HPO terms contain a lay person synonym.

From the paper:
In total, 36% of the HPO terms from the (2017-06-21) release have at least one layperson synonym (4547 of 12,623). 89% of the diseases (8666 of 9657) in the HPO database have at least one HPO annotation with a layperson synonym and 60% of all disease annotations (73,932 of 122,120) are referring to HPO terms with lay translations.

This analysis is designed to identify diseases that are diagnosable given the set of terms with lay person synonyms, or a subset provided in a survey platform, such genome connect or phenotypr.

Link to preliminary lay person paper:

https://docs.google.com/document/d/1hx7iA-1g5tf9-OybTeCd1eSJ1bB-v8KDg2jWgrKrp7g/edit


#### Files

__lay_person_analysis.tsv__

Fields:

A. MONDO Disease ID

B. Disease Label

C. Annotation Sufficiency of phenotypes annotated to disease, Range 0-1

D. Annotation Sufficiency of subset of annotated phenotypes that have lay person syns, Range 0-1

E. D/C (normalized annotation suff) 0-1

F. OwlSim Score (PhenoDigm) Range 0-100

G. AVG(E,(F/100)) Range 0-1

H. Number of phenotypes annotated to disease

I. Number of intersecting phenotypes


__genome_connect_analysis.tsv__

Fields:

same as lay person analysis

#### Summary of results

We analyzed 7717 rare diseases and their direct phenotype annotations. In this preliminary analysis, only exact phenotype overlaps were evaluated; in cases where no exact matches were found, fuzzy matches such as parent/sibling terms were not evaluated.

Diseases were binned into 6 groups:
1. No overlapping phenotypes
2. Full coverage of phenotypes
3. .75 - .99 normalized diagnosability score (G)
4. .50 - .74
5. .25 - .49
6. .01 - .24

__Lay Person__

1. No overlapping phenotypes: 753, 10.8%
2. Full coverage of phenotypes: 642, 8.3%
3. .75-.99: 4691, 60.8%
4. .50-.74: 1456, 18.9%
5. .25-.49: 139,  1.8%
6. .01-.24: 36, .4%


__Genome Connect__

The genome connect survey contains a subset of 216 HPO terms.

1. No overlapping phenotypes: 2995, 38.8%
2. Full coverage of phenotypes: 183, 2.4%
3. .75-.99: 179, 2.3%
4. .50-.74: 1851, 23.99%
5. .25-.49: 2311,  29.95%
6. .01-.24: 198, 2.6%

#### Disease Enrichment
During the analysis, the question came up if a disease group(s) are enriched in the set of non-annotated terms.  It was hypothesized that certain phenotype terms may be harder to phrase in lay person terms due to their nature (e.g. lab values, psychological, etc.)

Given the set disease-phenotype annotations in Monarch as a background data set and enriching on disease terms

1. Analysis of terms with no lay synonyms

Out of the 8312 terms without a lay person synonym, 5215 terms are directly associated with a disease.  Running an enrichment on these terms yielded the disease groups:
nervous system disease
rare neurologic disease
neoplasm (disease)
cardiovascular disease
peripheral nervous system disease
rare genetic neurological disorder
rare tumor
rare neoplastic disease
rare circulatory system disease
rare peripheral neuropathy
peripheral neuropathy

A full list is available in the enriched-disease.tsv file

2. Analysis of terms with lay synonyms

The analysis of 4555 terms with a lay person synonym yielded some borderline significant disease groups:

| Disease | Label | P value (corrected) |
|------|---------|---------|
|MONDO:0019054| congenital limb malformation | 0.06 |
|HP:0040064| Abnormality of limbs  | 0.0733 |
|MONDO:0018455 | dysostosis of genetic origin with limb anomaly as a major feature | 0.0874 |
|MONDO:0018235 | dysostosis with limb anomaly as a major feature | 0.0874 |

The set of diseases that are enriched are related to limb abnormalities, which is not unsurprising since skeletal abnormalities and limb abnormalities have the largest coverage of lay person annotations: 60% and 69% respectively.

