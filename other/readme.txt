Question: What genes are haploinsufficient in Mice?

First pass will look at heterozygous simple genotypes (only one variant)
Second pass could include hemizygous variants

# path is intrinsic->gvc/vslc->zygosity
MATCH (genotype:genotype)-[hasPhenotype:RO:0002200]->(phenotype:Phenotype)
MATCH (genotype:genotype)-[:BFO:0000051!]->(background)-[tax:RO:0002162]->(mouse{iri:'http://purl.obolibrary.org/obo/NCBITaxon_10090'})
MATCH (genotype)-[:BFO:0000051!*1]->(vslc)-[:type]->(vslc_type{iri:'http://purl.obolibrary.org/obo/GENO_0000030'})
MATCH (vslc)-[:GENO:0000608]->(zygosity{iri:'http://purl.obolibrary.org/obo/GENO_0000135'})
OPTIONAL MATCH (vslc)-[:BFO:0000051!]->(allele)-[:GENO:0000408]->(gene)
RETURN DISTINCT genotype.iri as genotype, zygosity.label as zygosity, allele.iri as allele, gene.iri as gene

# Ideally use json but service cannot process nulls
curl -g 'https://scigraph-data.monarchinitiative.org/scigraph/cypher/execute?cypherQuery=MATCH+(genotype%3Agenotype)-[hasPhenotype%3ARO%3A0002200]-%3E(phenotype%3APhenotype)+MATCH+(genotype%3Agenotype)-[%3ABFO%3A0000051!]-%3E(background)-[tax%3ARO%3A0002162]-%3E(mouse{iri%3A%27http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FNCBITaxon_10090%27})+MATCH+(genotype)-[%3ABFO%3A0000051!*1]-%3E(vslc)-[%3Atype]-%3E(vslc_type{iri%3A%27http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FGENO_0000030%27})+MATCH+(vslc)-[%3AGENO%3A0000608]-%3E(zygosity{iri%3A%27http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FGENO_0000135%27})+OPTIONAL+MATCH+(vslc)-[%3ABFO%3A0000051!]-%3E(allele)-[%3AGENO%3A0000408]-%3E(gene)+RETURN+DISTINCT+genotype.iri+as+genotype%2C+zygosity.label+as+zygosity%2C+allele.iri+as+allele%2C+gene.iri+as+gene&limit=1000000' | grep '|' | sed -e 's/|\s//' | sed -e 's/\s\+|\s\+/\t/g' | sed -e 's/\s\+|$//g' | sed -e 's/"//g' >het_genotypes.tsv


Some Zygosity Classes:
GENO_0000135: heterozygous
GENO_0000606: hemizygous insertion-linked
GENO_0000604: hemizygous X-linked
GENO_0000136: homozygous

Lethal parent classes:
MP:0011400 lethality complete penetrance
MP:0010831 lethality, incomplete penetrance 
MP:0008569 lethality at weaning
MP:0010770 preweaning lethality