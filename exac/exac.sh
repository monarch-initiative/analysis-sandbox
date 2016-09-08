# Analysis of protein-coding genetic variation in 60,706 humans
# http://www.nature.com/nature/journal/v536/n7616/full/nature19057.html
#
# We have used this catalogue to calculate objective metrics of pathogenicity for sequence variants,
# and to identify genes subject to strong selection against various classes of mutation;
# identifying 3,230 genes with near-complete depletion of predicted protein-truncating variants,
# with 72% of these genes having no currently established human disease phenotype.

# How many of these genes have phenotypes in model orgs?

# Download supplementary files:
# http://www.nature.com/nature/journal/v536/n7616/extref/nature19057-s2.zip
# Convert nature19057-SI Table 13.xlsx, third spreadsheet to tsv >lof_intolerant.tsv


grep '0$' lof_intolerant.tsv > no_phenotype.tsv

# 2311/3230 = 71.5%

cut -f1 no_phenotype.tsv > est_list.txt
sed 's/\.[0-9]\+//' est_list.txt > temp1.txt && mv temp1.txt est_list.txt


virtualenv -p /usr/bin/python3 env
source env/bin/activate
pip3 install requests

./enst2entrez.py --input ./est_list.txt >gene_list.txt 2>ambiguous_ids.txt

# Check ids that didnt map to entrez
grep -v NCBIGene gene_list.txt > temp1.txt
while read line; do
  grep $line no_phenotype.tsv | cut -f2
done <temp1.txt > unmatched_genes.txt

# Try finding entrez id based on gene symbol with mygene
pip3 install mygene
./fetch-gene-ids.py --input unmatched_genes.txt --output mygene_fetch.txt

# 1 input query terms found no hit:
#	['PHF16']

# PHF16 is 9767


# combine files:
grep NCBIGene gene_list.txt | cut -f2 >genelist2.txt
cut -f2 mygene_fetch.txt >genelist1.txt

cat genelist1.txt genelist2.txt >final_gene_list.txt
wc -l final_gene_list.txt
# 2311 final_gene_list.txt

# clean up
rm est_list.txt genelist2.txt genelist1.txt temp1.txt unmatched_genes.txt ambiguous_ids.txt est_list.txt no_phenotype.tsv mygene_fetch.txt

# Run cypher and get stats
python3 ortholog-phenotype-stats.py --input ./final_gene_list.txt --output phenotypes.tsv

# Total orthologs with phenotypes: 2045
# Number of genes with ortholog-phenotypes from 1 taxon: 747
# Number of genes with ortholog-phenotypes from 2 taxa: 913
# Number of genes with ortholog-phenotypes from 3 taxa: 348
# Number of genes with ortholog-phenotypes from 4 taxa: 33
# Number of genes with ortholog-phenotypes from 5 taxa: 3


# 101434 phenotypes.tsv

