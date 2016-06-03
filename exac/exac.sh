# (From email)
# ExAc paper http://biorxiv.org/content/early/2015/10/30/030338
# in particular (page 10), they found 79% of loss-of-function-intolerant genes have no known associated disease (~2500).  
# we should add analysis of what model organisms can tell us about them!

# Download file 
curl 'https://docs.google.com/spreadsheets/d/18isx46crTMeeDif05BwBD37veko8Hyrvbf9WuACoubo/export?gid=0&format=tsv' > all_exac.tsv

# file all_exac.tsv
# all_exac.tsv: ASCII text, with CRLF line terminators

dos2unix all_exac.tsv

# Just get rows that do not have a path to clinvar
grep '0$' all_exac.tsv > no_clinvar.tsv

wc -l no_clinvar.tsv 
#2557 no_clinvar.tsv

# Get transcript IDs, remove exon number
cut -f1 no_clinvar.tsv > est_list.txt
sed 's/\.[0-9]\+//' est_list.txt > temp1.txt && mv temp1.txt est_list.txt

virtualenv -p /usr/bin/python3 env
source env/bin/activate
pip3 install requests

./enst2entrez.py --input ~/exac/est_list.txt >gene_list.txt 2>ambiguous_ids.txt

wc -l gene_list.txt
# 2569 gene_list.txt

wc -l ambiguous_ids.txt
# 12 ambiguous_ids.txt

# Check ambiguous IDs, use monarch search to see what's in our db to disambiguate
# removed following lines:
# ENST00000367242	NCBIGene:653659	TMEM183A
# ENST00000545588	NCBIGene:441155	ZC3H11A
# ENST00000376656	NCBIGene:202658	TRIM39
# ENST00000265437	NCBIGene:93655	ST7
# ENST00000373647	NCBIGene:2844	RABGAP1
# ENST00000433855	NCBIGene:3742	GALNT8
# ENST00000220429	NCBIGene:55056	GABPB1
# ENST00000328114	NCBIGene:284194	LGALS9C
# ENST00000337665	NCBIGene:100505585	ARHGEF1
# ENST00000291552	NCBIGene:102724594	U2AF1
# ENST00000291565	NCBIGene:105372824	PDXK
# ENST00000396832	NCBIGene:102800317	CSNK1E


# Check ids that didnt map to entrez
grep -v NCBIGene gene_list.txt > temp1.txt
while read line; do
  grep $line no_clinvar.tsv | cut -f2
done <temp1.txt > unmatched_genes.txt

# Try finding entrez id based on gene symbol with mygene
pip3 install mygene
./fetch-gene-ids.py --input unmatched_genes.txt --output mygene_fetch.txt

# 1 input query terms found dup hits:
#	[('FAT3', 2)]
# 1 input query terms found no hit:
#	['PHF16']

# FAT3 is right, PHF16 is 9767

# combine files:
grep NCBIGene gene_list.txt | cut -f2 >genelist2.txt
cut -f2 mygene_fetch.txt >genelist1.txt

cat genelist1.txt genelist2.txt >final_gene_list.txt
wc -l final_gene_list.txt
# 2557 final_gene_list.txt

# clean up
rm est_list.txt genelist2.txt genelist1.txt temp1.txt unmatched_genes.txt ambiguous_ids.txt est_list.txt no_clinvar.tsv mygene_fetch.txt

# Run cypher and get stats
python3 ortholog-phenotype-stats.py --input ./final_gene_list.txt --output phenotypes.tsv
# Total orthologs with phenotypes: 1551
# Number of genes with ortholog-phenotypes from 1 taxon: 1386
# Number of genes with ortholog-phenotypes from 2 taxa: 163
# Number of genes with ortholog-phenotypes from 3 taxa: 2
# Number of genes with ortholog-phenotypes from 4 taxa: 0

# 38972 gene-phenotype assocs

# However, as someone commented biorxip.org, if they are only considering clinvar and hgmd for gene-disease associations, they are likely missing a lot gene to disease assocs.  For example, they include the huntingtin gene (HTT), and JAK2 (which has 26 disease assocs in our db) as not have having a known associated disease

# Going back to check for gene - disease associations in our database
./get-gene-disease.py --input ./final_gene_list.txt --output gene-disease.tsv >gene-with-disease.list

wc -l gene-with-disease.list
# 934 gene-with-disease.list

# So total is 3229, genes without disease is 2557 - 934, 1623/3230 = 50.23%
# This is consistent with the biorxip.org commenter who said 850 genes were in disgenet
# The headline should read: 50% of loss-of-function-intolerant genes have no known associated disease

# Get the genes without an associated disease in monarch

grep -Fxv -f gene-with-disease.list final_gene_list.txt > genes-minus-monarch.txt

#Rerun script
python3 ortholog-phenotype-stats.py --input genes-minus-monarch.txt --output phenotypes.tsv

#Total orthologs with phenotypes: 836
#Number of genes with othorlog-phenotypes from 1 taxon: 772
#Number of genes with othorlog-phenotypes from 2 taxa: 63
#Number of genes with othorlog-phenotypes from 3 taxa: 1
#Number of genes with othorlog-phenotypes from 4 taxa: 0
#14012 gene-phenotype assocs

# % coverage with models = 76%








