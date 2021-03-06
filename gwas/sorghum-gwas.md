HapMap file: https://datacommons.cyverse.org/browse/iplant/home/shared/terraref/genomics/derived_data/bap/resequencing/danforth_center/version1/hapmap/all_combined_genotyped_lines_SNPS_052217.recalibrated_vcftools.filtered.recode.hmp.txt

VCF files: https://datacommons.cyverse.org/browse/iplant/home/shared/terraref/genomics/derived_data/bap/resequencing/danforth_center/version1/gvcf

Phenotype data:
https://drive.google.com/drive/folders/1K4OrHbaDvao7vrN0_V2yc0gD_wmvFXdC with more incoming via https://github.com/genophenoenvo/terraref-datasets/issues/24

https://github.com/genophenoenvo/terraref-datasets/files/4565996/short_format_traits_season_4.txt

Tassel: https://www.maizegenetics.net/tassel

Endpoints:
- SNP to phenotype associations with p-values and effect sizes
- Manhattan plots to determine p-value cutoff
- Phylogenetic tree
- Documentation on this process

Follow up analyses (to become future tickets):
- Look for conservation of gwas snps in closely related species (looking at syntenic regions)
- Look for similar phenotype data in qtl datasets in planteome

-------------------------------------------------

Some notes on the data:

Tassel requires one of the following formats for phenotypes:
Trait format: https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Load/Load#markdown-header-trait-format
Phenotype format: https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Load/Load#markdown-header-numerical-data

-------------------------------------------------

##### Downloading the hapmap and vcf files

Cyverse only downloading 1.1gb/31 then erroring, trying with irods
Irods not supported for OS, dockerizing...

Irods Dockerfile:

    FROM ubuntu:14.04

    VOLUME /data

    RUN apt-get -y update && apt-get install -y wget lsb-core

    RUN wget -qO - https://packages.irods.org/irods-signing-key.asc | apt-key add -
    RUN echo "deb [arch=amd64] https://packages.irods.org/apt/ $(lsb_release -sc) main" | tee /etc/apt/sources.list.d/renci-irods.list
    RUN apt-get -y update && apt-get install -y irods-server irods-database-plugin-postgres

    ENTRYPOINT ["/bin/sh"]

build and run

    docker build . --tag irods

    docker run -it -v /home/kshefchek/git/terraref-datasets:/data irods

set up iicommands per https://wiki.cyverse.org/wiki/display/DS/Setting+Up+iCommands#SettingUpiCommands-co  

```iinit```

hostname: data.cyverse.org  
port: 1247  
zone: iplant  

```cd data```

Relevant iget options:
-K  verify the checksum
-r  recursive

For the hapmap file:

    iget -K /iplant/home/shared/terraref/genomics/derived_data/bap/resequencing/danforth_center/version1/hapmap/all_combined_genotyped_lines_SNPS_052217.recalibrated_vcftools.filtered.recode.hmp.txt

checksum should match adb69f50a3cf1122278f317c3e052e85

    md5sum all_combined_genotyped_lines_SNPS_052217.recalibrated_vcftools.filtered.recode.hmp.txt
    $ adb69f50a3cf1122278f317c3e052e85  all_combined_genotyped_lines_SNPS_052217.recalibrated_vcftools.filtered.recode.hmp.txt

Hapmap header contains extra characters in comparison to trait file (eg IDUE_PI452692 vs PI452692)

    head -1 all_combined_genotyped_lines_SNPS_052217.recalibrated_vcftools.filtered.recode.hmp.txt | sed 's/[A-Z]\{4\}_\(PI[0-9]\{5,6\}\)/\1/g' > all_combined_genotyped_lines_SNPS_052217.hmp.txt
    tail -n +2 all_combined_genotyped_lines_SNPS_052217.recalibrated_vcftools.filtered.recode.hmp.txt >> all_combined_genotyped_lines_SNPS_052217.hmp.txt
    
For the vcf files:

```iget -K -r /iplant/home/shared/terraref/genomics/derived_data/bap/resequencing/danforth_center/version1/gvcf/```

grab some coffee...

##### Merging the VCF files

Unzip, bgzip and tabix index

```
gunzip *
docker pull biocontainers/vcftools:v0.1.16-1-deb_cv1
docker run \
    --volume `pwd`:/data \
    biocontainers/vcftools:v0.1.16-1-deb_cv1 \
    /bin/sh -c 'for F in *.vcf ; do bgzip ${F} ; tabix -f -p vcf ${F}.gz ; done'
```

Download Sorghum reference, https://phytozome.jgi.doe.gov/pz/portal.html#!info?alias=Org_Sbicolor

Note that NCBI download doesn't work, headers need to match (eg >Chr01)

Make fai, and .dict indexes of fasta file
```
docker run \
    --volume `pwd`:/data \
    --user 1001 \
    -d \
    biocontainers/samtools:v1.9-4-deb_cv1 \
    /bin/sh -c 'samtools faidx Sbicolor_313_v3.0.fa'

docker run \
    --volume `pwd`:/data \
    --user 1001 \
    -d \
    biocontainers/picard:v2.3.0_cv2 \
    /bin/sh -c 'picard CreateSequenceDictionary R=Sbicolor_313_v3.0.fa O=Sbicolor_313_v3.0.dict'
```

Merge VCF files using bcftools (takes 2 days and change
```
docker run \
    --volume `pwd`:/data \
    --user 1001 \
    -d \
    biocontainers/bcftools:v1.9-1-deb_cv1 \
    /bin/sh -c 'bcftools merge --merge both --threads 64 --output-type z \
                --gvcf Sbicolor_313_v3.0.fa  *vcf.gz > bcf-merged.vcf.gz'
```
Apply SNP filters, see 
https://docs.terraref.org/protocols/genomic-data#applying-hard-snp-filters-with-gatk-variantfiltration


Run GATK3 GenotypeGVCFs (est 70 hours)

```
docker run \
  -v `pwd`:/gatk/my_data \
  -d \
  broadinstitute/gatk3:3.5-0 \
  /bin/sh -c 'java -Xmx100g -jar GenomeAnalysisTK.jar -T GenotypeGVCFs -R /gatk/my_data/Sbicolor_313_v3.0.fa \
  -V /gatk/my_data/bcf-merged.vcf.gz -o /gatk/my_data/all_combined_Genotyped_lines.vcf --max_alternate_alleles 10'
```

Run GATK3 VariantFiltration

```
docker run \
  -v `pwd`:/gatk/my_data \
  -d \
  broadinstitute/gatk3:3.5-0 \
  /bin/sh -c 'java -Xmx100g -jar GenomeAnalysisTK.jar -T VariantFiltration -R /gatk/my_data/Sbicolor_313_v3.0.fa \
  -V /gatk/my_data/all_combined_Genotyped_lines.vcf -o /gatk/my_data/all_combined_Genotyped_lines_filtered.vcf \
  --filterExpression "QD < 2.0" --filterName "QD" --filterExpression "FS > 60.0" \
  --filterName "FS" --filterExpression "MQ < 40.0" --filterName "MQ" --filterExpression "MQRankSum < -12.5" \
  --filterName "MQRankSum" --filterExpression "ReadPosRankSum < -8.0" --filterName "ReadPosRankSum"'
```

Extract only cultivars with phenotype data (for this run season 4)
```
docker run \
    --volume `pwd`:/data \
    --user 1001 \
    biocontainers/bcftools:v1.9-1-deb_cv1 \
    /bin/sh -c 'bcftools view --output-file all_combined_Genotyped_lines_filtered.season4.vcf.gz --output-type z --samples-file season_4_cultivars.txt  all_combined_Genotyped_lines_filtered.vcf.gz'
```

Keep only biallelic sites, filter and recode
```
docker run \
    --volume `pwd`:/data \
    --user 1001 \
    -d \
    biocontainers/vcftools:v0.1.16-1-deb_cv1 \
    /bin/sh -c 'vcftools --gzvcf all_combined_Genotyped_lines_filtered.season4.vcf.gz --min-alleles 2 --max-alleles 2 \
        --out all_combined_Genotyped_lines.filtered.season4.recode.vcf --max-missing 0.2 --recode'
```

##### Setup and run Tassel 5

    mkdir tassel && cd tassel

Dockerfile, TODO figure out how to accept license

    FROM openjdk:8

    VOLUME /data
    WORKDIR /tassel

    RUN wget https://tassel.bitbucket.io/installer/TASSEL_5_unix.sh -O /tassel/TASSEL_5_unix.sh

    ENTRYPOINT ["/bin/sh"]

Build and run

    docker build . --tag tassel
    docker run -it -v /home/kshefchek/git/terraref-rd2/tassel:/data tassel
    
Inside the container install tassel via:

    sh TASSEL_5_unix.sh

Sort VCF 
```
cd /data
/usr/local/TASSEL5/run_pipeline.pl -Xms75G -Xmx400G -SortGenotypeFilePlugin -inputFile ./all_combined_Genotyped_lines.filtered.season4.recode.vcf -outputFile ./all_combined_Genotyped_lines.filtered.season4.recode.sorted.vcf -fileType VCF
```

Run GWAS
```
cd /data

/usr/local/TASSEL5/run_pipeline.pl -Xms512m -Xmx400g -debug /tassel/debug.txt -fork1 -vcf ./all_combined_Genotyped_lines.filtered.season4.recode.sorted.vcf  -FilterSiteBuilderPlugin -siteMinAlleleFreq 0.01 -endPlugin -fork2  -t ./traits_season_4_vcf.tsv -combine3 -input1 -input2 -intersect -FixedEffectLMPlugin -siteStatsOut -siteStatFile siteStatFile.tsv -biallelicOnly true -endPlugin -export glm_output
```

Bonferroni Correction

    #!/usr/bin/env python3
    import csv
    
    # change denominator depending on number of comparisons
    correct_p = .05/292010910

    output = open('./glm_output1_bfcorrect.tsv', 'w')

    with open('./glm_output1.txt', newline='') as csvfile:
        csvreader = csv.DictReader(csvfile, delimiter='\t')
        csvwriter = csv.DictWriter(output, delimiter='\t', fieldnames=csvreader.fieldnames)
        csvwriter.writeheader()
        for row in csvreader:
            if float(row['p']) <= correct_p:
                csvwriter.writerow(row)


get all significant markers for trimming vcf

    tail -n +2 glm_output1_bfcorrect.tsv | cut -f2 | sort -u > markers-bf.txt

get trait marker combinations to filter allele affects (glm_output2.txt)

    tail -n +2 glm_output1_bfcorrect.tsv | cut -f1,2 > trait-markers.txt


Use join and awk to create filtered glm_output2.txt

    export TMPDIR=/home/kshefchek/git/terraref-datasets/tassel/tmp
    join -t $'\t' -j1 -o1.2,1.3,1.4,1.5,1.6,1.7,1.8 <(<glm_output2.txt awk '{print $1"-"$2" "$0}' | sort -k1,1) <(<trait-markers.txt awk '{print $1"-"$2" "$0}' | sort -k1,1) > glm_output2_bfcorrect.tsv


Generate slim vcf file
```
LANG=en_EN sort -s -k1 markers-bf.txt -o markers-bf.txt
head -11 all_combined_Genotyped_lines.filtered.season4.recode.sorted.vcf >header.txt
tail --lines=+12 all_combined_Genotyped_lines.filtered.season4.recode.sorted.vcf | LANG=en_EN sort -s -k3 > vcf-sorted.txt
LANG=en_EN join -t $'\t' -1 3 -2 1  vcf-sorted.txt markers-bf.txt >vcf-slim.vcf
awk ' BEGIN {OFS = "\t"}; {t = $1; $1 = $2; $2 = $3; $3 = t; print; } ' vcf-slim.vcf | sort --version-sort -k3,1 > vcf-slim.pos.sorted
cat header.txt vcf-slim.pos.sorted >vcf-slim.sorted.vcf
```

Re-sort with tassel
```
docker exec -it ${container_id} /bin/bash
/usr/local/TASSEL5/run_pipeline.pl -fork1 -vcf /data/vcf-slim.sorted.vcf -export -exportType VCF -sortPositions
```

Output files: https://data.monarchinitiative.org/genophenoenvo/tassel5/season4/

Try web app:
https://plants.ensembl.org/Sorghum_bicolor/Tools/VEP

Reference: Sorghum_bicolor_NCBIv3

Constructed cmd

    ./vep --appris --biotype --buffer_size 5000 --check_existing --distance 5000 --mane --sift b --species sorghum_bicolor --symbol --transcript_version --tsl --cache --input_file [input_data] --output_file [output_file]

https://plants.ensembl.org/Sorghum_bicolor/Tools/VEP/Results?tl=1PQiIV1EFQqwCEV7-19684844

Try with up/down distance of 100k:
https://plants.ensembl.org/Sorghum_bicolor/Tools/VEP/Results?tl=nY6uoetDeReR1KgH-19684846


Generate distance matrix and newick trees

```
/usr/local/TASSEL5/run_pipeline.pl -Xms75G -Xmx400G -vcf ./all_combined_Genotyped_lines_vcftools.filtered.recode.sorted.vcf -DistanceMatrixPlugin -endPlugin -export all_cultivars_distance.txt

/usr/local/TASSEL5/run_pipeline.pl -Xms75G -Xmx400G -vcf ./all_combined_Genotyped_lines_vcftools.filtered.recode.sorted.vcf -tree Neighbor -treeSaveDistance false -export tree_newick.nj.txt -exportType Text

/usr/local/TASSEL5/run_pipeline.pl -Xms75G -Xmx400G -vcf ./all_combined_Genotyped_lines_vcftools.filtered.recode.sorted.vcf -tree UPGMA -treeSaveDistance false -export tree_newick.upgma.txt -exportType Text

/usr/local/TASSEL5/run_pipeline.pl -Xms75G -Xmx400G -importGuess all_cultivars_distance.txt -RemoveNaNFromDistanceMatrixPlugin -endPlugin -export all_cultivars_distance_nonan.txt
```

Generate distance matrix for 76-subset

```
docker run \
    --volume `pwd`:/data \
    --user 1001 \
    biocontainers/bcftools:v1.9-1-deb_cv1 \
    /bin/sh -c 'bcftools view --output-file all_combined_Genotyped_lines_filtered_subset76.vcf --output-type v --samples-file subset76.txt all_combined_Genotyped_lines_vcftools.filtered.recode.vcf'
    
/usr/local/TASSEL5/run_pipeline.pl -Xms75G -Xmx400G -SortGenotypeFilePlugin -inputFile ./all_combined_Genotyped_lines_filtered_subset76.vcf -outputFile ./all_combined_Genotyped_lines_filtered_subset76.recode.sorted.vcf -fileType VCF

/usr/local/TASSEL5/run_pipeline.pl -Xms75G -Xmx400G -vcf ./all_combined_Genotyped_lines_filtered_subset76.recode.sorted.vcf -DistanceMatrixPlugin -endPlugin -export subset76_distance_withnan.txt

/usr/local/TASSEL5/run_pipeline.pl -Xms75G -Xmx400G -importGuess subset76_distance_withnan.txt -RemoveNaNFromDistanceMatrixPlugin -endPlugin -export subset76_distance.txt
```

Work on https://github.com/genophenoenvo/terraref-datasets/issues/120

```
wget ftp://ftp.ensemblgenomes.org/pub/plants/current/variation/vep/sorghum_bicolor_vep_48_Sorghum_bicolor_NCBIv3.tar.gz

mv sorghum_bicolor /home/kshefchek/.vep/

./vep --appris --biotype --buffer_size 5000 --check_existing --distance 5000 --mane --cache --genomes --cache_version 48 --species sorghum_bicolor --uniprot --plugin GO --symbol --transcript_version --tsl --input_file /home/kshefchek/git/terraref-rd2/tassel/data/season4-slim.vcf --output_file ./vep_plus_go.tsv
```

##############
Outdated (analysis from hapmap file)

Generate slim vcf file from hapmap file

    sort -s -k1,1 all_combined_genotyped_lines_SNPS_052217.hmp.txt > hmp-sorted.txt
    join -t $'\t' markers-bf.txt hmp-sorted.txt >hapmap-slim.hmp.txt
    sort --version-sort -k1,1 hapmap-slim.hmp.txt > hapmap-slim.pos.sorted
    cat header.txt hapmap-slim.pos.sorted >hapmap-slim.hmp.txt

apparently version-sort is not the right way to sort a vcf, apply -sortPositions and convert to VCF

    /tassel/run_pipeline.pl -fork1 -h /data/hapmap-slim.hmp.txt -export -exportType VCF -sortPositions
