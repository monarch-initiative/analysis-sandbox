HapMap file: https://datacommons.cyverse.org/browse/iplant/home/shared/terraref/genomics/derived_data/bap/resequencing/danforth_center/version1/hapmap/all_combined_genotyped_lines_SNPS_052217.recalibrated_vcftools.filtered.recode.hmp.txt

Phenotype data:
https://drive.google.com/drive/folders/1K4OrHbaDvao7vrN0_V2yc0gD_wmvFXdC with more incoming via https://github.com/genophenoenvo/terraref-datasets/issues/24

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

- Tall format contains time series data which I'm fairly certain tassel does not support.  For now we'll take the final measurement
- Some measurements have the same time, eg leaf temperature.  @MagicMilly doesn't think this is a resolution issue but rather there were multiple iterative measurements taken at once.
- Some measurements were only taken for a small subset of cultivars, eg leaf temperature was only measured for 12 culitvars.  Should we add some filter for traits with < N cultivars (50?)
- Cultivars in the hapmap file are missing from the tall_format and vice versa.  67 cultivars from the hapmap file do not have phenotypes, bringing the cultivar sample size down to 295.
- Some measurements do not look like traits, eg lodging_present, harvest_lodging_rating.  Should these be treated as factors and covariates or removed?

Tassel requires one of the following formats for phenotypes:
Trait format: https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Load/Load#markdown-header-trait-format
Phenotype format: https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Load/Load#markdown-header-numerical-data

-------------------------------------------------

The gist of all this

![meme](https://i.imgflip.com/3p1hfz.jpg)

##### Downloading the hapmap file

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

    iinit

    iget -K /iplant/home/shared/terraref/genomics/derived_data/bap/resequencing/danforth_center/version1/hapmap/all_combined_genotyped_lines_SNPS_052217.recalibrated_vcftools.filtered.recode.hmp.txt

checksum should match adb69f50a3cf1122278f317c3e052e85

    md5sum all_combined_genotyped_lines_SNPS_052217.recalibrated_vcftools.filtered.recode.hmp.txt
    $ adb69f50a3cf1122278f317c3e052e85  all_combined_genotyped_lines_SNPS_052217.recalibrated_vcftools.filtered.recode.hmp.txt

Hapmap header contains extra characters in comparison to trait file (eg IDUE_PI452692 vs PI452692)

    head -1 all_combined_genotyped_lines_SNPS_052217.recalibrated_vcftools.filtered.recode.hmp.txt | sed 's/[A-Z]\{4\}_\(PI[0-9]\{5,6\}\)/\1/g' > all_combined_genotyped_lines_SNPS_052217.hmp.txt
    tail -n +2 all_combined_genotyped_lines_SNPS_052217.recalibrated_vcftools.filtered.recode.hmp.txt >> all_combined_genotyped_lines_SNPS_052217.hmp.txt

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
    docker run -it -v /home/kshefchek/git/terraref-datasets/tassel:/data tassel

    /tassel/run_pipeline.pl -Xms512m -Xmx425g -debug /tassel/debug.txt -fork1 -h ./all_combined_genotyped_lines_SNPS_052217.hmp.txt  -FilterSiteBuilderPlugin -siteMinAlleleFreq 0.01 -endPlugin -fork2  -t ./short_format_traits_season_4.tsv -combine3 -input1 -input2 -intersect -FixedEffectLMPlugin -siteStatsOut -siteStatFile siteStatFile.tsv -endPlugin -export glm_output


Bonferroni Correction

    #!/usr/bin/env python3
    import csv

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

    sort -s -k1,1 all_combined_genotyped_lines_SNPS_052217.hmp.txt > hmp-sorted.txt
    join -t $'\t' markers-bf.txt hmp-sorted.txt >hapmap-slim.hmp.txt
    sort --version-sort -k1,1 hapmap-slim.hmp.txt > hapmap-slim.pos.sorted
    cat header.txt hapmap-slim.pos.sorted >hapmap-slim.hmp.txt

Looks sorted to me but Tassel doesn't think so, apply -sortPositions and convert to VCF

    /tassel/run_pipeline.pl -fork1 -h /data/hapmap-slim.hmp.txt -export -exportType VCF -sortPositions


vcf file: https://data.monarchinitiative.org/tassel5/hapmap-slim.vcf

Try web app:
https://plants.ensembl.org/Sorghum_bicolor/Tools/VEP

Reference: Sorghum_bicolor_NCBIv3

Constructed cmd

    ./vep --appris --biotype --buffer_size 5000 --check_existing --distance 5000 --mane --sift b --species sorghum_bicolor --symbol --transcript_version --tsl --cache --input_file [input_data] --output_file [output_file]

https://plants.ensembl.org/Sorghum_bicolor/Tools/VEP/Results?tl=1PQiIV1EFQqwCEV7-19684844

Try with up/down distance of 100k:
https://plants.ensembl.org/Sorghum_bicolor/Tools/VEP/Results?tl=nY6uoetDeReR1KgH-19684846
