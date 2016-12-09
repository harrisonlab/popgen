#!/bin/bash
#Functionally annotates all the variants in the given VCF file as well as
#creates a custom subset of different types of polymorphisms (genic, coding,
#non-synonymous, synonymous, as well as four fold degenerate sites) useful for
#downstream analyses of population demography and selection.

#Input:
#First argument: VCF file to be annotated. Genome database must be pre-built.
#(see build_genome_database.sh script)
#Second argument: Name of the genome database to be used for generating the VCF file annotation
#(see build_genome_database.sh script)

#Output:
#Various VCF files with suffix "_annotated"

vcf=$1
genome_name=$2

snpeff=/home/sobczm/bin/snpEff
scripts=/home/sobczm/bin/popgen

filename=$(basename "$vcf")
java -Xmx4g -jar $snpeff/snpEff.jar -v -ud 0 $genome_name $vcf > ${filename%.vcf}_annotated.vcf
mv snpEff_genes.txt snpEff_genes_${filename%.vcf}.txt
mv snpEff_summary.html  snpEff_summary__${filename%.vcf}.html

#Create subsamples of SNPs containing those in a given category

#genic (includes 5', 3' UTRs)
java -jar $snpeff/SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant') || (ANN[*].EFFECT has 'nonsense_variant') || (ANN[*].EFFECT has 'synonymous_variant') || (ANN[*].EFFECT has 'intron_variant') || (ANN[*].EFFECT has '5_prime_UTR_variant') || (ANN[*].EFFECT has '3_prime_UTR_variant')" ${filename%.vcf}_annotated.vcf > ${filename%.vcf}_gene.vcf
#coding
java -jar $snpeff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'synonymous_variant')" ${filename%.vcf}_annotated.vcf > ${filename%.vcf}_coding.vcf
#non-synonymous
java -jar $snpeff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant')" ${filename%.vcf}_annotated.vcf > ${filename%.vcf}_nonsyn.vcf
#synonymous
java -jar $snpeff/SnpSift.jar filter "(ANN[0].EFFECT has 'synonymous_variant')" ${filename%.vcf}_annotated.vcf > ${filename%.vcf}_syn.vcf
#Four-fold degenrate sites (output file suffix: 4fd)
python $scripts/summary_stats/parse_snpeff_synonymous.py ${filename%.vcf}_syn.vcf
