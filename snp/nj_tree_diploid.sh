#!/bin/bash
vcf_file=$1
vcf_tools=/home/sobczm/bin/vcftools/bin
scripts=/home/sobczm/bin/popgen/snp

gzip $vcf_file
zcat "$vcf_file".gz | $vcf_tools/vcf-to-tab > "${vcf_file%.vcf}.tab"
gzip -d "$vcf_file".gz
perl $scripts/vcf_tab_to_fasta_alignment.pl -i "${vcf_file%.vcf}.tab" >"${vcf_file%.vcf}.fasta"
seqtk=/home/sobczm/bin/seqtk/seqtk
#Pick one random allele at heterozygous sites to obtain haploid sequences
$seqtk randbase "${vcf_file%.vcf}.fasta" >"${vcf_file%.vcf}_haploidised.fasta"
Rscript --vanilla $scripts/nj.R "${vcf_file%.vcf}_haploidised.fasta"
