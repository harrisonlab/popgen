#!/bin/bash
#Input: VCF file containing bi-allelic SNPs. Output: a neighbour joining tree, suffix: _nj.nwk, and PDF file with the tree, suffix: _nj.pdf
vcf_file=$1
vcf_tools=/home/sobczm/bin/vcftools/bin
scripts=/home/sobczm/bin/popgen/snp

gzip $vcf_file
zcat "$vcf_file".gz | $vcf_tools/vcf-to-tab > "${vcf_file%.vcf}.tab"
gzip -d "$vcf_file".gz
perl $scripts/vcf_tab_to_fasta_alignment.pl -i "${vcf_file%.vcf}.tab" >"${vcf_file%.vcf}.fasta"
Rscript --vanilla $scripts/nj.R "${vcf_file%.vcf}.fasta"
