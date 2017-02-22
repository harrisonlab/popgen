#!/bin/bash
scripts=/home/sobczm/bin/popgen/summary_stats
input_hap=/home/sobczm/popgen/other/passey/Ash_farm_172_pacbio_contigs_unmasked_filtered.vcf
input_dip=/home/sobczm/popgen/other/phytophthora/Ash_farm_172_pacbio_contigs_unmasked_no118.vcf

Rscript --vanilla $scripts/run_pcadapt.R $input_hap 1
Rscript --vanilla $scripts/run_pcadapt.R $input_dip 2