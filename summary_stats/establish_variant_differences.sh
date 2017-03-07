#!/bin/bash
scripts=/home/sobczm/bin/popgen/summary_stats
#After calling SNPs and stuctural variants, it may be useful to contrast variant population
#frequencies among different groups of individuals (populations) to e.g. further zoom in on
#candidate effector genes likely involved in resistance in different cultivars.

### V. inequalis SNPs
##F irst, create a cut-down VCF file containing only individuals of interest 
## (Worcester and Bramley populations) and then filter to remove SNPs with too many ##missing genotypes.
input=/home/sobczm/popgen/other/passey
cd $input

vcflib=/home/sobczm/bin/vcflib/bin
$vcflib/vcfremovesamples Ash_farm_172_pacbio_contigs_unmasked.vcf 049 083 096 097 098 101 106 118 119 199 >Ash_farm_172_pacbio_contigs_unmasked_bw.vcf 

vcftools=/home/sobczm/bin/vcftools/bin
$vcftools/vcftools --vcf Ash_farm_172_pacbio_contigs_unmasked_bw.vcf  --max-missing 0.95 --recode --out Ash_farm_172_pacbio_contigs_unmasked_bw_filtered

## Important: check script options below, in order to use the correct ones in
## a given analysis (e.g. "ply" argument for ploidy).
python $scripts/vcf_find_difference_pop.py --vcf Ash_farm_172_pacbio_contigs_unmasked_bw_filtered.recode.vcf --out Ash_farm_172_pacbio_contigs_unmasked_bw_filtered_fixed.vcf --ply 1 --pop1 202,,182,,173,,190,,172,,197,,196 --pop2 057,,024,,030,,007,,025,,044 --thr 0.95

###V. inequalis structural variants
input=/home/sobczm/popgen/snp/sv_calling
python $scripts/vcf_find_difference_pop.py --vcf $input/vinequalis/vinequalis_struc_variants.vcf --out $input/vinequalis/vinequalis_struc_variants_fixed.vcf --ply 1 --pop1 202,,182,,173,,190,,172,,197,,196 --pop2 057,,024,,030,,007,,025,,044 --thr 0.95
