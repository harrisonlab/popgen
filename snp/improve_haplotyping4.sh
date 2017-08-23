#Follow on from Rob's analysis and test different run settings for shapeit haplotyping.
datadir=/home/sobczm/popgen/snp/snp_chip/haplotypes/comparison/rgxha_tweaking_original
#based on 130_compare_map2imputed.sh scripts
SCRIPT_PATH=/home/sobczm/bin/test_haplotyping
export PATH=${SCRIPT_PATH}:${PATH}
export PATH=/home/vicker/git_repos/crosslink/bin:${PATH}
source ${SCRIPT_PATH}/haplotyping_all_samples_funcs.sh
export NODES="-N2,5,6,7,8,9,10"
export MAXJOBS="-L30"
scripts=/home/sobczm/bin/popgen/snp
export PATH=${PATH}:/home/vicker/programs/vcftools-0.1.15/bin
export PERL5LIB=${PERL5LIB}:/home/vicker/programs/vcftools-0.1.15/src/perl

#Use VCF input files with map-based haplotypes to allow comparison across different mapping populations.
cd $datadir/rgxha_vcf
#Convert map-based haplotypes in the VCF format to ABs:
python $scripts/vcf2AB.py rgxha_vcf_test_withcoords.vcf
#Split the input VCF by chromosome.
java -jar /home/sobczm/bin/snpEff/SnpSift.jar split rgxha_vcf_test_withcoords.vcf
#Copy *pmap and *rec files needed for Shapeit run.
cp -r $datadir/rgxha_test_ms_shapeit/*.pmap $datadir/rgxha_test_ms_shapeit/*.rec ./
#Convert VCF to PED
for number in 1 2 3 4 5 6 7
do
for subgenome in A B C D 
do
python $scripts/vcf2ped.py rgxha_vcf_test_withcoords.${number}${subgenome}.vcf ${number}${subgenome}.pmap random
done
done
#Convert VCF to Impute2 ref panel format.
for number in 1 2 3 4 5 6 7
do
for subgenome in A B C D 
do
python $scripts/vcf2haps.py rgxha_vcf_test_withcoords.${number}${subgenome}.vcf ${number}${subgenome}.pmap rgxha
done
done