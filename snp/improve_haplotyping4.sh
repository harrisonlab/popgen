#Follow on from Rob's analysis and test different run settings for shapeit haplotyping.
datadir=/home/sobczm/popgen/snp/snp_chip/haplotypes/comparison/rgxha_tweaking_original
#based on 130_compare_map2imputed.sh scripts
SCRIPT_PATH=/home/sobczm/bin/test_haplotyping
export PATH=${SCRIPT_PATH}:${PATH}
export PATH=/home/vicker/git_repos/crosslink/bin:${PATH}
source ${SCRIPT_PATH}/haplotyping_all_samples_funcs.sh
export NODES="-N1,2,5,6,7,8,9,10,11"
export MAXJOBS="-L30"
scripts=/home/sobczm/bin/popgen/snp
export PATH=${PATH}:/home/vicker/programs/vcftools-0.1.15/bin
export PERL5LIB=${PERL5LIB}:/home/vicker/programs/vcftools-0.1.15/src/perl

#Use VCF input files with map-based haplotypes to allow comparison across different mapping populations.
cd $datadir/rgxha_vcf
#Split the input VCF by chromosome.
java -jar /home/sobczm/bin/snpEff/SnpSift.jar split rgxha_vcf_test_withcoords.vcf
#Copy *pmap and *rec files needed for Shapeit run. Use the new version which features no missing
#variants in the rec file.
cp -r $datadir/rgxha_fixed/rgxha_test_shapeit/*.pmap $datadir/rgxha_fixed/rgxha_test_shapeit/*.rec ./
#Convert VCF to PED
for number in 1 2 3 4 5 6 7
do
for subgenome in A B C D 
do
python $scripts/vcf2ped.py rgxha_vcf_test_withcoords.${number}${subgenome}.vcf ${number}${subgenome}.pmap random
done
done
rm *pmap *rec
#Substitute rec2 for rec and pmap2 for pmap.
for r in $datadir/rgxha_vcf/*.rec2 
do
mv $r ${r%.rec2}.rec
done
for p in $datadir/rgxha_vcf/*.pmap2 
do
mv $p ${p%.pmap2}.pmap
done

#Convert map-based haplotypes in the VCF format to alleles:
for number in 1 2 3 4 5 6 7
do
for subgenome in A B C D 
do
python $scripts/vcf2AB.py rgxha_vcf_test_withcoords.$number$subgenome.vcf $number$subgenome
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

cp -r $datadir/rgxha_vcf $datadir/rgxha_vcf_phasing

#First phase using best settings.
export SHAPEIT_EFF_POPN_SIZE=1000
export SHAPEIT_SEED=${RANDOM}${RANDOM}
export SHAPEIT_PRUNE=8
export SHAPEIT_BURN=30
export SHAPEIT_MAIN=20
export SHAPEIT_STATES=500
export SHAPEIT_WINDOW=1.5
#Phase without reference panels - default settings
impute_haplotypes rgxha_vcf_phasing
cp -r $datadir/rgxha_vcf_phasing $datadir/rgxha_vcf_phasing_with_ref
#Phase with reference panel
impute_haplotypes2 rgxha_vcf_phasing_with_ref

for hapfile in $datadir/rgxha_vcf_phasing/*.haps
do
    lg=$(basename ${hapfile} .phased.haps)
    conv_phased2AB.py  ${hapfile}  > ${hapfile/.phased.haps/_imputed_AB}   
    compare_haplotypes.py ${hapfile/.phased.haps/_imputed_AB} $datadir/rgxha_vcf_phasing/${lg}_original_vcf  > ${hapfile/.phased.haps/_out}
done

for hapfile in $datadir/rgxha_vcf_phasing_with_ref/*.haps
do
    lg=$(basename ${hapfile} .phased.haps)
    conv_phased2AB.py  ${hapfile}  > ${hapfile/.phased.haps/_imputed_AB}   
    compare_haplotypes.py ${hapfile/.phased.haps/_imputed_AB} $datadir/rgxha_vcf_phasing_with_ref/${lg}_original_vcf  > ${hapfile/.phased.haps/_out}
done

Rscript --vanilla $scripts/plots_phasing.R RGxHA_vcf $datadir/rgxha_test_ms_shapeit $datadir/rgxha_vcf_phasing
Rscript --vanilla $scripts/plots_phasing.R RGxHA_vcf_with_ref $datadir/rgxha_test_ms_shapeit $datadir/rgxha_vcf_phasing_with_ref

#Use bcftools to convert from VCF to Impute2.
#bcftools=/home/sobczm/bin/bcftools-1.5/bin/bcftools
#$bcftools convert --gensample bcf rgxha_vcf_test_withcoords.1A.vcf

#Confirm that the VCF file and converted hapfiles do match.
cd $datadir/rgxha_vcf_phasing_with_ref_test
for number in 1 2 3 4 5 6 7
do
for subgenome in A B C D 
do
python $scripts/compare_impute2vcf.py rgxha_vcf_test_withcoords.${number}${subgenome}.vcf ${number}${subgenome}.hap
done 
done

#Confirm that the VCF file and the previously converted (and definietely working) hapfiles do match.
#Generate the Impute2 files with the matching sample names.
cp -r $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_rgxha_impute2
list_inds=/home/vicker/octoploid_mapping/consensus_map4/popn_RGxHA/map/conf/popn_order
for number in 1 2 3 4 5 6 7
do
for subgenome in A B C D 
do
python $scripts/convert_to_impute2_general.py $datadir/rgxha_test_ms_${number}${subgenome}_orig_AB $datadir/rgxha_test_ms_shapeit_rgxha_impute2/${number}${subgenome}.pmap $list_inds rgxha
done
done

#Convert B to Cs in the legend file.
for leg in $datadir/rgxha_test_ms_shapeit_rgxha_impute2/*.legend
do
sed -i 's/B/C/g' $leg
done

for number in 1 2 3 4 5 6 7
do
for subgenome in A B C D 
do
python $scripts/compare_impute2vcf.py $datadir/rgxha_vcf_phasing_with_ref_test/rgxha_vcf_test_withcoords.${number}${subgenome}.vcf $datadir/rgxha_test_ms_shapeit_rgxha_impute2/${number}${subgenome}.hap
done 
done