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

#See if use of EMxFE haplotypes improves phasing of RGxHA 
cp -r $datadir/rgxha_vcf_rob_impute2 $datadir/rgxha_vcf_phasing_with_emxfe
cd $datadir/emxfe_vcf
#Split the input VCF by chromosome.
java -jar /home/sobczm/bin/snpEff/SnpSift.jar split emxfe_vcf_test_all_withcoords.vcf
#Convert to Impute2
for number in 1 2 3 4 5 6 7
do
for subgenome in A B C D 
do
vcf2impute2fmt.py emxfe_vcf_test_all_withcoords.${number}${subgenome}.vcf ${number}${subgenome}
done 
done

#Gunzip legend and hap files.
for f in *.legend *.hap; do gzip $f; done

#Copy the Impute2 files over to the Shapeit2 run directory.
cp -r *legend.gz *hap.gz *sample $datadir/rgxha_vcf_phasing_with_emxfe

#Phase
cd ..
export SHAPEIT_EFF_POPN_SIZE=1000
export SHAPEIT_SEED=${RANDOM}${RANDOM}
export SHAPEIT_PRUNE=8
export SHAPEIT_BURN=30
export SHAPEIT_MAIN=20
export SHAPEIT_STATES=500
export SHAPEIT_WINDOW=1.5
#Phase without reference panels - default settings
impute_haplotypes2 rgxha_vcf_phasing_with_emxfe

#See if use of FLxCH haplotypes improves phasing of RGxHA 
#OK, so still some markers missing in EMxFE. Need to remove those markers from
#RGxHA input for now.
cd $datadir/rgxha_vcf_phasing_with_emxfe
for number in 1 2 3 4 5 6 7
do
for subgenome in A B C D 
do
python $scripts/vcfremove.py $datadir/emxfe_vcf/emxfe_vcf_test_all_withcoords.${number}${subgenome}.vcf $datadir/rgxha_vcf_phasing_with_emxfe/rgxha_vcf_test_withcoords.${number}${subgenome}.vcf 
done 
done

#Convert VCF to PED
for number in 1 2 3 4 5 6 7
do
for subgenome in A B C D 
do
python $scripts/vcf2ped.py rgxha_vcf_test_withcoords.${number}${subgenome}_filtered.vcf ${number}${subgenome}.pmap random
done
done
rm *pmap *rec
#Substitute rec2 for rec and pmap2 for pmap.
for r in *.rec2 
do
mv $r ${r%.rec2}.rec
done
for p in *.pmap2 
do
mv $p ${p%.pmap2}.pmap
done

#Remove old impute2 files
rm *.hap* *.legend* *.sample
for number in 1 2 3 4 5 6 7
do
for subgenome in A B C D 
do
vcf2impute2fmt.py $datadir/emxfe_vcf/emxfe_vcf_test_all_withcoords.${number}${subgenome}.vcf ${number}${subgenome}
done 
done

#Gunzip legend and hap files.
for f in *.legend *.hap; do gzip $f; done

#Convert map-based haplotypes in the VCF format to alleles:
for number in 1 2 3 4 5 6 7
do
for subgenome in A B C D 
do
python $scripts/vcf2AB.py rgxha_vcf_test_withcoords.${number}${subgenome}_filtered.vcf $number$subgenome
done
done


#Phase
cd ..
export SHAPEIT_EFF_POPN_SIZE=1000
export SHAPEIT_SEED=${RANDOM}${RANDOM}
export SHAPEIT_PRUNE=8
export SHAPEIT_BURN=30
export SHAPEIT_MAIN=20
export SHAPEIT_STATES=500
export SHAPEIT_WINDOW=1.5
#Phase without reference panels - default settings
impute_haplotypes2 rgxha_vcf_phasing_with_emxfe

for hapfile in $datadir/rgxha_vcf_phasing_with_emxfe/*.haps
do
    lg=$(basename ${hapfile} .phased.haps)
    conv_phased2AB.py  ${hapfile}  > ${hapfile/.phased.haps/_imputed_AB}   
    compare_haplotypes.py ${hapfile/.phased.haps/_imputed_AB} $datadir/rgxha_vcf_phasing_with_emxfe/${lg}_original_vcf  > ${hapfile/.phased.haps/_out}
done
Rscript --vanilla $scripts/plots_phasing.R RGxHA_EMxFE_ref $datadir/rgxha_test_ms_shapeit $datadir/rgxha_vcf_phasing_with_emxfe
