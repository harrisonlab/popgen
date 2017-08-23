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

cd $datadir
mkdir rgxha_vcf
cd rgxha_vcf
#convert crosslink files into vcf format
name=rgxha_vcf_test
mapfile=/home/vicker/fragaria_x_ananassa/octoseq/maps/consensus/vesca2consensus_map_noambiguous_2017-08-22.csv
locfiledir=/home/vicker/octoploid_mapping/consensus_map4/popn_RGxHA/map/popnordergrps
maternal="Redgauntlet"
paternal="Hapil"
orderfile=/home/vicker/octoploid_mapping/consensus_map4/popn_RGxHA/map/conf/popn_order
affycalls=/home/vicker/octoploid_mapping/consensus_map4/popn_RGxHA/mhr_phr_nmh_merged.tsv
orig_parental_calls=/home/vicker/octoploid_mapping/consensus_map4/popn_RGxHA/map/parental_genotypes.tsv

#collect all phased loc files together
cat ${locfiledir}/??.loc > all.loc

#create an ordered list of sample names, including maternal and paternal parents as first two items
echo ${maternal}  >  sample_names
echo ${paternal}  >> sample_names
cat ${orderfile}  >> sample_names

mapbased_haplotypes2vcf.py   all.loc   ${orig_parental_calls}   sample_names >   ${name}_nocoords.vcf

#convert affy calls sample names into same format as crosslink sample names
echo -ne 'probeset_id\t' > affy.tsv
head -n 1 ${affycalls} | tr '\t' '\n' | tail -n +2 | tr '._' ' ' | cut -d' ' -f4 | tr '\n' '\t' >> affy.tsv
echo >> affy.tsv

#append just the monomorphic calls
tail -n +2 ${affycalls} | grep 'MHR' | probe_to_snp.py >> affy.tsv

#convert unphased affy genotypes into "phased" monomorphic vcf
unphased_affy2vcf.py affy.tsv | sed 's/\//|/g'> ${name}_monomorphic.vcf

#reorder columns so that the two files match then concat the markers
vcf-shuffle-cols -t ${name}_nocoords.vcf ${name}_monomorphic.vcf > ${name}_monoordered.vcf
vcf-concat ${name}_monoordered.vcf ${name}_nocoords.vcf > ${name}_all_nocoords.vcf

#add map positions
specify_snp_coords.py ${name}_all_nocoords.vcf ${mapfile} > ${name}_all_withcoords.vcf

#For EMxFE
cd $datadir
mkdir emxfe_vcf
cd emxfe_vcf
#convert crosslink files into vcf format
name=emxfe_vcf_test
mapfile=/home/vicker/fragaria_x_ananassa/octoseq/maps/consensus/vesca2consensus_map_noambiguous_2017-08-22.csv
locfiledir=/home/vicker/octoploid_mapping/consensus_map4/popn_EMxFE/map/popnordergrps
maternal="Emily"
paternal="Fenella"
orderfile=/home/vicker/octoploid_mapping/consensus_map4/popn_EMxFE/map/conf/popn_order
affycalls=/home/vicker/octoploid_mapping/consensus_map4/popn_EMxFE/mhr_phr_nmh_merged.tsv
orig_parental_calls=/home/vicker/octoploid_mapping/consensus_map4/popn_EMxFE/map/parental_genotypes.tsv
cat ${locfiledir}/??.loc > all.loc

#create an ordered list of sample names, including maternal and paternal parents as first two items
echo ${maternal}  >  sample_names
echo ${paternal}  >> sample_names
cat ${orderfile}  >> sample_names
mapbased_haplotypes2vcf.py   all.loc   ${orig_parental_calls}   sample_names >   ${name}_nocoords.vcf

#convert affy calls sample names into same format as crosslink sample names
echo -ne 'probeset_id\t' > affy.tsv
head -n 1 ${affycalls} | tr '\t' '\n' | tail -n +2 | tr '._' ' ' | cut -d' ' -f4 | tr '\n' '\t' >> affy.tsv
echo >> affy.tsv

#append just the monomorphic calls
tail -n +2 ${affycalls} | grep 'MHR' | probe_to_snp.py >> affy.tsv

#convert unphased affy genotypes into "phased" monomorphic vcf
unphased_affy2vcf.py affy.tsv | sed 's/\//|/g'> ${name}_monomorphic.vcf

#reorder columns so that the two files match then concat the markers
vcf-shuffle-cols -t ${name}_nocoords.vcf ${name}_monomorphic.vcf > ${name}_monoordered.vcf
vcf-concat ${name}_monoordered.vcf ${name}_nocoords.vcf > ${name}_all_nocoords.vcf

#add map positions
specify_snp_coords.py ${name}_all_nocoords.vcf ${mapfile} > ${name}_all_withcoords.vcf

#For FLxCH
cd $datadir
mkdir flxch_vcf
cd flxch_vcf
#convert crosslink files into vcf format
name=flxch_vcf_test
mapfile=/home/vicker/fragaria_x_ananassa/octoseq/maps/consensus/vesca2consensus_map_noambiguous_2017-08-22.csv
locfiledir=/home/vicker/octoploid_mapping/consensus_map4/popn_FLxCH/map/popnordergrps
maternal="Flamenco"
paternal="Chandler"
orderfile=/home/vicker/octoploid_mapping/consensus_map4/popn_FLxCH/map/conf/popn_order
affycalls=/home/vicker/octoploid_mapping/consensus_map4/popn_FLxCH/mhr_phr_nmh_merged.tsv
orig_parental_calls=/home/vicker/octoploid_mapping/consensus_map4/popn_FLxCH/map/parental_genotypes.tsv
cat ${locfiledir}/??.loc > all.loc

#create an ordered list of sample names, including maternal and paternal parents as first two items
echo ${maternal}  >  sample_names
echo ${paternal}  >> sample_names
cat ${orderfile}  >> sample_names
mapbased_haplotypes2vcf.py   all.loc   ${orig_parental_calls}   sample_names >   ${name}_nocoords.vcf

#convert affy calls sample names into same format as crosslink sample names
echo -ne 'probeset_id\t' > affy.tsv
head -n 1 ${affycalls} | tr '\t' '\n' | tail -n +2 | tr '._' ' ' | cut -d' ' -f4 | tr '\n' '\t' >> affy.tsv
echo >> affy.tsv

#append just the monomorphic calls
tail -n +2 ${affycalls} | grep 'MHR' | probe_to_snp.py >> affy.tsv

#convert unphased affy genotypes into "phased" monomorphic vcf
unphased_affy2vcf.py affy.tsv | sed 's/\//|/g'> ${name}_monomorphic.vcf

#reorder columns so that the two files match then concat the markers
vcf-shuffle-cols -t ${name}_nocoords.vcf ${name}_monomorphic.vcf > ${name}_monoordered.vcf
vcf-concat ${name}_monoordered.vcf ${name}_nocoords.vcf > ${name}_all_nocoords.vcf

#add map positions
specify_snp_coords.py ${name}_all_nocoords.vcf ${mapfile} > ${name}_all_withcoords.vcf
specify_snp_coords.py ${name}_nocoords.vcf ${mapfile} > ${name}_withcoords.vcf