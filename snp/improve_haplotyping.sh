#!/bin/bash
#Follow on from Rob's analysis and test different run settings for shapeit haplotyping.
set -eu
datadir=/home/sobczm/popgen/snp/snp_chip/haplotypes/comparison/rgxha_tweaking

#based on 130_compare_map2imputed.sh scripts
SCRIPT_PATH=/home/sobczm/bin/test_haplotyping
export PATH=${SCRIPT_PATH}:${PATH}
export PATH=/home/vicker/git_repos/crosslink/bin:${PATH}
source ${SCRIPT_PATH}/haplotyping_all_samples_funcs.sh
export NODES="-N2,5,6,7,8,9,10"
export MAXJOBS="-L30"
script=/home/sobczm/bin/popgen/snp
name=rgxha_test_ms

cd $datadir
mapfile=/home/vicker/fragaria_x_ananassa/octoseq/maps/consensus/vesca2consensus_map_integerposns_2017-07-09.csv

#list of markers present in the mapfile
cat ${mapfile} | cut -d, -f1,2 > marker_list
rm -f filtered_map.csv
mkdir popnordergrps && cd popnordergrps
original=/home/vicker/octoploid_mapping/test_haplotyping/rgxha
cp -r $original/popnordergrps/*loc $original/popnordergrps/*loc2 $original/popnordergrps/*map2 $original/popnordergrps/*tmp ./
cd ..
cp $original/parental_genotypes.tsv ./

for locfile in popnordergrps/*.loc
do
    lg=$(basename ${locfile} .loc)
    tmpfile=${locfile/loc/tmp}
    loc2file=${locfile/loc/loc2}
    map2file=${locfile/loc/map2}

    #map markers in current lg
    cat marker_list | grep ${lg} | cut -d, -f 1 > marker_list_${lg}

    #remove any loc file markers not in the mapfile and correct LG
    cat ${locfile} | fgrep -f marker_list_${lg} > ${tmpfile}

    #remove markers not in locfile
    cat ${tmpfile} | cut -d' ' -f1 > locmarker_list_${lg}
    cat ${mapfile} | grep ${lg} | fgrep -f locmarker_list_${lg} > ${map2file}
    cat ${map2file} >> filtered_map.csv

    #order locfile by mapfile
    orderlocbymap.py ${tmpfile} ${map2file} > ${loc2file}

 #this version modified to not require progeny names
    #also using bp position from consensus map not rgxha cm map
    #seems to list parental haplotypes in reverse order but progeny not affected by this
    mapbased_haplotypes2.py \
        ${loc2file} \
        ${map2file} \
        parental_genotypes.tsv \
        > ${locfile/loc/hap}

    #alternative implementation
    mapbased_haplotypes2_alt.py \
        ${loc2file} \
        ${map2file} \
        parental_genotypes.tsv \
        > ${locfile/loc/alt}

    conv_cl2affy3.py ${loc2file} parental_genotypes.tsv \
        | gzip \
        > ${locfile/loc/csv}.gz

    cat ${map2file} > ${locfile/\.loc/}_map.csv
done

zcat popnordergrps/1A.csv.gz | head -n 1 > ${name}.csv 
#Need to edit this to add the parents at the end of the line
sed -i 's/$/,RG(168),HA(169)/' ${name}.csv 


#Convert from .alt to .csv and overwrite the .csv files produced by conv_cl2affy3.py 
for altfile in popnordergrps/*.alt
do
python $scripts/convert_to_affy.py $altfile | gzip >${altfile/alt/csv}.gz
done 

rm -f ${name}.csv.gz
zcat popnordergrps/*.csv.gz >> ${name}.csv
gzip ${name}.csv

#reformat ready for shapeit
convert2shapeit filtered_map.csv ${name}.csv.gz ${name}_shapeit

#Default run (all default apart from eff pop size)
export SHAPEIT_EFF_POPN_SIZE=1000
export SHAPEIT_SEED=${RANDOM}${RANDOM}
export SHAPEIT_PRUNE=8
export SHAPEIT_BURN=7
export SHAPEIT_MAIN=20
export SHAPEIT_STATES=100
export SHAPEIT_WINDOW=2
impute_haplotypes ${name}_shapeit

##Start changing values
#Prune
#6
cp -r ${name}_shapeit ${name}_shapeit_prune_6
export SHAPEIT_PRUNE=6
impute_haplotypes ${name}_shapeit_prune_6
#10
cp -r ${name}_shapeit ${name}_shapeit_prune_10
export SHAPEIT_PRUNE=10
impute_haplotypes ${name}_shapeit_prune_10
#12
cp -r ${name}_shapeit ${name}_shapeit_prune_12
export SHAPEIT_PRUNE=12
impute_haplotypes ${name}_shapeit_prune_12

#Burn
export SHAPEIT_PRUNE=8
#5
cp -r ${name}_shapeit ${name}_shapeit_burn_5
export SHAPEIT_BURN=5
impute_haplotypes ${name}_shapeit_burn_5
#10
cp -r ${name}_shapeit ${name}_shapeit_burn_10
export SHAPEIT_BURN=10
impute_haplotypes ${name}_shapeit_burn_10
#15
cp -r ${name}_shapeit ${name}_shapeit_burn_15
export SHAPEIT_BURN=15
impute_haplotypes ${name}_shapeit_burn_15


#Main
export SHAPEIT_BURN=7
#10
cp -r ${name}_shapeit ${name}_shapeit_main_10
export SHAPEIT_MAIN=10
impute_haplotypes ${name}_shapeit_main_10
#50
cp -r ${name}_shapeit ${name}_shapeit_main_50
export SHAPEIT_MAIN=50
impute_haplotypes ${name}_shapeit_main_50
#75
cp -r ${name}_shapeit ${name}_shapeit_main_75
export SHAPEIT_MAIN=75
impute_haplotypes ${name}_shapeit_main_75

#Conditioning
export SHAPEIT_MAIN=20
#50
cp -r ${name}_shapeit ${name}_shapeit_states_50
export SHAPEIT_STATES=50
impute_haplotypes ${name}_shapeit_states_50
#200
cp -r ${name}_shapeit ${name}_shapeit_states_200
export SHAPEIT_STATES=200
impute_haplotypes ${name}_shapeit_states_200
#300
cp -r ${name}_shapeit ${name}_shapeit_states_300
export SHAPEIT_STATES=300
impute_haplotypes ${name}_shapeit_states_300

#Window
export SHAPEIT_STATES=100
#0.5
cp -r ${name}_shapeit ${name}_shapeit_window_05
export SHAPEIT_WINDOW=0.5
impute_haplotypes ${name}_shapeit_window_05
#3
cp -r ${name}_shapeit ${name}_shapeit_window_3
export SHAPEIT_WINDOW=3
impute_haplotypes ${name}_shapeit_window_3
#4
cp -r ${name}_shapeit ${name}_shapeit_window_4
export SHAPEIT_WINDOW=4
impute_haplotypes ${name}_shapeit_window_4

#For comparison purposes with map-based dataset, need to print out everything
#per line but the last 4 haplotypes (i.e. the parents)
for hap in *_shapeit*/*.haps
do
cat $hap | rev | cut -c 10- | rev >${hap/haps/haps2}
done

for locfile in popnordergrps/*.loc2
do
    lg=$(basename ${locfile} .loc2)
    conv_cl2AB2.py     ${locfile}  ${locfile/loc2/csv.gz}  >  ${name}_${lg}_orig_AB    #original
done

#Compare haplotypes
for hapfile in *_shapeit*/*.haps2
do
    lg=$(basename ${hapfile} .phased.haps)
    conv_phased2AB.py  ${hapfile}  > ${name}_${lg}_imputed_AB         
    compare_haplotypes.py ${name}_${lg}_imputed_AB ${name}_${lg}_orig_AB  > ${name}_${lg}_out
done