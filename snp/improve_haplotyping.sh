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
scripts=/home/sobczm/bin/popgen/snp
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
cat $hap | rev | cut -c 8- | rev >${hap/haps/haps2}
done

for locfile in popnordergrps/*.loc2
do
    lg=$(basename ${locfile} .loc2)
    conv_cl2AB2.py     ${locfile}  ${locfile/loc2/csv.gz}  >  ${name}_${lg}_orig_AB    #original
done

#Compare haplotypes
for hapfile in *_shapeit*/*.haps2
do
    lg=$(basename ${hapfile} .phased.haps2)
    conv_phased2AB.py  ${hapfile}  > ${hapfile/.phased.haps2/_imputed_AB}   
    compare_haplotypes.py ${hapfile/.phased.haps2/_imputed_AB} ${name}_${lg}_orig_AB  > ${hapfile/.phased.haps2/_out}
done

#Print comparison figures for shapeit runs using original dataset
#(without parents, just progeny)
datadir=/home/sobczm/popgen/snp/snp_chip/haplotypes/comparison/rgxha_tweaking_original
cd $datadir
Rscript --vanilla $scripts/plots_phasing.R burn_10 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_burn_10

Rscript --vanilla $scripts/plots_phasing.R burn_15 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_burn_15

Rscript --vanilla $scripts/plots_phasing.R burn_5 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_burn_5

Rscript --vanilla $scripts/plots_phasing.R main_10 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_main_10

Rscript --vanilla $scripts/plots_phasing.R main_50 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_main_50

Rscript --vanilla $scripts/plots_phasing.R main_75 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_main_75

Rscript --vanilla $scripts/plots_phasing.R prune_10 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_prune_10

Rscript --vanilla $scripts/plots_phasing.R prune_12 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_prune_12

Rscript --vanilla $scripts/plots_phasing.R prune_6 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_prune_6

Rscript --vanilla $scripts/plots_phasing.R states_200 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_states_200

Rscript --vanilla $scripts/plots_phasing.R states_300 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_states_300

Rscript --vanilla $scripts/plots_phasing.R states_50 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_states_50

Rscript --vanilla $scripts/plots_phasing.R window_05 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_window_05

Rscript --vanilla $scripts/plots_phasing.R window_05 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_window_05

Rscript --vanilla $scripts/plots_phasing.R window_3 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_window_3

Rscript --vanilla $scripts/plots_phasing.R window_4 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_window_4

#Copy all the PDFs to a singledir
for my_pdf in *_shapeit*/*.pdf 
do
cp $my_pdf $datadir/figures
done

#Compare phasing with parents added to that without
Rscript --vanilla $scripts/plots_phasing.R parents_added $datadir/rgxha_test_ms_shapeit /home/sobczm/popgen/snp/snp_chip/haplotypes/comparison/rgxha_tweaking/rgxha_test_ms_shapeit

cp /home/sobczm/popgen/snp/snp_chip/haplotypes/comparison/rgxha_tweaking/rgxha_test_ms_shapeit/*.pdf $datadir/figures

#Convert to png 
for my_pdf in $datadir/figures/*.pdf
do
convert -verbose -density 500 "${my_pdf}" "${my_pdf%.*}.png"
done

#In the original dataset directory, fine-tuning additional options.
#burn 
#15
cp -r ${name}_shapeit ${name}_shapeit_burn_15
export SHAPEIT_BURN=15
impute_haplotypes ${name}_shapeit_burn_15
#20
cp -r ${name}_shapeit ${name}_shapeit_burn_20
export SHAPEIT_BURN=20
impute_haplotypes ${name}_shapeit_burn_20
#25
cp -r ${name}_shapeit ${name}_shapeit_burn_25
export SHAPEIT_BURN=25
impute_haplotypes ${name}_shapeit_burn_25
#30
cp -r ${name}_shapeit ${name}_shapeit_burn_30
export SHAPEIT_BURN=30
impute_haplotypes ${name}_shapeit_burn_30

#states 
export SHAPEIT_BURN=7
#300
cp -r ${name}_shapeit ${name}_shapeit_states_300
export SHAPEIT_STATES=300
impute_haplotypes ${name}_shapeit_states_300
#375
cp -r ${name}_shapeit ${name}_shapeit_states_375
export SHAPEIT_STATES=375
impute_haplotypes ${name}_shapeit_states_375
#450
cp -r ${name}_shapeit ${name}_shapeit_states_450
export SHAPEIT_STATES=450
impute_haplotypes ${name}_shapeit_states_450
#525
cp -r ${name}_shapeit ${name}_shapeit_states_525
export SHAPEIT_STATES=525
impute_haplotypes ${name}_shapeit_states_525

export SHAPEIT_STATES=100
#1.8
cp -r ${name}_shapeit ${name}_shapeit_window_1_8
export SHAPEIT_WINDOW=1.8
impute_haplotypes ${name}_shapeit_window_1_8
#1.5
cp -r ${name}_shapeit ${name}_shapeit_window_1_5
export SHAPEIT_WINDOW=1.5
impute_haplotypes ${name}_shapeit_window_1_5
#1.2
cp -r ${name}_shapeit ${name}_shapeit_window_1_2
export SHAPEIT_WINDOW=1.2
impute_haplotypes ${name}_shapeit_window_1_2
#0.9
cp -r ${name}_shapeit ${name}_shapeit_window_0_9
export SHAPEIT_WINDOW=0.9
impute_haplotypes ${name}_shapeit_window_0_9

#Plot 
Rscript --vanilla $scripts/plots_phasing.R burn_15 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_burn_15
Rscript --vanilla $scripts/plots_phasing.R burn_20 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_burn_20
Rscript --vanilla $scripts/plots_phasing.R burn_25 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_burn_25
Rscript --vanilla $scripts/plots_phasing.R burn_30 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_burn_30

Rscript --vanilla $scripts/plots_phasing.R states_300 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_states_300
Rscript --vanilla $scripts/plots_phasing.R states_375 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_states_375
Rscript --vanilla $scripts/plots_phasing.R states_450 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_states_450
Rscript --vanilla $scripts/plots_phasing.R states_525 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_states_525

Rscript --vanilla $scripts/plots_phasing.R window_1_8 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_window_1_8
Rscript --vanilla $scripts/plots_phasing.R window_1_5 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_window_1_5
Rscript --vanilla $scripts/plots_phasing.R window_1_2 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_window_1_2
Rscript --vanilla $scripts/plots_phasing.R window_0_9 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_window_0_9
#Test all the selected changed options combined, for state, window_size and burn.
export SHAPEIT_BURN=30
export SHAPEIT_STATES=500
export SHAPEIT_WINDOW=1.5
cp -r ${name}_shapeit ${name}_shapeit_combined_1
impute_haplotypes ${name}_shapeit_combined_1

cp -r ${name}_shapeit ${name}_shapeit_combined_1_rep1
impute_haplotypes ${name}_shapeit_combined_1_rep1

cp -r ${name}_shapeit ${name}_shapeit_combined_1_rep2
impute_haplotypes ${name}_shapeit_combined_1_rep2

#To convert map-based haplotypes to the IMPUTE2 format needed for reference panel option,
#Use the orig_A_B file to get haplotypes and the shapeit input files to get other info
Rscript --vanilla $scripts/plots_phasing.R combined_rep1 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_combined_1
Rscript --vanilla $scripts/plots_phasing.R combined_rep2 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_combined_1_rep1
Rscript --vanilla $scripts/plots_phasing.R combined_rep3 $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_combined_1_rep2

