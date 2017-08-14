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

#Convert from 