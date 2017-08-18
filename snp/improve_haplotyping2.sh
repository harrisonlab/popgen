#Follow on from Rob's analysis and test different run settings for shapeit haplotyping.
set -eu
datadir=/home/sobczm/popgen/snp/snp_chip/haplotypes/comparison/rgxha_tweaking_original

#based on 130_compare_map2imputed.sh scripts
SCRIPT_PATH=/home/sobczm/bin/test_haplotyping
export PATH=${SCRIPT_PATH}:${PATH}
export PATH=/home/vicker/git_repos/crosslink/bin:${PATH}
source ${SCRIPT_PATH}/haplotyping_all_samples_funcs.sh
export NODES="-N2,5,6,7,8,9,10"
export MAXJOBS="-L30"
scripts=/home/sobczm/bin/popgen/snp
name=rgxha_test_ms

#See what happens to phasing of RGxHA when map-based haplotypes from EMxFE or FLxCH fed as reference panel.

cd $datadir

original=/home/vicker/octoploid_mapping/consensus_map4/popn_EMxFE/map/popnordergrps
cp -r $original/* $datadir/popnordergrps_emxfe

for locfile in popnordergrps_emxfe/*.loc
do
    lg=$(basename ${locfile} .loc)
    tmpfile=${locfile/loc/tmp}
    loc2file=${locfile/loc/loc2}
    map2file=${locfile/loc/map2}

    #map markers in current lg
    cat marker_list | grep ${lg} | cut -d, -f 1 > popnordergrps_emxfe/marker_list_${lg}

    #remove any loc file markers not in the mapfile and correct LG
    cat ${locfile} | fgrep -f popnordergrps_emxfe/marker_list_${lg} > ${tmpfile}

    #remove markers not in locfile
    cat ${tmpfile} | cut -d' ' -f1 > popnordergrps_emxfe/locmarker_list_${lg}
    cat ${mapfile} | grep ${lg} | fgrep -f popnordergrps_emxfe/locmarker_list_${lg} > ${map2file}
    cat ${map2file} >> popnordergrps_emxfe/filtered_map.csv

    #order locfile by mapfile
    orderlocbymap.py ${tmpfile} ${map2file} > ${loc2file}

    #this version modified to not require progeny names
    #also using bp position from consensus map not rgxha cm map
    #seems to list parental haplotypes in reverse order but progeny not affected by this

    conv_cl2affy3.py ${loc2file} parental_genotypes.tsv \
        | gzip \
        > ${locfile/loc/csv}.gz

    cat ${map2file} > ${locfile/\.loc/}_map.csv
done

for locfile in popnordergrps_emxfe/*.loc2
do
    lg=$(basename ${locfile} .loc2)
    conv_cl2AB2.py     ${locfile}  ${locfile/loc2/csv.gz}  >  ${name}_${lg}_orig_AB_emxfe    #original
done

original=/home/vicker/octoploid_mapping/consensus_map4/popn_FLxCH/map/popnordergrps
cp -r $original/* $datadir/popnordergrps_flxch

for locfile in popnordergrps_flxch/*.loc
do
    lg=$(basename ${locfile} .loc)
    tmpfile=${locfile/loc/tmp}
    loc2file=${locfile/loc/loc2}
    map2file=${locfile/loc/map2}

    #map markers in current lg
    cat marker_list | grep ${lg} | cut -d, -f 1 > popnordergrps_flxch/marker_list_${lg}

    #remove any loc file markers not in the mapfile and correct LG
    cat ${locfile} | fgrep -f popnordergrps_flxch/marker_list_${lg} > ${tmpfile}

    #remove markers not in locfile
    cat ${tmpfile} | cut -d' ' -f1 > popnordergrps_flxch/locmarker_list_${lg}
    cat ${mapfile} | grep ${lg} | fgrep -f popnordergrps_flxch/locmarker_list_${lg} > ${map2file}
    cat ${map2file} >> popnordergrps_flxch/filtered_map.csv

    #order locfile by mapfile
    orderlocbymap.py ${tmpfile} ${map2file} > ${loc2file}

    conv_cl2affy3.py ${loc2file} parental_genotypes.tsv \
        | gzip \
        > ${locfile/loc/csv}.gz

    cat ${map2file} > ${locfile/\.loc/}_map.csv
done

for locfile in popnordergrps_flxch/*.loc2
do
    lg=$(basename ${locfile} .loc2)
    conv_cl2AB2.py     ${locfile}  ${locfile/loc2/csv.gz}  >  ${name}_${lg}_orig_AB_flxch    #original
done

cp -r $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_emxfe
list_inds=/home/vicker/octoploid_mapping/consensus_map4/popn_EMxFE/map/conf/popn_order
for number in 1 2 3 4 5 6 7
do
for subgenome in A B C D 
do
python $scripts/convert_to_impute2_general.py $datadir/rgxha_test_ms_${number}${subgenome}_orig_AB_emxfe $datadir/rgxha_test_ms_shapeit_emxfe/${number}${subgenome}.pmap $list_inds emxfe
done
done

cp -r $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_flxch
list_inds2=/home/vicker/octoploid_mapping/consensus_map4/popn_FLxCH/map/conf/popn_order
for number in 1 2 3 4 5 6 7
do
for subgenome in A B C D 
do
python $scripts/convert_to_impute2_general.py $datadir/rgxha_test_ms_${number}${subgenome}_orig_AB_flxch $datadir/rgxha_test_ms_shapeit_flxch/${number}${subgenome}.pmap $list_inds2 flxch
done
done

export SHAPEIT_EFF_POPN_SIZE=1000
export SHAPEIT_SEED=${RANDOM}${RANDOM}
export SHAPEIT_PRUNE=8
export SHAPEIT_BURN=7
export SHAPEIT_MAIN=20
export SHAPEIT_STATES=100
export SHAPEIT_WINDOW=2
#Phase using the reference panels - default settings
impute_haplotypes2 ${name}_shapeit_emxfe
impute_haplotypes2 ${name}_shapeit_flxch

#The same but with the best settings (combined)
cp -r ${name}_shapeit_emxfe ${name}_shapeit_emxfe_combined
cp -r ${name}_shapeit_flxch ${name}_shapeit_flxch_combined

export SHAPEIT_BURN=30
export SHAPEIT_STATES=500
export SHAPEIT_WINDOW=1.5

impute_haplotypes2 ${name}_shapeit_emxfe_combined
impute_haplotypes2 ${name}_shapeit_flxch_combined

for hapfile in *_shapeit_emxfe*/*.haps *_shapeit_flxch*/*.haps
do
    lg=$(basename ${hapfile} .phased.haps)
    conv_phased2AB.py  ${hapfile}  > ${hapfile/.phased.haps/_imputed_AB}   
    compare_haplotypes.py ${hapfile/.phased.haps/_imputed_AB} ${name}_${lg}_orig_AB  > ${hapfile/.phased.haps/_out}
done

Rscript --vanilla $scripts/plots_phasing.R EMxFE_reference $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_emxfe
Rscript --vanilla $scripts/plots_phasing.R EMxFE_reference_combined $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_emxfe_combined
Rscript --vanilla $scripts/plots_phasing.R FLxCH_reference $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_flxch
Rscript --vanilla $scripts/plots_phasing.R FLxCH_reference_combined $datadir/rgxha_test_ms_shapeit $datadir/rgxha_test_ms_shapeit_flxch_combined

for my_pdf in *_shapeit_emxfe*/*.pdf *_shapeit_flxch*/*.pdf 
do
cp $my_pdf $datadir/figures2
done

for my_pdf in $datadir/figures2/*.pdf
do
convert -verbose -density 500 "${my_pdf}" "${my_pdf%.*}.png"
done
