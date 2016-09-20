#!/bin/bash
input=/home/sobczm/popgen/codon/codonW/input
scripts=/home/sobczm/bin/popgen/codon
effectors=/home/sobczm/popgen/input/effectors

#Get a selection of codonw outputs containing different
#effector subsets
cd $input/inhouse/only_annotated
cd 125_final_genes_combined.cdna_pass_one_notrans_only_an

################## Only genomes with annotation tables
#125
python $scripts/match_ids.py $effectors/CAZY/all/125_CAZY_gene_table.txt \
125_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/CAZY/secreted/125_CAZY_secreted_gene_table.txt \
125_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/effectorP/secreted/F.oxysporum_fsp_cepae_125_EffectorP_secreted_gene_table.txt \
125_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/mimp/secreted/125_genes_in_2kb_mimp_secreted_gene_table.txt \
125_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/mimp/all/125_genes_in_2kb_mimp_gene_table.txt \
125_final_genes_combined.cdna_pass_one_notrans_only_an.out

cp $scripts/codonw_analysis.R ./
Rscript --vanilla codonw_analysis.R

cd ../A13_final_genes_combined.cdna_pass_one_notrans_only_an
#A13
python $scripts/match_ids.py $effectors/CAZY/all/A13_CAZY_gene_table.txt \
A13_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/CAZY/secreted/A13_CAZY_secreted_gene_table.txt \
A13_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/effectorP/secreted/F.oxysporum_fsp_cepae_A13_EffectorP_secreted_gene_table.txt \
A13_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/mimp/secreted/A13_genes_in_2kb_mimp_secreted_gene_table.txt \
A13_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/mimp/all/A13_genes_in_2kb_mimp_gene_table.txt \
A13_final_genes_combined.cdna_pass_one_notrans_only_an.out

cp $scripts/codonw_analysis.R ./
sed -i -e 's/125_/A13_/g' codonw_analysis.R
Rscript --vanilla codonw_analysis.R

#A23
cd ../A23_final_genes_combined.cdna_pass_one_notrans_only_an
python $scripts/match_ids.py $effectors/CAZY/all/A23_CAZY_gene_table.txt \
A23_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/CAZY/secreted/A23_CAZY_secreted_gene_table.txt \
A23_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/effectorP/secreted/F.oxysporum_fsp_cepae_A23_EffectorP_secreted_gene_table.txt \
A23_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/mimp/secreted/A23_genes_in_2kb_mimp_secreted_gene_table.txt \
A23_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/mimp/all/A23_genes_in_2kb_mimp_gene_table.txt \
A23_final_genes_combined.cdna_pass_one_notrans_only_an.out

cp $scripts/codonw_analysis.R ./
sed -i -e 's/125_/A23_/g' codonw_analysis.R
Rscript --vanilla codonw_analysis.R


#A28
cd ../A28_final_genes_combined.cdna_pass_one_notrans_only_an
python $scripts/match_ids.py $effectors/CAZY/all/A28_CAZY_gene_table.txt \
A28_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/CAZY/secreted/A28_CAZY_secreted_gene_table.txt \
A28_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/effectorP/secreted/F.oxysporum_fsp_cepae_A28_EffectorP_secreted_gene_table.txt \
A28_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/mimp/secreted/A28_genes_in_2kb_mimp_secreted_gene_table.txt \
A28_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/mimp/all/A28_genes_in_2kb_mimp_gene_table.txt \
A28_final_genes_combined.cdna_pass_one_notrans_only_an.out

cp $scripts/codonw_analysis.R ./
sed -i -e 's/125_/A28_/g' codonw_analysis.R
Rscript --vanilla codonw_analysis.R



#CB3
cd ../CB3_final_genes_combined.cdna_pass_one_notrans_only_an
python $scripts/match_ids.py $effectors/CAZY/all/CB3_CAZY_gene_table.txt \
CB3_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/CAZY/secreted/CB3_CAZY_secreted_gene_table.txt \
CB3_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/effectorP/secreted/F.oxysporum_fsp_cepae_CB3_EffectorP_secreted_gene_table.txt \
CB3_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/mimp/secreted/CB3_genes_in_2kb_mimp_secreted_gene_table.txt \
CB3_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/mimp/all/CB3_genes_in_2kb_mimp_gene_table.txt \
CB3_final_genes_combined.cdna_pass_one_notrans_only_an.out

cp $scripts/codonw_analysis.R ./
sed -i -e 's/125_/CB3_/g' codonw_analysis.R
Rscript --vanilla codonw_analysis.R

#Fus2_canu_new
cd ../Fus2_final_genes_combined.cdna_pass_one_notrans_only_an
python $scripts/match_ids.py $effectors/CAZY/all/Fus2_canu_new_CAZY_gene_table.txt \
Fus2_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/CAZY/secreted/Fus2_canu_new_CAZY_secreted_gene_table.txt \
Fus2_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/effectorP/secreted/F.oxysporum_fsp_cepae_Fus2_canu_new_EffectorP_secreted_gene_table.txt \
Fus2_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/mimp/secreted/Fus2_canu_new_genes_in_2kb_mimp_secreted_gene_table.txt \
Fus2_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/mimp/all/Fus2_canu_new_genes_in_2kb_mimp_gene_table.txt \
Fus2_final_genes_combined.cdna_pass_one_notrans_only_an.out

cp $scripts/codonw_analysis.R ./
sed -i -e 's/125_/Fus2_canu_new_/g' codonw_analysis.R
Rscript --vanilla codonw_analysis.R

#PG
cd ../PG_final_genes_combined.cdna_pass_one_notrans_only_an
python $scripts/match_ids.py $effectors/CAZY/all/PG_CAZY_gene_table.txt \
PG_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/CAZY/secreted/PG_CAZY_secreted_gene_table.txt \
PG_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/effectorP/secreted/F.oxysporum_fsp_cepae_PG_EffectorP_secreted_gene_table.txt \
PG_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/mimp/secreted/PG_genes_in_2kb_mimp_secreted_gene_table.txt \
PG_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/mimp/all/PG_genes_in_2kb_mimp_gene_table.txt \
PG_final_genes_combined.cdna_pass_one_notrans_only_an.out

cp $scripts/codonw_analysis.R ./
sed -i -e 's/125_/PG_/g' codonw_analysis.R
Rscript --vanilla codonw_analysis.R

################## Rest of the in-houe genomes
cd $input/inhouse/all

#55
cd ../55_final_genes_combined.cdna_pass_one
python $scripts/match_ids.py $effectors/CAZY/all/55_CAZY_gene_table.txt \
55_final_genes_combined.cdna_pass_one.out

python $scripts/match_ids.py $effectors/CAZY/secreted/55_CAZY_secreted_gene_table.txt \
55_final_genes_combined.cdna_pass_one.out

python $scripts/match_ids.py $effectors/effectorP/secreted/F.oxysporum_fsp_cepae_55_EffectorP_secreted_gene_table.txt \
55_final_genes_combined.cdna_pass_one.out

python $scripts/match_ids.py $effectors/mimp/secreted/55_genes_in_2kb_mimp_secreted_gene_table.txt \
55_final_genes_combined.cdna_pass_one.out

python $scripts/match_ids.py $effectors/mimp/all/55_genes_in_2kb_mimp_gene_table.txt \
55_final_genes_combined.cdna_pass_one.out

cp $scripts/codonw_analysis.R ./
sed -i -e 's/125_/55_/g' codonw_analysis.R
sed -i -e 's/pass_one_notrans_only_an/pass_one/g' codonw_analysis.R
Rscript --vanilla codonw_analysis.R

#D2
cd ../D2_final_genes_combined.cdna_pass_one
python $scripts/match_ids.py $effectors/CAZY/all/D2_CAZY_gene_table.txt \
D2_final_genes_combined.cdna_pass_one.out

python $scripts/match_ids.py $effectors/CAZY/secreted/D2_CAZY_secreted_gene_table.txt \
D2_final_genes_combined.cdna_pass_one.out

python $scripts/match_ids.py $effectors/effectorP/secreted/F.oxysporum_fsp_cepae_D2_EffectorP_secreted_gene_table.txt \
D2_final_genes_combined.cdna_pass_one.out

python $scripts/match_ids.py $effectors/mimp/secreted/D2_genes_in_2kb_mimp_secreted_gene_table.txt \
D2_final_genes_combined.cdna_pass_one.out

python $scripts/match_ids.py $effectors/mimp/all/D2_genes_in_2kb_mimp_gene_table.txt \
D2_final_genes_combined.cdna_pass_one.out

cp $scripts/codonw_analysis.R ./
sed -i -e 's/125_/D2_/g' codonw_analysis.R
sed -i -e 's/pass_one_notrans_only_an/pass_one/g' codonw_analysis.R
Rscript --vanilla codonw_analysis.R

#proliferatum
cd ../proliferatum_final_genes_combined.cdna_pass_one
python $scripts/match_ids.py $effectors/CAZY/all/A8_CAZY_gene_table.txt \
proliferatum_final_genes_combined.cdna_pass_one.out

python $scripts/match_ids.py $effectors/CAZY/secreted/A8_CAZY_secreted_gene_table.txt \
proliferatum_final_genes_combined.cdna_pass_one.out

python $scripts/match_ids.py $effectors/effectorP/secreted/F.proliferatum_A8_EffectorP_secreted_gene_table.txt \
proliferatum_final_genes_combined.cdna_pass_one.out

python $scripts/match_ids.py $effectors/mimp/secreted/A8_genes_in_2kb_mimp_secreted_gene_table.txt \
proliferatum_final_genes_combined.cdna_pass_one.out

python $scripts/match_ids.py $effectors/mimp/all/A8_genes_in_2kb_mimp_gene_table.txt \
proliferatum_final_genes_combined.cdna_pass_one.out

cp $scripts/codonw_analysis.R ./
sed -i -e 's/125_/A8_/g' codonw_analysis.R
sed -i -e 's/pass_one_notrans_only_an/pass_one/g' codonw_analysis.R
Rscript --vanilla codonw_analysis.R
