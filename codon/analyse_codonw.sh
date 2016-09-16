#!/bin/bash
input=/home/sobczm/popgen/codon/codonW/input
scripts=/home/sobczm/bin/popgen/codon
effectors=/home/sobczm/popgen/input/effectors

#Get a selection of codonw output containing different
#effector subsets
cd $input/inhouse/only_annotated
cd 125_final_genes_combined.cdna_pass_one_notrans_only_an

python $scripts/match_ids.py $effectors/CAZY/all/125_CAZY_gene_table.txt \
125_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/CAZY/secreted/125_CAZY_secreted_gene_table.txt \
125_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/effectorP/secreted/F.oxysporum_fsp_cepae_125_EffectorP_secreted_gene_table.txt \
125_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/mimp/secreted/125_genes_in_2kb_mimp_secreted_gene_table.txt \
125_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/CAZY/all/55_CAZY_gene_table.txt \
55_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/CAZY/secreted/55_CAZY_secreted_gene_table.txt \
55_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/effectorP/secreted/F.oxysporum_fsp_cepae_55_EffectorP_secreted_gene_table.txt \
55_final_genes_combined.cdna_pass_one_notrans_only_an.out

python $scripts/match_ids.py $effectors/mimp/secreted/55_genes_in_2kb_mimp_secreted_gene_table.txt \
55_final_genes_combined.cdna_pass_one_notrans_only_an.out
