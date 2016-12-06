#!/bin/bash
input=/home/sobczm/popgen/codon/codonW/input
scripts=/home/sobczm/bin/popgen/codon
effectors=/home/sobczm/popgen/input/effectors

#Further dissection of Fus2 codon bias analysis results into lineage-specific and core region.

#Prepare independent annotation tables for lineage-specific and core region genes
cd $input/inhouse/notrans/Fus2_final_genes_combined.cdna_pass_one_notrans_only_an
python $scripts/separate_c_ls.py $effectors/CAZY/all/Fus2_canu_new_CAZY_gene_table.txt
python $scripts/separate_c_ls.py $effectors/CAZY/secreted/Fus2_canu_new_CAZY_secreted_gene_table.txt
python $scripts/separate_c_ls.py $effectors/effectorP/secreted/F.oxysporum_fsp_cepae_Fus2_canu_new_EffectorP_secreted_gene_table.txt
python $scripts/separate_c_ls.py $effectors/mimp/secreted/Fus2_canu_new_genes_in_2kb_mimp_secreted_gene_table.txt
python $scripts/separate_c_ls.py $effectors/mimp/all/Fus2_canu_new_genes_in_2kb_mimp_gene_table.txt
python $scripts/separate_c_ls.py Fus2_final_genes_appended_gene_table.txt

python $scripts/match_ids.py $effectors/CAZY/all/Fus2_canu_new_CAZY_gene_table_c.txt \
Fus2_final_genes_combined.cdna_pass_one_notrans.out
python $scripts/match_ids.py $effectors/CAZY/all/Fus2_canu_new_CAZY_gene_table_ls.txt \
Fus2_final_genes_combined.cdna_pass_one_notrans.out

python $scripts/match_ids.py $effectors/CAZY/secreted/Fus2_canu_new_CAZY_secreted_gene_table_c.txt \
Fus2_final_genes_combined.cdna_pass_one_notrans.out
python $scripts/match_ids.py $effectors/CAZY/secreted/Fus2_canu_new_CAZY_secreted_gene_table_ls.txt \
Fus2_final_genes_combined.cdna_pass_one_notrans.out

python $scripts/match_ids.py $effectors/effectorP/secreted/F.oxysporum_fsp_cepae_Fus2_canu_new_EffectorP_secreted_gene_table_c.txt \
Fus2_final_genes_combined.cdna_pass_one_notrans.out
python $scripts/match_ids.py $effectors/effectorP/secreted/F.oxysporum_fsp_cepae_Fus2_canu_new_EffectorP_secreted_gene_table_ls.txt \
Fus2_final_genes_combined.cdna_pass_one_notrans.out

python $scripts/match_ids.py $effectors/mimp/secreted/Fus2_canu_new_genes_in_2kb_mimp_secreted_gene_table_c.txt \
Fus2_final_genes_combined.cdna_pass_one_notrans.out
python $scripts/match_ids.py $effectors/mimp/secreted/Fus2_canu_new_genes_in_2kb_mimp_secreted_gene_table_ls.txt \
Fus2_final_genes_combined.cdna_pass_one_notrans.out

python $scripts/match_ids.py $effectors/mimp/all/Fus2_canu_new_genes_in_2kb_mimp_gene_table_c.txt \
Fus2_final_genes_combined.cdna_pass_one_notrans.out
python $scripts/match_ids.py $effectors/mimp/all/Fus2_canu_new_genes_in_2kb_mimp_gene_table_ls.txt \
Fus2_final_genes_combined.cdna_pass_one_notrans.out

python $scripts/match_ids.py Fus2_final_genes_appended_gene_table_c.txt \
Fus2_final_genes_combined.cdna_pass_one_notrans.out
python $scripts/match_ids.py Fus2_final_genes_appended_gene_table_ls.txt \
Fus2_final_genes_combined.cdna_pass_one_notrans.out

#And now for all the regions CB3_final_genes_combined
python $scripts/match_ids.py $effectors/CAZY/all/Fus2_canu_new_CAZY_gene_table.txt \
Fus2_final_genes_combined.cdna_pass_one_notrans.out

python $scripts/match_ids.py $effectors/CAZY/secreted/Fus2_canu_new_CAZY_secreted_gene_table.txt \
Fus2_final_genes_combined.cdna_pass_one_notrans.out

python $scripts/match_ids.py $effectors/effectorP/secreted/F.oxysporum_fsp_cepae_Fus2_canu_new_EffectorP_secreted_gene_table.txt \
Fus2_final_genes_combined.cdna_pass_one_notrans.out

python $scripts/match_ids.py $effectors/mimp/secreted/Fus2_canu_new_genes_in_2kb_mimp_secreted_gene_table.txt \
Fus2_final_genes_combined.cdna_pass_one_notrans.out

python $scripts/match_ids.py $effectors/mimp/all/Fus2_canu_new_genes_in_2kb_mimp_gene_table.txt \
Fus2_final_genes_combined.cdna_pass_one_notrans.out

#And append gene expression rank to each gene in the results file (add add. column):
python $scripts/match_rank.py gene_expression_rank_fus2.txt Fus2_final_genes_combined.cdna_pass_one_notrans.out

#Execute separately for lineage-specific (_ls) and core (_c) regions as well as all regions.
cp $scripts/codonw_analysis.R ./
mv codonw_analysis.R codonw_analysis_c.R
Rscript --vanilla codonw_analysis_c.R

cp $scripts/codonw_analysis.R ./
mv codonw_analysis.R codonw_analysis_ls.R
Rscript --vanilla codonw_analysis_ls.R

cp $scripts/codonw_analysis.R ./
Rscript --vanilla codonw_analysis.R
