#!/bin/bash
#Obtain lists of effector genes from Fusarium oxysporum genomes
fus=/home/groups/harrisonlab/project_files/fusarium/gene_pred
input=/home/sobczm/popgen/input

mkdir -p $input/effectors/CAZY/secreted
cd $input/effectors/CAZY/secreted
# Secreted CAZY genes
cp $fus/CAZY/F.proliferatum/A8/A8_CAZY_secreted.gff ./
cp $fus/CAZY/F.oxysporum_fsp_lycopersici/4287/4287_CAZY_secreted.gff ./
cp $fus/CAZY/F.oxysporum/fo47/fo47_CAZY_secreted.gff ./
cp $fus/CAZY/F.oxysporum_fsp_cepae/125/125_CAZY_secreted.gff ./
cp $fus/CAZY/F.oxysporum_fsp_cepae/D2/D2_CAZY_secreted.gff ./
cp $fus/CAZY/F.oxysporum_fsp_cepae/55/55_CAZY_secreted.gff ./
cp $fus/CAZY/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_CAZY_secreted.gff ./
cp $fus/CAZY/F.oxysporum_fsp_cepae/A1-2/A1-2_CAZY_secreted.gff ./
cp $fus/CAZY/F.oxysporum_fsp_cepae/HB6/HB6_CAZY_secreted.gff ./
cp $fus/CAZY/F.oxysporum_fsp_cepae/A13/A13_CAZY_secreted.gff ./
cp $fus/CAZY/F.oxysporum_fsp_cepae/PG/PG_CAZY_secreted.gff ./
cp $fus/CAZY/F.oxysporum_fsp_cepae/A23/A23_CAZY_secreted.gff ./
cp $fus/CAZY/F.oxysporum_fsp_cepae/A28/A28_CAZY_secreted.gff ./
cp $fus/CAZY/F.oxysporum_fsp_cepae/CB3/CB3_CAZY_secreted.gff ./

mkdir -p $input/effectors/CAZY/all
cd $input/effectors/CAZY/all
# All CAZY genes
cp $fus/CAZY/F.proliferatum/A8/A8_CAZY.gff ./
cp $fus/CAZY/F.oxysporum_fsp_lycopersici/4287/4287_CAZY.gff ./
cp $fus/CAZY/F.oxysporum/fo47/fo47_CAZY.gff ./
cp $fus/CAZY/F.oxysporum_fsp_cepae/125/125_CAZY.gff ./
cp $fus/CAZY/F.oxysporum_fsp_cepae/D2/D2_CAZY.gff ./
cp $fus/CAZY/F.oxysporum_fsp_cepae/55/55_CAZY.gff ./
cp $fus/CAZY/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_CAZY.gff ./
cp $fus/CAZY/F.oxysporum_fsp_cepae/A1-2/A1-2_CAZY.gff ./
cp $fus/CAZY/F.oxysporum_fsp_cepae/HB6/HB6_CAZY.gff ./
cp $fus/CAZY/F.oxysporum_fsp_cepae/A13/A13_CAZY.gff ./
cp $fus/CAZY/F.oxysporum_fsp_cepae/PG/PG_CAZY.gff ./
cp $fus/CAZY/F.oxysporum_fsp_cepae/A23/A23_CAZY.gff ./
cp $fus/CAZY/F.oxysporum_fsp_cepae/A28/A28_CAZY.gff ./
cp $fus/CAZY/F.oxysporum_fsp_cepae/CB3/CB3_CAZY.gff ./

mkdir -p $input/effectors/effectorP/secreted
cd $input/effectors/effectorP/secreted
# Secreted EffectorP proteins
fa=/home/groups/harrisonlab/project_files/fusarium/analysis
cp $fa/effectorP/F.proliferatum/A8/F.proliferatum_A8_EffectorP_secreted.gff ./
#cp $fa/effectorP/F.oxysporum_fsp_lycopersici/4287/4287_EffectorP_secreted.gff ./ #NA
#cp $fa/effectorP/F.oxysporum/fo47/fo47_EffectorP_secreted.gff ./ #NA
cp $fa/effectorP/F.oxysporum_fsp_cepae/125/F.oxysporum_fsp_cepae_125_EffectorP_secreted.gff ./
cp $fa/effectorP/F.oxysporum_fsp_cepae/D2/F.oxysporum_fsp_cepae_D2_EffectorP_secreted.gff ./
cp $fa/effectorP/F.oxysporum_fsp_cepae/55/F.oxysporum_fsp_cepae_55_EffectorP_secreted.gff ./
cp $fa/effectorP/F.oxysporum_fsp_cepae/Fus2_canu_new/F.oxysporum_fsp_cepae_Fus2_canu_new_EffectorP_secreted.gff ./
cp $fa/effectorP/F.oxysporum_fsp_cepae/A1-2/F.oxysporum_fsp_cepae_A1-2_EffectorP_secreted.gff ./
cp $fa/effectorP/F.oxysporum_fsp_cepae/HB6/F.oxysporum_fsp_cepae_HB6_EffectorP_secreted.gff ./
cp $fa/effectorP/F.oxysporum_fsp_cepae/A13/F.oxysporum_fsp_cepae_A13_EffectorP_secreted.gff ./
cp $fa/effectorP/F.oxysporum_fsp_cepae/PG/F.oxysporum_fsp_cepae_PG_EffectorP_secreted.gff ./
cp $fa/effectorP/F.oxysporum_fsp_cepae/A23/F.oxysporum_fsp_cepae_A23_EffectorP_secreted.gff ./
cp $fa/effectorP/F.oxysporum_fsp_cepae/A28/F.oxysporum_fsp_cepae_A28_EffectorP_secreted.gff ./
cp $fa/effectorP/F.oxysporum_fsp_cepae/CB3/F.oxysporum_fsp_cepae_CB3_EffectorP_secreted.gff ./

#Get rid of the prefix in the filename
for f in *.gff; do echo $f | sed -r 's/.*(F.oxysporum_fsp_cepae_)(.*)/\2/'; done

mkdir -p $input/effectors/mimp/all
cd $input/effectors/mimp/all
#All genes within 2kbp of MIMPs
cp $fa/mimps/F.proliferatum/A8/A8_genes_in_2kb_mimp.gff ./
cp $fa/mimps/F.oxysporum_fsp_lycopersici/4287/4287_genes_in_2kb_mimp.gff ./
cp $fa/mimps/F.oxysporum/fo47/fo47_genes_in_2kb_mimp.gff ./
cp $fa/mimps/F.oxysporum_fsp_cepae/125/125_genes_in_2kb_mimp.gff ./
cp $fa/mimps/F.oxysporum_fsp_cepae/D2/D2_genes_in_2kb_mimp.gff ./
cp $fa/mimps/F.oxysporum_fsp_cepae/55/55_genes_in_2kb_mimp.gff ./
cp $fa/mimps/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_genes_in_2kb_mimp.gff ./
cp $fa/mimps/F.oxysporum_fsp_cepae/A1-2/A1-2_genes_in_2kb_mimp.gff ./
cp $fa/mimps/F.oxysporum_fsp_cepae/HB6/HB6_genes_in_2kb_mimp.gff ./
cp $fa/mimps/F.oxysporum_fsp_cepae/A13/A13_genes_in_2kb_mimp.gff ./
cp $fa/mimps/F.oxysporum_fsp_cepae/PG/PG_genes_in_2kb_mimp.gff ./
cp $fa/mimps/F.oxysporum_fsp_cepae/A23/A23_genes_in_2kb_mimp.gff ./
cp $fa/mimps/F.oxysporum_fsp_cepae/A28/A28_genes_in_2kb_mimp.gff ./
cp $fa/mimps/F.oxysporum_fsp_cepae/CB3/CB3_genes_in_2kb_mimp.gff ./


mkdir -p $input/effectors/mimp/secreted
cd $input/effectors/mimp/secreted
#All secreted proteins within 2kbp of MIMPs
cp $fa/mimps/F.proliferatum/A8/A8_genes_in_2kb_mimp_secreted.gff ./
cp $fa/mimps/F.oxysporum_fsp_lycopersici/4287/4287_genes_in_2kb_mimp_secreted.gff ./
cp $fa/mimps/F.oxysporum/fo47/fo47_genes_in_2kb_mimp_secreted.gff ./
cp $fa/mimps/F.oxysporum_fsp_cepae/125/125_genes_in_2kb_mimp_secreted.gff ./
cp $fa/mimps/F.oxysporum_fsp_cepae/D2/D2_genes_in_2kb_mimp_secreted.gff ./
cp $fa/mimps/F.oxysporum_fsp_cepae/55/55_genes_in_2kb_mimp_secreted.gff ./
cp $fa/mimps/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_genes_in_2kb_mimp_secreted.gff ./
cp $fa/mimps/F.oxysporum_fsp_cepae/A1-2/A1-2_genes_in_2kb_mimp_secreted.gff ./
cp $fa/mimps/F.oxysporum_fsp_cepae/HB6/HB6_genes_in_2kb_mimp_secreted.gff ./
cp $fa/mimps/F.oxysporum_fsp_cepae/A13/A13_genes_in_2kb_mimp_secreted.gff ./
cp $fa/mimps/F.oxysporum_fsp_cepae/PG/PG_genes_in_2kb_mimp_secreted.gff ./
cp $fa/mimps/F.oxysporum_fsp_cepae/A23/A23_genes_in_2kb_mimp_secreted.gff ./
cp $fa/mimps/F.oxysporum_fsp_cepae/A28/A28_genes_in_2kb_mimp_secreted.gff ./
cp $fa/mimps/F.oxysporum_fsp_cepae/CB3/CB3_genes_in_2kb_mimp_secreted.gff ./


#Obtain a table of effector genes with cds_to_chromosome_coords.py
find /home/sobczm/popgen/input/effectors -type f -name '*.gff' -exec /home/sobczm/bin/popgen/codon/cds_to_chromosome_coords.py {} \;
#Can't do that as need to feed prefix, so have to run the script manually -.-
