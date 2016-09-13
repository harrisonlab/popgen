#!/bin/bash
input=/home/sobczm/popgen/input
scripts=/home/sobczm/bin/popgen/codon
#Create a directory for all codonW analyses (copy codonW binary folder)
cd /home/sobczm/popgen/codon/
cp -r /home/sobczm/bin/codonW /home/sobczm/popgen/codon
mkdir codonW/input && cd codonW/input
#Fetch CDS sequences for all external Fusarium oxysporum genomes available as input for analysis
#And clean them up before codonW run_codonw
$scripts/check_proper_cds.py $input/conglutians/Fusarium_oxysporum_f_sp_conglutinans_race_2_54008.FO_PHW808_V1.cds.all.fa
cp $input/conglutians/*pass.fasta ./ && rm $input/conglutians/*fail.fasta
$scripts/check_proper_cds.py $input/cubense/Fusarium_oxysporum_f_sp_cubense_race_1.Foc1_1.0.cds.all.fa
cp $input/cubense/*pass.fasta ./ && rm $input/cubense/*fail.fasta
$scripts/check_proper_cds.py $input/fo47/Fusarium_oxysporum_fo47.FO_Fo47_V1.cds.all.fa
cp $input/fo47/*pass.fasta ./ && rm $input/fo47/*fail.fasta
$scripts/check_proper_cds.py $input/fujikuroi/Fusarium_fujikuroi.EF1.cds.all.fa
cp $input/fujikuroi/*pass.fasta ./ && rm $input/fujikuroi/*fail.fasta
$scripts/check_proper_cds.py $input/lycopersici/Fusarium_oxysporum_f_sp_lycopersici_mn25.FO_MN25_V1.cds.all.fa
cp $input/lycopersici/*pass.fasta ./ && rm $input/lycopersici/*fail.fasta
$scripts/check_proper_cds.py $input/melonis/Fusarium_oxysporum_f_sp_melonis_26406.FO_melonis_V1.cds.all.fa
cp $input/melonis/*pass.fasta ./ && rm $input/melonis/*fail.fasta
$scripts/check_proper_cds.py $input/pisi/Fusarium_oxysporum_f_sp_pisi_hdv247.FO_HDV247_V1.cds.all.fa
cp $input/pisi/*pass.fasta ./ && rm $input/pisi/*fail.fasta
$scripts/check_proper_cds.py $input/radices-lycopersici/Fusarium_oxysporum_f_sp_radicis_lycopersici_26381.FO_CL57_V1.cds.all.fa
cp $input/radices-lycopersici/*pass.fasta ./ && rm $input/radices-lycopersici/*fail.fasta
$scripts/check_proper_cds.py $input/raphani/Fusarium_oxysporum_f_sp_raphani_54005.FO_PHW815_V1.cds.all.fa
cp $input/raphani/*pass.fasta ./ && rm $input/raphani/*fail.fasta
pw $input/vasinifectum/Fusarium_oxysporum_f_sp_vasinfectum_25433.FO_Cotton_V1.cds.all.fa
cp $input/vasinifectum/*pass.fasta ./ && rm $input/vasinifectum/*fail.fasta
#Does not work for this species, as contigs do not contain a stop codon
$scripts/check_proper_cds.py $input/verticillioides/Fusarium_verticillioides.ASM14955v1.cds.all.fa
cp $input/verticillioides/*pass.fasta ./ && rm $input/verticillioides/*fail.fasta

#Fetch additional Fusarium genomes
##F. graminearum  (assembly RR)
cd $input
mkdir graminearum
cd graminearum
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-32/fasta/fusarium_graminearum/cds/Fusarium_graminearum.RR.cds.all.fa.gz
gunzip -d Fusarium_graminearum.RR.cds.all.fa.gz
for fasta in *.fa; do sed -i -e 's/>/>graminearum_/' $fasta; done;
$scripts/check_proper_cds.py $input/graminearum/Fusarium_graminearum.RR.cds.all.fa
cp $input/graminearum/*pass.fasta ./ && rm $input/graminearum/*fail.fasta

##F. pseudograminearum (assembly GCA_000303195)
cd $input
mkdir pseudograminearum
cd pseudograminearum
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-32/fasta/fusarium_pseudograminearum/cds/Fusarium_pseudograminearum.GCA_000303195.1.cds.all.fa.gz
gunzip -d Fusarium_pseudograminearum.GCA_000303195.1.cds.all.fa.gz
for fasta in *.fa; do sed -i -e 's/>/>pseudograminearum_/' $fasta; done;
$scripts/check_proper_cds.py $input/pseudograminearum/Fusarium_pseudograminearum.GCA_000303195.1.cds.all.fa
cp $input/pseudograminearum/*pass.fasta ./ && rm $input/pseudograminearum/*fail.fasta
##F. langsethiae (assembly ASM129263v1)
cd $input
mkdir langsethiae
cd langsethiae
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-32/fasta/fungi_ascomycota3_collection/fusarium_langsethiae/cds/Fusarium_langsethiae.ASM129263v1.cds.all.fa.gz
gunzip -d Fusarium_langsethiae.ASM129263v1.cds.all.fa.gz
for fasta in *.fa; do sed -i -e 's/>/>langsethiae_/' $fasta; done;
$scripts/check_proper_cds.py $input/langsethiae/Fusarium_langsethiae.ASM129263v1.cds.all.fa
cp $input/langsethiae/*pass.fasta ./ && rm $input/langsethiae/*fail.fasta

#In-house Fusarium oxysporum genomes
$scripts/check_proper_cds.py $input/125/final/125_final_genes_combined.cdna.fasta
cp $input/125/final/*pass.fasta ./ && rm $input/125/final/*fail.fasta
$scripts/check_proper_cds.py $input/55/final/55_final_genes_combined.cdna.fasta
cp $input/55/final/*pass.fasta ./ && rm $input/55/final/*fail.fasta
$scripts/check_proper_cds.py $input/A1-2/final/A1-2_final_genes_combined.cdna.fasta
cp $input/A1-2/final/*pass.fasta ./ && rm $input/A1-2/final/*fail.fasta
$scripts/check_proper_cds.py $input/A13/final/A13_final_genes_combined.cdna.fasta
cp $input/A13/final/*pass.fasta ./ && rm $input/A13/final/*fail.fasta
$scripts/check_proper_cds.py $input/A23/final/A23_final_genes_combined.cdna.fasta
cp $input/A23/final/*pass.fasta ./ && rm $input/A23/final/*fail.fasta
$scripts/check_proper_cds.py $input/A28/final/A28_final_genes_combined.cdna.fasta
cp $input/A28/final/*pass.fasta ./ && rm $input/A28/final/*fail.fasta
$scripts/check_proper_cds.py $input/CB3/final/CB3_final_genes_combined.cdna.fasta
cp $input/CB3/final/*pass.fasta ./ && rm $input/CB3/final/*fail.fasta
$scripts/check_proper_cds.py $input/D2/final/D2_final_genes_combined.cdna.fasta
cp $input/D2/final/*pass.fasta ./ && rm $input/D2/final/*fail.fasta
$scripts/check_proper_cds.py $input/Fus2_canu_new/final/Fus2_final_genes_combined.cdna.fasta
cp $input/Fus2_canu_new/final/*pass.fasta ./ && rm $input/Fus2_canu_new/final/*fail.fasta
$scripts/check_proper_cds.py $input/HB6/final/HB6_final_genes_combined.cdna.fasta
cp $input/HB6/final/*pass.fasta ./ && rm $input/HB6/final/*fail.fasta
$scripts/check_proper_cds.py $input/PG/final/PG_final_genes_combined.cdna.fasta
cp $input/PG/final/*pass.fasta ./ && rm $input/PG/final/*fail.fasta
$scripts/check_proper_cds.py $input/proliferatum/final/proliferatum_final_genes_combined.cdna.fasta
cp $input/proliferatum/final/*pass.fasta ./ && rm $input/proliferatum/final/*fail.fasta
#Fusarium veneatum
cd $input
mkdir veneatum
cd veneatum
cp -r /home/groups/harrisonlab/project_files/fusarium_venenatum/gene_pred/final/F.venenatum/strain1/final ./
for filename in *; do mv "$filename" "veneatum_$filename"; done;
for fasta in *.fasta; do sed -i -e 's/>/>veneatum_/' $fasta; done;
$scripts/check_proper_cds.py  $input/veneatum/final/veneatum_final_genes_Braker.cds.fasta
cp $input/veneatum/final/*pass.fasta ./ && rm $input/veneatum/final/*fail.fasta
#In the future, may want to include F. oxysporum f. pisi and narcissi



#################
# CODONW RUNS
#################
#- 1) All genes
## IN-HOUSE GENOMES
#- 2) All genes minus alternative transcripts
#- 3) All genes minus alternative transcripts minus transposons
#- 4) All genes minus alternative transcripts minus transposons minus genes with no annotation

cw=/home/sobczm/popgen/codon/codonW
#1) All genes
cd $cw
for f in $cw/input/all/*.fasta;
do
cd $cw/input/all/
mkdir ${f%.*}
cp $f ${f%.*}/
cd ${f%.*}
filename=$(basename "$f")
# -cutot: tabulation of total codon usage
$cw/codonw $filename -nomenu -silent -cutot
mv ${filename%.*}.blk ${filename%.*}.cutot
# -cutab: tabulation of codon usage by gene
$cw/codonw $filename -nomenu -silent -cutab
mv ${filename%.*}.blk ${filename%.*}.cutab
# Correspondence analyscdis
# coa_num: percentage of genes taken from each tail of Nec distribution
$cw/codonw $filename -nomenu -silent -coa_cu -coa_num 5%
# Calculate all codon bias indices
# Use "cai.coa", "cbi.coa" and "fop.coa".  generated during correspondence
# Analysis to calculate the indices CAI, CBI and Fop
$cw/codonw $filename -all_indices -nomenu -silent -fop_file fop.coa -cai_file cai.coa -cbi_file cbi.coa
done

#2) All genes minus alternative transcripts
mkdir $cw/input/inhouse
cd $cw/input/all
cp 125_final_genes_combined.cdna_pass.fasta $cw/input/inhouse
cp 55_final_genes_combined.cdna_pass.fasta $cw/input/inhouse
cp A1-2_final_genes_combined.cdna_pass.fasta $cw/input/inhouse
cp A13_final_genes_combined.cdna_pass.fasta $cw/input/inhouse
cp A23_final_genes_combined.cdna_pass.fasta $cw/input/inhouse
cp A28_final_genes_combined.cdna_pass.fasta $cw/input/inhouse
cp CB3_final_genes_combined.cdna_pass.fasta $cw/input/inhouse
cp D2_final_genes_combined.cdna_pass.fasta $cw/input/inhouse
cp Fus2_final_genes_combined.cdna_pass.fasta $cw/input/inhouse
cp HB6_final_genes_combined.cdna_pass.fasta $cw/input/inhouse
cp PG_final_genes_combined.cdna_pass.fasta $cw/input/inhouse
cp proliferatum_final_genes_combined.cdna_pass.fasta $cw/input/inhouse
cp veneatum_final_genes_Braker.cds_pass.fasta $cw/input/inhouse
cd $cw/input/inhouse

for f in $cw/input/inhouse/*.fasta
do
$scripts/keep_one_gene.py $f
done

for f in $cw/input/inhouse/*one.fasta;
do
cd $cw/input/inhouse/
mkdir ${f%.*}
cp $f ${f%.*}/
cd ${f%.*}
filename=$(basename "$f")
# -cutot: tabulation of total codon usage
$cw/codonw $filename -nomenu -silent -cutot
mv ${filename%.*}.blk ${filename%.*}.cutot
# -cutab: tabulation of codon usage by gene
$cw/codonw $filename -nomenu -silent -cutab
mv ${filename%.*}.blk ${filename%.*}.cutab
# Correspondence analyscdis
# coa_num: percentage of genes taken from each tail of Nec distribution
$cw/codonw $filename -nomenu -silent -coa_cu -coa_num 5%
# Calculate all codon bias indices
# Use "cai.coa", "cbi.coa" and "fop.coa".  generated during correspondence
# Analysis to calculate the indices CAI, CBI and Fop
$cw/codonw $filename -all_indices -nomenu -silent -fop_file fop.coa -cai_file cai.coa -cbi_file cbi.coa
done

#- 3) All genes minus alternative transcripts minus transposons
