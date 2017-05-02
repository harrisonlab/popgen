#!/bin/bash
input=/home/sobczm/popgen/summary_stats
scripts=/home/sobczm/bin/popgen/summary_stats
#Obtain the VCF file with SNPs.
cd $input/noA13
#Find only SNPs specific to samples A23, Fus2, 55 and 125 and polarized against non-pathogens
python $scripts/vcf_find_difference_pop.py --vcf Fus2_canu_contigs_unmasked_noA13_filtered.vcf --out Fus2_canu_contigs_unmasked_pathogen.vcf --ply 1 --pop1 FOC125,,FOCFus2,,FOC55,,FOCA23 --pop2 FOCA1-2,,FOCA28,,FOCCB3,,FOCD2,,FOCHB6,,FOCPG --thr 0.95

#Find only SNPs specific to samples A23, Fus2, 55 and 125 and polarized against non-pathogens, allowing low frequency of that variant amongst the non-pathogens.
#However, the result identical
python $scripts/vcf_find_difference_pop.py --vcf Fus2_canu_contigs_unmasked_noA13_filtered.vcf --out Fus2_canu_contigs_unmasked_pathogen2.vcf --ply 1 --pop1 FOC125,,FOCFus2,,FOC55,,FOCA23 --pop2 FOCA1-2,,FOCA28,,FOCCB3,,FOCD2,,FOCHB6,,FOCPG --thr 0.25

#Check for overlap with different GFF files.
fus=/home/groups/harrisonlab/project_files/fusarium/gene_pred
cp $fus/CAZY/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_CAZY_secreted.gff ./
cp $fus/CAZY/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_CAZY.gff ./
fa=/home/groups/harrisonlab/project_files/fusarium/analysis
cp $fa/effectorP/F.oxysporum_fsp_cepae/Fus2_canu_new/F.oxysporum_fsp_cepae_Fus2_canu_new_EffectorP_secreted.gff ./
cp $fa/mimps/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_genes_in_2kb_mimp.gff ./
cp $fa/mimps/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_genes_in_2kb_mimp_secreted.gff ./
intersectBed -wa -a Fus2_canu_contigs_unmasked_pathogen2.vcf -b Fus2_canu_new_CAZY_secreted.gff >test.vcf
intersectBed -wa -a Fus2_canu_contigs_unmasked_pathogen2.vcf -b F.oxysporum_fsp_cepae_Fus2_canu_new_EffectorP_secreted.gff >test2.vcf
#OK, now only focus on FOC125,,FOCFus2,,FOC55,,FOCA23
cd $input/patho
for a in *.gff
do 
intersectBed -wa -a Fus2_canu_contigs_unmasked_patho_filtered.recode_annotated.vcf -b $a >${a%.gff}.vcf
done
#Slim down the file to retain only those, or minus 55.


#Check if any SNPs retained.
cd $input/patho-no55
for a in *.gff
do 
intersectBed -wa -a Fus2_canu_contigs_unmasked_patho_no55_filtered.recode.vcf -b $a >${a%.gff}.vcf
done