#!/bin/bash
wdir=/home/sobczm/popgen/codon/blast/dagchainer/testing
scripts=/home/sobczm/bin/popgen/codon

############################################################################
# Analysis of duplication levels in the Fus2 genome, as it is the only non-frag
# mented genome available.
#############################################################################
#Seperating duplications into two classes: segmental and tandem

#First of all, using a maximum number of 5 intermittent genes to classify a
#gene as being a result of tandem duplication
cd $wdir
gene_table=Fus2_final_genes_appended_gene_table.txt
dagchainer_blast=Fus2_final_genes_combined.cdna_one.fasta_vs_Fus2_final_genes_combined.cdna_one_nucl.db_filtered_dagchainer.aligncoordsf
mkdir 5genes && cd 5genes
$scripts/detect_duplications.py --b ../$dagchainer_blast --g ../$gene_table --o gene --t 5

#Now more, lenient with a maximum number of 10 intermittent genes to classify a
#gene as being a result of tandem duplication
cd $wdir
mkdir 10genes && cd 10genes
$scripts/detect_duplications.py --b ../$dagchainer_blast --g ../$gene_table --o gene --t 10

#This time around, using the distance between the genes as the marker to classify the
#duplications as tandem or segmental
cd $wdir
mkdir 10kbp && cd 10kbp
$scripts/detect_duplications.py --b ../$dagchainer_blast --g ../$gene_table --o distance --t 10000

#And 20kbp for tandem duplications
cd $wdir
mkdir 20kbp && cd 20kbp
$scripts/detect_duplications.py --b ../$dagchainer_blast --g ../$gene_table --o distance --t 20000

#Conclusions
#High number of duplicated contigs at the ends of contigs: is this due to
#Contigs being shorter there, and therefore more likely to find a well-covered hit?
#Contigs 10 and 14, 16 are full of duplications, as other lineage-specific contigs.

#Now, focussing on the output summary table from the second run, where tandem duplications
#defined by max. 10 gene distance. How many potential effector genes are duplicated?

#Get a selection of summary table outputs containing different
#effector subsets
effectors=/home/sobczm/popgen/input/effectors
summary_table=Fus2_final_genes_combined.cdna_one.fasta_vs_Fus2_final_genes_combined.cdna_one_nucl.db_filtered_dagchainer_summary

cd $wdir/10genes
#13 out of 940 CAZY genes duplicated. All segmental duplications (in total 1355 seg. duplications detected across the genome).
#So seeing slight under-enrichment compared to the rest of the genome.
python $scripts/match_ids.py $effectors/CAZY/all/Fus2_canu_new_CAZY_gene_table.txt $summary_table
#9 out of 386 secreted CAZY genes segmentally duplicated.
python $scripts/match_ids.py $effectors/CAZY/secreted/Fus2_canu_new_CAZY_secreted_gene_table.txt $summary_table
#11 out of 355 genes duplicated. Genes appearing in 5, 6, 7, 9 copies.
python $scripts/match_ids.py $effectors/effectorP/secreted/F.oxysporum_fsp_cepae_Fus2_canu_new_EffectorP_secreted_gene_table.txt \
$summary_table

#75 genes duplicated, 5 of them tandemly (out of 150 genes annotated in this category, so enriched for duplication). Up to 10-16 copies of some
#genes seen.
python $scripts/match_ids.py $effectors/mimp/all/Fus2_canu_new_genes_in_2kb_mimp_gene_table.txt $summary_table
#7 genes found, all segmentally duplicated, out of 31 genes annotated.
python $scripts/match_ids.py $effectors/mimp/secreted/Fus2_canu_new_genes_in_2kb_mimp_secreted_gene_table.txt \
$summary_table

#Having investigated INTERPROSCAN and PFAM annotations of the putatively duplicated genes, it appears most of them are transposons.
#Going to exclude transposons now as well as unnanotated genes and plot the rest of the duplicated genes this time around.
#Excluded terms and domains:
#IPR000477 # Reverse transcriptase (RT) catalytic domain profile
#IPR012337 # Ribonuclease H domain
#IPR018289 # MULE transposase domain
#PF03221 # Tc5 transposase DNA-binding domain
#PF00078 # Reverse transcriptase (RNA-dependent)
#IPR025476 # Helitron helicase-like domain
#IPR008906 # hAT family C-terminal dimerisation region
#transpos*
#integrase

annotations=/home/sobczm/popgen/input/annotations
cd $annotations
cat Fus2_canu_new_gene_annotations.tab | sort > temp
grep 'IPR000477\|IPR004875\|IPR025476\|IPR012337\|IPR018289\|PF03221\|PF00078\|PF00078\|IPR008906\|transpos*\|integrase' temp | sort >temp2
comm -3 temp temp2 | cut -f 1 | awk '$0="Fus2_"$0'>Fus2_no_transposons.txt

#Keep only non-transposon genes
python $scripts/match_ids.py Fus2_no_transposons.txt Fus2_duplicated_genes_annotations.txt
cp Fus2_no_transposons_Fus2_duplicated_genes_annotations $wdir
cd $wdir
mkdir no_transposon
python $scripts/match_ids.py Fus2_no_transposons_Fus2_duplicated_genes_annotations Fus2_final_genes_combined.cdna_one.fasta_vs_Fus2_final_genes_combined.cdna_one_nucl.db_filtered_dagchainer.aligncoordsf
#Generates a file containing BLAST results for duplicated non-transposon and transposon genes each, respectively.
mv Fus2_no_transposons_Fus2_duplicated_genes_annotatio_Fus2_final_genes_combined.cdna_one.fasta_vs_Fus2_final_genes_combined.cdna_one_nucl.db_filtered_dagchainer_remainder transposon_duplications.aligncoordsf
mv Fus2_no_transposons_Fus2_duplicated_genes_annotatio_Fus2_final_genes_combined.cdna_one.fasta_vs_Fus2_final_genes_combined.cdna_one_nucl.db_filtered_dagchainer non-transposon_duplications.aligncoordsf
cp non-transposon_duplications.aligncoordsf ./no_transposon
cp Fus2_final_genes_appended_gene_table.txt ./no_transposon
cd $wdir/no_transposon

#Run duplication analysis, using max. 10 genes distance to define tandem duplications
$scripts/detect_duplications.py --b non-transposon_duplications.aligncoordsf --g Fus2_final_genes_appended_gene_table.txt --o gene --t 10
#The same pattern of duplication distribution detected as when using the transposon-containing dataset.
