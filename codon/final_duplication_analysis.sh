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
$popgen/codon/detect_duplications.py --b ../$dagchainer_blast --g ../$gene_table --o gene --t 5

#Now more, lenient with a maximum number of 10 intermittent genes to classify a
#gene as being a result of tandem duplication
cd $wdir
mkdir 10genes && cd 10genes
$popgen/codon/detect_duplications.py --b ../$dagchainer_blast --g ../$gene_table --o gene --t 10

#This time around, using the distance between the genes as the marker to classify the
#duplications as tandem or segmental
cd $wdir
mkdir 10kbp && cd 10kbp
$popgen/codon/detect_duplications.py --b ../$dagchainer_blast --g ../$gene_table --o distance --t 10000

#And 20kbp for tandem duplications
cd $wdir
mkdir 20kbp && cd 20kbp
$popgen/codon/detect_duplications.py --b ../$dagchainer_blast --g ../$gene_table --o distance --t 20000

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
