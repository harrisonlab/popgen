#!/bin/bash
wdir=/home/sobczm/popgen/codon/blast
scripts=/home/sobczm/bin/popgen/codon
dag=/home/sobczm/bin/DAGCHAINER

for z in *filtered
do
cp $z ./dagchainer
done

for z in *gene_table.txt
do
cp $z ./dagchainer
done

#blast_to_dagchainer
cd ./dagchainer
python $scripts/blast_to_dagchainer.py 125_final_genes_combined.cdna_one.fasta_vs_125_final_genes_combined.cdna_one_nucl.db_filtered \
125_final_genes_appended_gene_table.txt
python $scripts/blast_to_dagchainer.py 55_final_genes_combined.cdna_one.fasta_vs_55_final_genes_combined.cdna_one_nucl.db_filtered \
55_final_genes_appended_gene_table.txt
python $scripts/blast_to_dagchainer.py A1-2_final_genes_combined.cdna_one.fasta_vs_A1-2_final_genes_combined.cdna_one_nucl.db_filtered \
A1-2_final_genes_appended_gene_table.txt
python $scripts/blast_to_dagchainer.py A13_final_genes_combined.cdna_one.fasta_vs_A13_final_genes_combined.cdna_one_nucl.db_filtered \
A13_final_genes_appended_gene_table.txt
python $scripts/blast_to_dagchainer.py A23_final_genes_combined.cdna_one.fasta_vs_A23_final_genes_combined.cdna_one_nucl.db_filtered \
A23_final_genes_appended_gene_table.txt
python $scripts/blast_to_dagchainer.py A28_final_genes_combined.cdna_one.fasta_vs_A28_final_genes_combined.cdna_one_nucl.db_filtered \
A28_final_genes_appended_gene_table.txt
python $scripts/blast_to_dagchainer.py CB3_final_genes_combined.cdna_one.fasta_vs_CB3_final_genes_combined.cdna_one_nucl.db_filtered \
CB3_final_genes_appended_gene_table.txt
python $scripts/blast_to_dagchainer.py D2_final_genes_combined.cdna_one.fasta_vs_D2_final_genes_combined.cdna_one_nucl.db_filtered \
D2_final_genes_appended_gene_table.txt
python $scripts/blast_to_dagchainer.py Fus2_final_genes_combined.cdna_one.fasta_vs_Fus2_final_genes_combined.cdna_one_nucl.db_filtered \
Fus2_final_genes_appended_gene_table.txt
python $scripts/blast_to_dagchainer.py HB6_final_genes_combined.cdna_one.fasta_vs_HB6_final_genes_combined.cdna_one_nucl.db_filtered \
HB6_final_genes_appended_gene_table.txt
python $scripts/blast_to_dagchainer.py PG_final_genes_combined.cdna_one.fasta_vs_PG_final_genes_combined.cdna_one_nucl.db_filtered \
PG_final_genes_appended_gene_table.txt
python $scripts/blast_to_dagchainer.py proliferatum_final_genes_combined.cdna_one.fasta_vs_proliferatum_final_genes_combined.cdna_one_nucl.db_filtered \
proliferatum_final_genes_appended_gene_table.txt

#Run DAGchainer to detect duplications
for a in *dagchainer
do
$dag/run_DAG_chainer.pl -i $a -s -M 1000 -D 10000000 -x 10 -A 1
done

#-M Maximum match score (default: 50) otherwise, -log(evalue)
#-D maximum distance allowed between two matches in basepairs. (default: 200000)
#-A Minium number of Aligned Pairs (default: 6)

#Remove redundant lines and DAGChainer comment lines from the output
for a in *.aligncoords
do
sort $a | uniq | grep -Ev '^#' > ${a}f
done
