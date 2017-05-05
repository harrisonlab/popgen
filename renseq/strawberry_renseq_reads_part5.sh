#!/bin/bash
scripts=/home/sobczm/bin/popgen/renseq
input=/home/sobczm/popgen/renseq/strawberry/reads/Helen_Bates_EMR.RH.ENQ-1704.A.01/analysis/Rpf2

#Final Rpf2 analysis to detect variants in Hapil and RG.
#Prefilter BLAST results to only retain sequences with 90% homology to Rpf2 in vesca
#and miniumum 1000 bp alignment.

for a in *.db
do
awk '$3>90 && $4 >=1000 {print $0}' $a >${a}_filtered
done

for a in *filtered
do
awk -F $"\t" '$15=="plus" {print $0}' $a | cut -f2 | sort | uniq >"${a%.*}_pos.txt"
awk -F $"\t" '$15=="minus" {print $0}' $a | cut -f2 | sort | uniq >"${a%.*}_neg.txt"
done

#Extract FASTA sequences and append cultivar id
#Plus strand contigs
python $scripts/keep_list_genes2.py F.vesca_Rpf2_vs_F06.contigs_nucl_pos.txt F06.contigs.fasta No
python $scripts/keep_list_genes2.py F.vesca_Rpf2_vs_G06.contigs_nucl_pos.txt G06.contigs.fasta No
python $scripts/keep_list_genes2.py F.vesca_Rpf2_vs_S1_ccs_3_99_nucl_pos.txt S1_ccs_3_99.fasta No
python $scripts/keep_list_genes2.py F.vesca_Rpf2_vs_S2_ccs_3_99_nucl_pos.txt S2_ccs_3_99.fasta No
python $scripts/keep_list_genes2.py F.vesca_Rpf2_vs_S1_HGAP_polished_assembly_nucl_pos.txt S1_HGAP_polished_assembly.fasta No
python $scripts/keep_list_genes2.py F.vesca_Rpf2_vs_S2_HGAP_polished_assembly_nucl_pos.txt S2_HGAP_polished_assembly.fasta No
#Negative strand contigs
python $scripts/keep_list_genes2.py F.vesca_Rpf2_vs_F06.contigs_nucl_neg.txt F06.contigs.fasta Yes
python $scripts/keep_list_genes2.py F.vesca_Rpf2_vs_G06.contigs_nucl_neg.txt G06.contigs.fasta Yes
python $scripts/keep_list_genes2.py F.vesca_Rpf2_vs_S1_ccs_3_99_nucl_neg.txt S1_ccs_3_99.fasta Yes
python $scripts/keep_list_genes2.py F.vesca_Rpf2_vs_S2_ccs_3_99_nucl_neg.txt S2_ccs_3_99.fasta Yes
python $scripts/keep_list_genes2.py F.vesca_Rpf2_vs_S1_HGAP_polished_assembly_nucl_neg.txt S1_HGAP_polished_assembly.fasta Yes
python $scripts/keep_list_genes2.py F.vesca_Rpf2_vs_S2_HGAP_polished_assembly_nucl_neg.txt S2_HGAP_polished_assembly.fasta Yes
#Append population id
for a in *S1*pos.fasta *S1*neg.fasta *F06*pos.fasta *F06*neg.fasta 
do
sed -i 's/>/>Hapil_/g' $a
done

for a in *S2*pos.fasta *S2*neg.fasta *G06*pos.fasta *G06*neg.fasta 
do
sed -i 's/>/>RG_/g' $a
done

####Align with mafft
for a in *_neg.fasta
do
cat F.vesca_Rpf2.fasta $a ${a%_neg.fasta}_pos.fasta >${a%neg.fasta}all.fasta
done

