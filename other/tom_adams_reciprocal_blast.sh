#!/bin/bash
scripts=/home/sobczm/bin/popgen/other
input=/home/sobczm/popgen/other/phytophthora/RBB
cd $input
##Set up directory with required fasta files
#Copy relevant assemblies
frag_dir=/home/groups/harrisonlab/project_files/phytophthora_fragariae
#Race1 genomes
for Strain in Bc1 Nov5
do
    Assembly=repeat_masked/P.fragariae/$Strain/deconseq_Paen_repmask/"$Strain"_contigs_unmasked.fa
    cp $frag_dir/$Assembly $input
done

#Race 3 genomes
for Strain in Nov71 Nov9 Nov27
do
    Assembly=repeat_masked/P.fragariae/$Strain/deconseq_Paen_repmask/"$Strain"_contigs_unmasked.fa
    cp $frag_dir/$Assembly $input
done
#Bc16 genome 
cp $frag_dir/gene_pred/annotation/P.fragariae/Bc16/Bc16_genes_incl_ORFeffectors.gene.fasta $input

#Copy the list of RXLRs.
cp $frag_dir/analysis/Reciprocal_BLAST/Bc16_expressed_RxLR_headers_parsed.txt $input

#Make nucleotide BLAST db of all genomes.
for assembly in *.fa*
do
makeblastdb -in $assembly -input_type fasta -dbtype nucl -title "${assembly%.*}"_nucl.db -parse_seqids -out "${assembly%.*}"_nucl.db
done

#Run BLAST for all race 1 and race 2 assemblies against BC16...
db=Bc16_genes_incl_ORFeffectors.gene_nucl.db
for query in *.fa
do
blastn -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand"  -num_threads 1 -max_target_seqs 100 -evalue 0.0000000001 -query $query -db $db >> ${query}_vs_${db}
done

#And vice versa
query=Bc16_genes_incl_ORFeffectors.gene.fasta
for db in Bc1_contigs_unmasked_nucl.db Nov27_contigs_unmasked_nucl.db Nov5_contigs_unmasked_nucl.db Nov71_contigs_unmasked_nucl.db Nov9_contigs_unmasked_nucl.db
do
blastn -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand"  -num_threads 1 -max_target_seqs 100 -evalue 0.0000000001 -query $query -db $db >> ${query}_vs_${db}
done