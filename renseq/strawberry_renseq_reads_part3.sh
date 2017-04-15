#!/bin/bash
scripts=/home/sobczm/bin/popgen/renseq
input=/home/sobczm/popgen/renseq/strawberry/reads/Helen_Bates_EMR.RH.ENQ-1704.A.01
cd $input/analysis
#####Evaluation of PacBio reads and contigs from RenSeq

###How many genes from Rob's baits are present in the Ren-seq output??
#GFF subset: mRNA
awk '$3 == "mRNA" {print $0}' vesca_v1.1_nblrrs_augustus.gff >vesca_v1.1_nblrrs_augustus_mrna.gff
awk '$3 == "CDS" {print $0}' vesca_v1.1_nblrrs_augustus.gff >vesca_v1.1_nblrrs_augustus_cds.gff
##Extract the sequences of the genes to be used as blast database (CDS + 5'/3' UTR)
genome_seq=/home/sobczm/popgen/renseq/strawberry/genome/fvesca_v1.1_all.fa
gffread=/home/sobczm/bin/gffread/gffread/gffread
$gffread -w vesca_v1.1_nblrrs_augustus.fasta -g $genome_seq vesca_v1.1_nblrrs_augustus.gff
$gffread -w vesca_v1.1_nblrrs_augustus_mrna.fasta -g $genome_seq vesca_v1.1_nblrrs_augustus_mrna.gff
$gffread -w vesca_v1.1_nblrrs_augustus_cds.fasta -g $genome_seq vesca_v1.1_nblrrs_augustus_cds.gff

#Make databases of mRNA and CDS sequences
makeblastdb -in vesca_v1.1_nblrrs_augustus_cds.fasta -input_type fasta -dbtype nucl \
-title vesca_v1.1_nblrrs_augustus_cds_nucl.db -parse_seqids -out vesca_v1.1_nblrrs_augustus_cds_nucl.db

makeblastdb -in vesca_v1.1_nblrrs_augustus_mrna.fasta -input_type fasta -dbtype nucl \
-title vesca_v1.1_nblrrs_augustus_mrna_nucl.db -parse_seqids -out vesca_v1.1_nblrrs_augustus_mrna_nucl.db

#Run blast search against the CDS database.
db=vesca_v1.1_nblrrs_augustus_cds_nucl.db
for ass in $input/assembly/G06.contigs.fasta $input/assembly/F06.contigs.fasta $input/assembly/F06_1_S1_roi_fl1_a90.fasta $input/assembly/G06_1_S2_roi_fl1_a90.fasta
do
blastn -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand"  -num_threads 1 -max_target_seqs 100 -evalue 0.0000000001 -query $ass -db $db >> $(basename $ass)_vs_$db
done