#!/bin/bash
input=/home/sobczm/popgen/rnaseq
scripts=/home/sobczm/bin/popgen/rnaseq

##Run Blastx annotation against TAIR.
#Create BLAST target database
cd $input/annotation
makeblastdb -in TAIR10_pep_20110103_representative_gene_model_updated.fasta -input_type fasta -dbtype prot \
-title TAIR10_cds_prot.db -parse_seqids -out TAIR10_cds_prot.db

perl $scripts/run_blast_reciprocal.pl Fragaria_vesca_v1.1.a2_cds.fasta TAIR10_cds_prot.db 
perl $scripts/run_blast_reciprocal.pl Fragaria_vesca_v1.1.a2_cds_ns.fasta TAIR10_cds_prot.db 

#InterProScan annotation of translated proteins
#python $HOME/bin/popgen/clock/dn_ds/remove_terminal_stop.py Fragaria_vesca_v1.1.a2_cds.fasta
emboss=/home/armita/prog/emboss/EMBOSS-4.0.0/bin
$emboss/transeq Fragaria_vesca_v1.1.a2_cds.fasta -outseq Fragaria_vesca_v1.1.a2_proteins.fasta
for a in Fragaria_vesca_v1.1.a2_proteins.fasta
do
b=$( echo $a | sed -e 's/.fasta/\_nostop.fasta/' )
sed 's/\*//g' $a > $b
qsub /home/sobczm/bin/popgen/renseq/sub_interproscan.sh $b
done

#CDS lengths to use for FPKM calculations
python /home/sobczm/bin/popgen/renseq/write_seq_length.py Fragaria_vesca_v1.1.a2_cds.fasta
