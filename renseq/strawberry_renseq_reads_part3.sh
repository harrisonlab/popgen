#!/bin/bash
scripts=/home/sobczm/bin/popgen/renseq
input=/home/sobczm/popgen/renseq/strawberry/reads/Helen_Bates_EMR.RH.ENQ-1704.A.01
cd $input/analysis
#####Evaluation of PacBio reads and contigs from RenSeq

###How many genes from Rob's baits are present in the Ren-seq output??
#GFF subset: mRNA and CDS
awk '$3 == "mRNA" {print $0}' vesca_v1.1_nblrrs_augustus.gff >vesca_v1.1_nblrrs_augustus_mrna.gff
awk '$3 == "CDS" {print $0}' vesca_v1.1_nblrrs_augustus.gff >vesca_v1.1_nblrrs_augustus_cds.gff
##Extract the sequences of the genes to be used as blast database (CDS + 5'/3' UTR)
genome_seq=/home/sobczm/popgen/renseq/strawberry/genome/fvesca_v1.1_all.fa
gffread=/home/sobczm/bin/gffread/gffread/gffread
$gffread -w vesca_v1.1_nblrrs_augustus.fasta -g $genome_seq vesca_v1.1_nblrrs_augustus.gff
$gffread -w vesca_v1.1_nblrrs_augustus_mrna.fasta -g $genome_seq vesca_v1.1_nblrrs_augustus_mrna.gff
$gffread -w vesca_v1.1_nblrrs_augustus_cds.fasta -g $genome_seq vesca_v1.1_nblrrs_augustus_cds.gff
##########Main vesca GFF annotation
genome_gff=/home/sobczm/popgen/renseq/strawberry/genome/Fragaria_vesca_v1.1.a2.gff3 
#GFF subset: mRNA and CDS
awk '$3 == "mRNA" {print $0}' $genome_gff >Fragaria_vesca_v1.1.a2_mrna.gff
awk '$3 == "CDS" {print $0}' $genome_gff >Fragaria_vesca_v1.1.a2_cds.gff
##Extract the sequences of the genes to be used as blast database (CDS + 5'/3' UTR)
$gffread -w Fragaria_vesca_v1.1.a2.fasta -g $genome_seq $genome_gff 
$gffread -w Fragaria_vesca_v1.1.a2_mrna.fasta -g $genome_seq Fragaria_vesca_v1.1.a2_mrna.gff
$gffread -w Fragaria_vesca_v1.1.a2_cds.fasta -g $genome_seq Fragaria_vesca_v1.1.a2_cds.gff

#Make databases of mRNA and CDS sequences (genes selected for RenSeq and all genes in vesca 1.1 genome) as well as baits sequences
makeblastdb -in vesca_v1.1_nblrrs_augustus_cds.fasta -input_type fasta -dbtype nucl \
-title vesca_v1.1_nblrrs_augustus_cds_nucl.db -parse_seqids -out vesca_v1.1_nblrrs_augustus_cds_nucl.db

makeblastdb -in vesca_v1.1_nblrrs_augustus_mrna.fasta -input_type fasta -dbtype nucl \
-title vesca_v1.1_nblrrs_augustus_mrna_nucl.db -parse_seqids -out vesca_v1.1_nblrrs_augustus_mrna_nucl.db

makeblastdb -in Fragaria_vesca_v1.1.a2_mrna.fasta -input_type fasta -dbtype nucl \
-title Fragaria_vesca_v1.1.a2_mrna.fasta_nucl.db -parse_seqids -out Fragaria_vesca_v1.1.a2_mrna.fasta_nucl.db

#Needed to remove redundant sequences showing more than 1 sequence per ID
grep ">" Fragaria_vesca_v1.1.a2_cds.fasta | cut -f1 -d" " | sort | uniq -c | sort -n -k1,1
# 2 >augustus_masked-LG5-processed-gene-63_24-mRNA-1
#     2 >maker-LG0-augustus-gene-47_114-mRNA-1
#      2 >maker-LG2-augustus-gene-159_157-mRNA-1
#      2 >maker-LG2-augustus-gene-180_162-mRNA-1
#      2 >maker-LG2-augustus-gene-207_206-mRNA-1
#      2 >maker-LG4-augustus-gene-186_102-mRNA-1
#      2 >maker-LG4-augustus-gene-76_95-mRNA-1
#      2 >maker-LG6-snap-gene-258_175-mRNA-1
#      3 >XS:temp5
#      3 >XS:temp6
#      4 >XS:temp2
#      4 >XS:temp4
#      5 >XS:temp1
#      6 >XS:temp3
grep ">" Fragaria_vesca_v1.1.a2_mrna.fasta | cut -f1 -d" " | sort | uniq -c | sort -n -k1,1
#3 >XS:temp5
#      3 >XS:temp6
#      4 >XS:temp2
#      4 >XS:temp4
#      5 >XS:temp1
#      6 >XS:temp3

makeblastdb -in Fragaria_vesca_v1.1.a2_cds.fasta -input_type fasta -dbtype nucl \
-title Fragaria_vesca_v1.1.a2_mrna.fasta_nucl.db -parse_seqids -out Fragaria_vesca_v1.1.a2_cds.fasta_nucl.db

#Make database using FASTA file with removed duplicates.
makeblastdb -in Fragaria_vesca_v1.1.a2_cds_removed.fasta -input_type fasta -dbtype nucl \
-title Fragaria_vesca_v1.1.a2_cds.fasta_nucl.db -parse_seqids -out Fragaria_vesca_v1.1.a2_cds_removed.fasta_nucl.db

makeblastdb -in Fragaria_vesca_v1.1.a2_mrna_removed.fasta -input_type fasta -dbtype nucl \
-title Fragaria_vesca_v1.1.a2_mrna.fasta_nucl.db -parse_seqids -out Fragaria_vesca_v1.1.a2_mrna_removed.fasta_nucl.db

##Copied over putative bait sequences from /home/vicker/ananassa_renseq/final_bait_library_versions/2015-07-20_1x_masked_20K/mycroarray_vescaandiinumae/MYbaits-150720-Vickerstaff-Strawberry/probes-R4-final.fas
#Make Blast database
makeblastdb -in probes-R4-final.fas -input_type fasta -dbtype nucl \
-title probes-R4-final.nucl.db -parse_seqids -out probes-R4-final.nucl.db

#Run blast search against the CDS database of chosen genes as well as the entire strawberry genome CDS database as well as baits database.
for db in probes-R4-final.nucl.db Fragaria_vesca_v1.1.a2_cds_removed.fasta_nucl.db vesca_v1.1_nblrrs_augustus_cds_nucl.db
do
for ass in $input/assembly/G06.contigs.fasta $input/assembly/F06.contigs.fasta $input/assembly/F06_1_S1_roi_fl1_a90.fasta $input/assembly/G06_1_S2_roi_fl1_a90.fasta
do
    blastn -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand"  -num_threads 1 -max_target_seqs 100 -evalue 0.0000000001 -query $ass -db $db >> $(basename $ass)_vs_$db
done
done
#Generate a file with sequence lengths for CDS input files.
#Note for future reference: BLAST can automatically output these values if fields slen and qlen specified.
for a in *.fasta
do
python $scripts/write_seq_length.py $a
done

