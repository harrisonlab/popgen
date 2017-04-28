#!/bin/bash
scripts=/home/sobczm/bin/popgen/renseq
input=/home/sobczm/popgen/renseq/strawberry/reads/Helen_Bates_EMR.RH.ENQ-1704.A.01

#Concatenate the raw reads
#Cultivar1
mkdir -p $input/assembly
cd $input/Raw_reads_S1/F06_1/Analysis_Results
for a in *subreads.fastqcd $input/assembly
do
cd $input/assembly
cat $a >> $input/assembly/F06_1_S1.fastq
done

#Cultivar2
cd $input/Raw_reads_S1/G06_1/Analysis_Results
for a in *subreads.fastq
do
cat $a >> $input/assembly/G06_1_S2.fastq
done

#Assemble each cultivar with Canu
cd $input/assembly
qsub $scripts/sub_canu_renseq.sh F06_1_S1.fastq 1m F06 F06_assembly
qsub $scripts/sub_canu_renseq.sh G06_1_S2.fastq 1m G06 G06_assembly

#Analyse HQ ROI, at 90% accuracy, min 1 full-length pass for NBS genes.
mkdir $input/assembly/nbs-parser
cp $input/Helen_Bates_EMR.RH.ENQ-1704.A.01_S1_MinFullPasses1_Accuracy90/data/reads_of_insert.fasta \
$input/assembly/nbs-parser/F06_1_S1_roi_fl1_a90.fasta
sh $scripts/sub_nlrparser.sh F06_1_S1_roi_fl1_a90.fasta
cp $input/Helen_Bates_EMR.RH.ENQ-1704.A.01_S2_MinFullPasses1_Accuracy90/data/reads_of_insert.fasta \
$input/assembly/nbs-parser/G06_1_S2_roi_fl1_a90.fasta
sh $scripts/sub_nlrparser.sh G06_1_S2_roi_fl1_a90.fasta

#Analyse the assemblies for NBS genes.
cd $input/assembly/F06_assembly
cat F06.contigs.fasta >> ../F06.assembly.fasta 
cat F06.unassembled.fasta >> ../F06.assembly.fasta 

cd $input/assembly/G06_assembly
cat G06.contigs.fasta >> ../G06.assembly.fasta 
cat G06.unassembled.fasta >> ../G06.assembly.fasta 
cp $input/assembly/*.assembly.fasta $input/assembly/nbs-parser

cd $input/assembly/nbs-parser
sh $scripts/sub_nlrparser.sh F06.assembly.fasta 
sh $scripts/sub_nlrparser.sh G06.assembly.fasta

#Cluster the assemblies at 98% ID
usearch=/home/sobczm/bin/usearch/usearch9.0.2132_i86linux32
for fasta in G06.assembly.fasta F06.assembly.fasta
do
id=0.98
$usearch -cluster_fast $fasta -id $id -sort length -centroids ${fasta%.fasta}_${id}_centroids.fasta
done 

#Create blast databases out of the assemblies
for assembly in G06.assembly.fasta F06.assembly.fasta
do
makeblastdb -in $assembly -input_type fasta -dbtype nucl \
-title "${assembly%.*}"_nucl.db -parse_seqids -out "${assembly%.*}"_nucl.db
done

#Search the assemblies for Rpf2
for db in F06.assembly_nucl.db G06.assembly_nucl.db
do
blastn -outfmt 6 -num_threads 1 -max_target_seqs 100 -evalue 0.0000000001 -query F.vesca_Rpf2.fasta -db $db >> F.vesca_Rpf2_vs_$db
done

#Make db of just the assembled contigs in each assembly
for assembly in G06.contigs.fasta F06.contigs.fasta
do
makeblastdb -in $assembly -input_type fasta -dbtype nucl \
-title "${assembly%.*}"_nucl.db -parse_seqids -out "${assembly%.*}"_nucl.db
done

for db in F06.contigs_nucl.db G06.contigs_nucl.db
do
blastn -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand"  -num_threads 1 -max_target_seqs 100 -evalue 0.0000000001 -query F.vesca_Rpf2.fasta -db $db >> F.vesca_Rpf2_vs_$db
done

#Fish out all individual ccs reads containing Rpf2 with blast
for reads in F06_1_S1_roi_fl1_a90.fasta G06_1_S2_roi_fl1_a90.fasta
do
makeblastdb -in $reads -input_type fasta -dbtype nucl \
-title "${reads%.*}"_nucl.db -parse_seqids -out "${reads%.*}"_nucl.db
done

for db in F06_1_S1_roi_fl1_a90_nucl.db G06_1_S2_roi_fl1_a90_nucl.db
do
blastn -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand" -num_threads 1 -max_target_seqs 100 -evalue 0.0000000001 -query -db $db >> F.vesca_Rpf2_vs_$db
done


###Alignment of individual Rpf2-like ccs reads
#Separate the gene list into those on the positive and negative strands.
for a in F.vesca_Rpf2_vs_F06_1_S1_roi_fl1_a90_nucl.db F.vesca_Rpf2_vs_G06_1_S2_roi_fl1_a90_nucl.db
do
awk -F $"\t" '$13=="plus" {print $0}' $a | cut -f2 | sort | uniq >"${a%.*}_pos.txt"
awk -F $"\t" '$13=="minus" {print $0}' $a | cut -f2 | sort | uniq >"${a%.*}_neg.txt"
done

#Extract FASTA sequences and append cultivar id
#Plus strand contigs
python $scripts/keep_list_genes2.py F.vesca_Rpf2_vs_F06_1_S1_roi_fl1_a90_nucl_pos.txt F06_1_S1_roi_fl1_a90.fasta No
python $scripts/keep_list_genes2.py F.vesca_Rpf2_vs_G06_1_S2_roi_fl1_a90_nucl_pos.txt G06_1_S2_roi_fl1_a90.fasta No

#Negative strand contigs - reverse complement
python $scripts/keep_list_genes2.py F.vesca_Rpf2_vs_F06_1_S1_roi_fl1_a90_nucl_neg.txt F06_1_S1_roi_fl1_a90.fasta Yes
python $scripts/keep_list_genes2.py F.vesca_Rpf2_vs_G06_1_S2_roi_fl1_a90_nucl_neg.txt G06_1_S2_roi_fl1_a90.fasta Yes

sed -i 's/>/>Hap_/g' F.vesca_Rpf2_vs_F06_1_S1_roi_fl1_a90_nucl_neg.fasta
sed -i 's/>/>Hap_/g' F.vesca_Rpf2_vs_F06_1_S1_roi_fl1_a90_nucl_pos.fasta

sed -i 's/>/>RG_/g' F.vesca_Rpf2_vs_G06_1_S2_roi_fl1_a90_nucl_neg.fasta
sed -i 's/>/>RG_/g' F.vesca_Rpf2_vs_G06_1_S2_roi_fl1_a90_nucl_pos.fasta

####Align the ccs with mafft
cat F.vesca_Rpf2.fasta *_neg.fasta *_pos.fasta >F.vesca_Rpf2_fl1_a90_reads.fasta
qsub $scripts/sub_mafft.sh F.vesca_Rpf2_fl1_a90_reads.fasta


###########Do the same for contigs
for a in F.vesca_Rpf2_vs_F06.contigs_nucl.db F.vesca_Rpf2_vs_G06.contigs_nucl.db
do
awk -F $"\t" '$13=="plus" {print $0}' $a | cut -f2 | sort | uniq >"${a%.*}_pos.txt"
awk -F $"\t" '$13=="minus" {print $0}' $a | cut -f2 | sort | uniq >"${a%.*}_neg.txt"
done

python $scripts/keep_list_genes2.py F.vesca_Rpf2_vs_F06.contigs_nucl_pos.txt F06.contigs.fasta No
python $scripts/keep_list_genes2.py F.vesca_Rpf2_vs_G06.contigs_nucl_pos.txt G06.contigs.fasta No

#Negative strand contigs - reverse complement
python $scripts/keep_list_genes2.py F.vesca_Rpf2_vs_F06.contigs_nucl_neg.txt F06.contigs.fasta Yes
python $scripts/keep_list_genes2.py F.vesca_Rpf2_vs_G06.contigs_nucl_neg.txt G06.contigs.fasta Yes

sed -i 's/>/>Hap_/g' F.vesca_Rpf2_vs_F06.contigs_nucl_neg.fasta
sed -i 's/>/>Hap_/g' F.vesca_Rpf2_vs_F06.contigs_nucl_pos.fasta

sed -i 's/>/>RG_/g' F.vesca_Rpf2_vs_G06.contigs_nucl_neg.fasta
sed -i 's/>/>RG_/g' F.vesca_Rpf2_vs_G06.contigs_nucl_pos.fasta

cat F.vesca_Rpf2.fasta F.vesca_Rpf2_vs_F06.contigs*.fasta F.vesca_Rpf2_vs_G06.contigs*.fasta >F.vesca_Rpf2_contigs.fasta
qsub $scripts/sub_mafft.sh F.vesca_Rpf2_contigs.fasta