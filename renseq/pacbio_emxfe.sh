###Copy the input files to triticum
datadir=/data/seq_data/external/20171110_EmilyFenellaRenseq_pacbio

##Carried out on blacklace11
###Trim Illumina adaptors from raw PacBio reads using a script from Giolai et al. (2016)
input=/home/sobczm/popgen/renseq/strawberry/reads/pacbio_emxfe
scripts=/home/sobczm/bin/popgen/renseq

tar -zvxf $datadir/Helen_Bates_EMR_RH_ENQ-1704_A_01_S3.tar.gz
tar -zxvf $datadir/Helen_Bates_EMR_RH_ENQ-1704_A_01_S4.tar.gz 

cd $input/Helen_Bates_EMR_RH_ENQ-1704_A_01_S3/2017_10_25_PSEQ1569_406/SAM32001_PRO1514_S3_HMWDNA_Emily/raw_reads/D01_1/Analysis_Results
for a in *.bax.h5; do $scripts/trim-h5.py $a; done

cd $input/Helen_Bates_EMR_RH_ENQ-1704_A_01_S4/2017_10_25_PSEQ1569_406/SAM32002_PRO1514_S4_HMWDNA_Fenella/raw_reads/E01_1/Analysis_Results
for a in *.bax.h5; do $scripts/trim-h5.py $a; done

scp -r $input/Helen_Bates_EMR_RH_ENQ-1704_A_01_S3 sobczm@10.1.10.170:/data/projects/sobczm/pacbio_emxfe
scp -r $input/Helen_Bates_EMR_RH_ENQ-1704_A_01_S4 sobczm@10.1.10.170:/data/projects/sobczm/pacbio_emxfe

#Prepare CCS at 99% accuracy from minimum 3 reads pass
bb=/home/sobczm/bin/pitchfork/workspace/bax2bam/bin
ccs=/home/sobczm/bin/unanimity/build/ccs 

#S3
cd /data/projects/sobczm/pacbio_emxfe/Helen_Bates_EMR_RH_ENQ-1704_A_01_S3/Analysis_Results
LD_LIBRARY_PATH=/home/sobczm/bin/pitchfork/deployment/lib $bb/bax2bam -o S3.bam m171026_073505_42165_c101248002550000001823283011021703_s1_p0.1.bax.h5 m171026_073505_42165_c101248002550000001823283011021703_s1_p0.2.bax.h5 m171026_073505_42165_c101248002550000001823283011021703_s1_p0.3.bax.h5
$ccs --minLength=1000 --minPredictedAccuracy=0.99 S3.bam.subreads.bam S3_ccs_3_99.bam

cd /data/projects/sobczm/pacbio_emxfe/Helen_Bates_EMR_RH_ENQ-1704_A_01_S4/Analysis_Results
LD_LIBRARY_PATH=/home/sobczm/bin/pitchfork/deployment/lib $bb/bax2bam -o S4.bam m171026_120100_42165_c101248002550000001823283011021704_s1_p0.1.bax.h5 m171026_120100_42165_c101248002550000001823283011021704_s1_p0.2.bax.h5 m171026_120100_42165_c101248002550000001823283011021704_s1_p0.3.bax.h5
$ccs --minLength=1000 --minPredictedAccuracy=0.99 S4.bam.subreads.bam S4_ccs_3_99.bam

#Copy the files over to the EMR cluster, convert to fastq and then fasta.
cd $input/analysis
scp -r sobczm@10.1.10.170:/data/projects/sobczm/pacbio_emxfe/Helen_Bates_EMR_RH_ENQ-1704_A_01_S3/Analysis_Results/S3_ccs_3_99.bam ./
scp -r sobczm@10.1.10.170:/data/projects/sobczm/pacbio_emxfe/Helen_Bates_EMR_RH_ENQ-1704_A_01_S4/Analysis_Results/S4_ccs_3_99.bam ./

bedtools bamtofastq -i S3_ccs_3_99.bam -fq S3_ccs_3_99.fastq
bedtools bamtofastq -i S4_ccs_3_99.bam -fq S4_ccs_3_99.fastq
/home/sobczm/bin/seqtk/seqtk seq -a S3_ccs_3_99.fastq > S3_ccs_3_99.fasta
/home/sobczm/bin/seqtk/seqtk seq -a S4_ccs_3_99.fastq > S4_ccs_3_99.fasta

#Concatenate raw reads and trim adaptors
for reads in $input/Helen_Bates_EMR_RH_ENQ-1704_A_01_S3/2017_10_25_PSEQ1569_406/SAM32001_PRO1514_S3_HMWDNA_Emily/raw_reads/D01_1/Analysis_Results/*subreads.fastq
do
cat $reads >> $input/assembly/D06_1_S3.fastq
done

for reads in $input/Helen_Bates_EMR_RH_ENQ-1704_A_01_S4/2017_10_25_PSEQ1569_406/SAM32002_PRO1514_S4_HMWDNA_Fenella/raw_reads/E01_1/Analysis_Results/*subreads.fastq
do
cat $reads >> $input/assembly/E06_1_S4.fastq
done

#Tar the reads and submit for QC check and adaptor trimming.
for a in *.fastq; do gzip $a; done
qsub $scripts/sub_read_qc_single_pacbio.sh $input/assembly/D06_1_S3.fastq.gz
qsub $scripts/sub_read_qc_single_pacbio.sh $input/assembly/E06_1_S4.fastq.gz

#Canu assembly
qsub $scripts/sub_canu_renseq.sh D06_1_S3_trim.fq.gz 1m D06_S3 D06_assembly
qsub $scripts/sub_canu_renseq.sh E06_1_S4_trim.fq.gz 1m E06_S4 E06_assembly

#copy assembly to the analysis folder
cp $input/assembly/D06_assembly/D06_S3.contigs.fasta $input/analysis/D06_S3_Emily.contigs.fasta
cp $input/assembly/E06_assembly/E06_S4.contigs.fasta $input/analysis/D06_S4_Fenella.contigs.fasta

cd $input/analysis

#Blast the raw reads to establish hits R genes.
#Copy the blast db and bait sequences
cp -r /home/sobczm/popgen/renseq/strawberry/reads/Helen_Bates_EMR.RH.ENQ-1704.A.01/analysis/vesca_v1.1_nblrrs_augustus_cds_nucl.db* ./
cp -r /home/sobczm/popgen/renseq/strawberry/reads/Helen_Bates_EMR.RH.ENQ-1704.A.01/analysis/probes-R4-final.fas ./
cp -r /home/sobczm/popgen/renseq/strawberry/reads/Helen_Bates_EMR.RH.ENQ-1704.A.01/analysis/vesca_v1.1_nblrrs_augustus_cds.fasta ./

for a in D06_S3_Emily.contigs.fasta D06_S4_Fenella.contigs.fasta S3_ccs_3_99.fasta S4_ccs_3_99.fasta
do 
#Remove spaces in FASTA header, and substitute commas with _, as required by BLAST suite.
sed 's, ,_,g' -i $a
sed 's/,/_/g' -i $a
perl /home/sobczm/bin/popgen/other/assemblathon_stats.pl $a > ${a%.fasta}.stat
qsub $scripts/sub_blastn_renseq.sh $a vesca_v1.1_nblrrs_augustus_cds_nucl.db
sh $scripts/sub_nlrparser.sh $(basename $a)
done

#Create blast databases out of the assemblies
for assembly in D06_S3_Emily.contigs.fasta D06_S4_Fenella.contigs.fasta S3_ccs_3_99.fasta S4_ccs_3_99.fasta
do
makeblastdb -in $assembly -input_type fasta -dbtype nucl \
-title "${assembly%.*}"_nucl.db -parse_seqids -out "${assembly%.*}"_nucl.db
done

#Search the baits sequences and NLR genes against the databases of sequences of Ren-Seq reads/assemblies.
for db in D06_S3_Emily.contigs_nucl.db D06_S4_Fenella.contigs_nucl.db S3_ccs_3_99_nucl.db S4_ccs_3_99_nucl.db
do
blastn -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand"  -num_threads 1 -max_target_seqs 1000000000 -evalue 0.0000000001 -query probes-R4-final.fas -db $db >> probes-R4-final.fas_vs_$db
qsub $scripts/sub_blastn_renseq.sh vesca_v1.1_nblrrs_augustus_cds.fasta $db
done

#Generate a file with sequence lengths for CDS input files.
#Note for future reference: BLAST can automatically output these values if fields slen and qlen specified.
for a in *.fasta
do
python $scripts/write_seq_length.py $a
done

for sample in D06_S3_Emily.contigs D06_S4_Fenella.contigs S3_ccs_3_99 S4_ccs_3_99
do
#Identify RBB matches.
python /home/sobczm/bin/popgen/other/rbb.py ${sample}.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db vesca_v1.1_nblrrs_augustus_cds.fasta_vs_${sample}_nucl.db >${sample}_vesca_v1.1_nblrrs_augustus_cds_rbb
#Count up the length of CDS, 5' UTR, 3' UTR using the original renseq vs vesca CDS BLAST output
python $scripts/blast_matches_summary2.py vesca_v1.1_nblrrs_augustus_cds_lengths.txt ${sample}_lengths.txt ${sample}.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db > ${sample}.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_cds
#Filter those results to only retain the RBB hits.
python $scripts/filter_by_two_columns.py ${sample}_vesca_v1.1_nblrrs_augustus_cds_rbb ${sample}.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_cds > ${sample}.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_cds_rbb_filt
python $scripts/filter_by_two_columns.py ${sample}_vesca_v1.1_nblrrs_augustus_cds_rbb ${sample}.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_5utr > ${sample}.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_5utr_rbb_filt 
python $scripts/filter_by_two_columns.py ${sample}_vesca_v1.1_nblrrs_augustus_cds_rbb ${sample}.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_3utr > ${sample}.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_3utr_rbb_filt 
done

#Plot figures summarising the results above in R.
Rscript --vanilla $scripts/stats_cds_emfe.R D06_S3_Emily.contigs.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_cds_rbb_filt D06_S4_Fenella.contigs.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_cds_rbb_filt EMFE_canu_vs_vesca_v1.1_nblrrs_augustus_cds
Rscript --vanilla $scripts/stats_cds_emfe.R S3_ccs_3_99.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_cds_rbb_filt S4_ccs_3_99.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_cds_rbb_filt EMFE_CCS_vs_vesca_v1.1_nblrrs_augustus_cds

for utr in 3 5
do
Rscript --vanilla $scripts/stats_utr_emfe.R D06_S3_Emily.contigs.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_${utr}utr_rbb_filt D06_S4_Fenella.contigs.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_${utr}utr_rbb_filt EMFE_canu_vs_vesca_v1.1_nblrrs_augustus_$utr $utr
Rscript --vanilla $scripts/stats_utr_emfe.R S3_ccs_3_99.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_${utr}utr_rbb_filt S4_ccs_3_99.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_${utr}utr_rbb_filt EMFE_CCS_vs_vesca_v1.1_nblrrs_augustus_$utr $utr
done

#Convert all to png
for my_pdf in *.pdf
do
convert -verbose -density 500 "${my_pdf}" "${my_pdf%.*}.png"
done

#Quick analysis of the NLR parser output.
for my_table in *.tsv
do
echo $my_table >>nlr_summary.txt
echo "complete" >>nlr_summary.txt
cat $my_table | awk '$3=="complete" {print $0}' | wc -l >>nlr_summary.txt
echo "partial" >>nlr_summary.txt
cat $my_table | awk '$3=="partial" {print $0}' | wc -l >>nlr_summary.txt
echo "pseudogene" >>nlr_summary.txt
cat $my_table | awk '$3=="pseudogene" {print $0}' | wc -l >>nlr_summary.txt
done