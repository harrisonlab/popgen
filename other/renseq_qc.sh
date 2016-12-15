#!/bin/bash
scripts=/home/sobczm/bin/popgen/other
input=/home/sobczm/popgen/other/renseq_redgauntlet
input_files=/home/miseq_data/minion/2016/MN18323/RENseqRG

cd $input

#Copying RedGauntlet genome assembly over
cp /home/groups/harrisonlab/project_files/fragaria_x_ananassa/octoseq_assembly/redgauntlet_contigs_2016-11-25_ei_version.fasta.gz $input

#Copying the contaminating P. syringae genome assembly over
#PacBio assembly
cp /home/hulinm/Richard_Harrison_EMR.RH.ENQ-933.01rev01_S7_S8_S9R1/circlator/5244.fasta $input
mv 5244.fasta Pacbio_5244.fasta
cp /home/hulinm/minion/psm5244/5244.fasta $input
mv 5244.fasta minion_5244.fasta

#Assembly stats
$scripts/assemblathon_stats.pl minion_5244.fasta >minion_5244.stats
$scripts/assemblathon_stats.pl Pacbio_5244.fasta >Pacbio_5244.stats
#7 vs 5 contigs in nanopore and PacBio assemblies. Similar total size >6 Mbp
#Will use the slightly better PacBio assembly as reference.
$scripts/assemblathon_stats.pl redgauntlet_contigs_2016-11-25_ei_version.fasta \
>redgauntlet.stats

#Copying nanopore reads for temporary storage
#downloaded - after basecalling with Metrichor
cp -r $input_files/downloaded $input

#"Pass" reads
sh $scripts/poretools.sh ./downloaded/pass

#"Fail" reads
sh $scripts/poretools.sh ./downloaded/fail

#Detection of possible overrepresented adapter sequences with FASTQC
sh $scripts/run_fastqc.sh --nogroup ./downloaded/pass.fastq
sh $scripts/run_fastqc.sh --nogroup ./downloaded/fail.fastq
#Found Illumina Single End PCR Primer 1 to be overrepresented - strawberry reads

#Carry out adapter trimming
fastq_mcf=/home/armita/all_idris/idris/prog/ea-utils.1.1.2-537/fastq-mcf
ILLUMINA_ADAPTERS=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
input_file=downloaded/pass.fastq
output=downloaded/pass_noadap.fastq
$fastq_mcf $ILLUMINA_ADAPTERS $input/$input_file -o $input/$output -0 -t 0.01 -p 5
#Result:
#Files: 1
#Total reads: 1314
#Too short after clip: 0
#Clipped 'end' reads: Count: 675, Mean: 885.97, Sd: 1336.50
input_file=downloaded/fail.fastq
output=downloaded/fail_noadap.fastq
$fastq_mcf $ILLUMINA_ADAPTERS $input/$input_file -o $input/$output -0 -t 0.01 -p 5

#However, first going to map all the reads to the P. syringae PacBio genome to
#eliminate contaminating reads. Takes about 10 hours at this rate.
sh $scripts/nanook.sh $input Pacbio_5244.fasta

#Next, align the reads to the Fragaria ananassa genome
sh $scripts/nanook.sh $input/strawberry redgauntlet_contigs_2016-11-25_ei_version.fasta

#Afterwards, six-frame translate the contigs and carry out NBS domain detection with
#HMMMER and NLR Parser.

##NLR Parser:
mkdir downloaded/nlrparser
cp downloaded/pass.fasta downloaded/nlrparser
cp downloaded/fail.fasta downloaded/nlrparser
sh /home/sobczm/bin/popgen/renseq/sub_nlrparser.sh $input/downloaded/nlrparser/pass.fasta
sort -k 1 fail_nlr.tsv >fail_nlr_sorted.tsv
#Check the number of unique reads with annotations
cut -f1 fail_nlr_sorted.tsv | cut -d"_" -f1 | sort | uniq | wc
#20 reads (unique)
sh /home/sobczm/bin/popgen/renseq/sub_nlrparser.sh $input/downloaded/nlrparser/fail.fasta
sort -k 1 pass_nlr.tsv >pass_nlr_sorted.tsv
#Check the number of unique reads with annotations
cut -f1 pass_nlr_sorted.tsv | cut -d"_" -f1 | sort | uniq | wc
#58 reads (unique)

##HMMMER
java -jar /home/sobczm/bin/popgen/renseq/Translate6Frame.jar \
-i downloaded/pass.fasta -o downloaded/pass_prot.fasta
java -jar /home/sobczm/bin/popgen/renseq/Translate6Frame.jar \
-i downloaded/fail.fasta -o downloaded/fail_prot.fasta

mkdir fail_hmmer
cp downloaded/fail_prot.fasta fail_hmmer
mkdir pass_hmmer
cp downloaded/pass_prot.fasta pass_hmmer

cd pass_hmmer
for file_i in pass_prot.fasta
do
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $file_i>temp && mv temp $file_i
awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%100==0){file=sprintf("myseq%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' <$file_i
for file in myseq*.fa
do
    Jobs=$(qstat | grep 'sub_hmmsca' | wc -l)
    while [ $Jobs -gt 100 ]
    do
        sleep 10
        printf "."
        Jobs=$(qstat | grep 'sub_hmmsca' | wc -l)
    done
qsub /home/sobczm/bin/popgen/renseq/sub_hmmscan.sh $file
done
done

input2=/home/sobczm/popgen/renseq/input/transcriptomes/plant_rgenes/common_onion
rgenes=/home/sobczm/bin/plant_rgenes

#Parse output
for a in *.out; do cat $a >> pass_hmmer.out; done
#Prepare input for domain analysis
mkdir domains
cp pass_prot.fasta pass_hmmer.out domains
cd domains
cp $input2/test2/db_descriptions.txt domains
mv pass_prot.fasta onion_167_TAIR10.protein.fa
mv pass_hmmer.out onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.out
#Carry out domain parsing
perl $rgenes/processing_scripts/K-parse_Pfam_domains_v3.1.pl --pfam onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.out \
--evalue 0.0001 --out onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.parsed.verbose --verbose T
perl $rgenes/processing_scripts/K-parse_Pfam_domains_NLR-fusions-v2.3.pl \
--indir ./ --evalue 0.0001 --outdir ./ --db_description db_descriptions.txt
#Sort the output
for a in onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.parsed.verbose.NLR*
do
sort -k 1 $a >${a%.*}_sorted.txt
done

cd fail_hmmer
for file_i in fail_prot.fasta
do
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $file_i>temp && mv temp $file_i
awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%100==0){file=sprintf("myseq%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' <$file_i
for file in myseq*.fa
do
    Jobs=$(qstat | grep 'sub_hmmsca' | wc -l)
    while [ $Jobs -gt 100 ]
    do
        sleep 10
        printf "."
        Jobs=$(qstat | grep 'sub_hmmsca' | wc -l)
    done
qsub /home/sobczm/bin/popgen/renseq/sub_hmmscan.sh $file
done
done

input2=/home/sobczm/popgen/renseq/input/transcriptomes/plant_rgenes/common_onion
rgenes=/home/sobczm/bin/plant_rgenes

#Parse output
for a in *.out; do cat $a >> fail_hmmer.out; done
#Prepare input for domain analysis
mkdir domains
cp fail_prot.fasta fail_hmmer.out domains
cd domains
cp $input2/test2/db_descriptions.txt domains
mv fail_prot.fasta onion_167_TAIR10.protein.fa
mv fail_hmmer.out onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.out
#Carry out domain parsing
perl $rgenes/processing_scripts/K-parse_Pfam_domains_v3.1.pl --pfam onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.out \
--evalue 0.0001 --out onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.parsed.verbose --verbose T
perl $rgenes/processing_scripts/K-parse_Pfam_domains_NLR-fusions-v2.3.pl \
--indir ./ --evalue 0.0001 --outdir ./ --db_description db_descriptions.txt
#Sort the output
for a in onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.parsed.verbose.NLR*
do
sort -k 1 $a >${a%.*}_sorted.txt
done
