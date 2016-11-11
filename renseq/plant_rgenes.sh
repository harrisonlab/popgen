#!/bin/bash
input=/home/sobczm/popgen/renseq/input/transcriptomes/plant_rgenes/onion
scripts=/home/sobczm/bin/popgen/renseq
rgenes=/home/sobczm/bin/plant_rgenes

#Download the following files from the Pfam ftp site (Pfam v. 30.0)
#Pfam-A.hmm
#Pfam-A.hmm.dat
#active_site.dat

#Generate binary files for Pfam-A.hmm by running the following commands:
hmmpress /home/sobczm/bin/hmmer-3.1b2/pfam30/Pfam-A.hmm

#Prepare 6 frame translations of all contigs in the assemblies
java -jar $scripts/Translate6Frame.jar -i $input/KIM/GBRQ01_1_fsa_nt_combined_kim.fasta \
-o onion_kim_protein.fa
java -jar $scripts/Translate6Frame.jar -i $input/NZ/GBGJ01_1_fsa_nt_nz.fasta \
-o onion_nz_protein.fa
java -jar $scripts/Translate6Frame.jar -i $input/RAJ/GBJZ01_1_fsa_nt_raj.fasta \
-o onion_kim_protein.fa
java -jar $scripts/Translate6Frame.jar -i $input/NZ/GBGJ01_1_fsa_nt_nz.fasta \
-o onion_nz_protein.fa
java -jar $scripts/Translate6Frame.jar -i $input/RAJ/GBJZ01_1_fsa_nt_raj.fasta \
-o onion_raj_protein.fa

#Execute the pipeline using all 6-frame translations of the assemblies
#Substitute the run_pfam_scan.sh with a simplified command
pfam_scan.pl -e_seq 1 -e_dom 1 -as -outfile onion_kim_protein.out -cpu 8 -fasta $input/KIM/onion_kim_protein.fa \
-dir /home/sobczm/bin/hmmer-3.1b2/pfam30
pfam_scan.pl -e_seq 1 -e_dom 1 -as -outfile onion_raj_protein.out -cpu 8 -fasta $input/RAJ/onion_raj_protein.fa \
-dir /home/sobczm/bin/hmmer-3.1b2/pfam30
pfam_scan.pl -e_seq 1 -e_dom 1 -as -outfile onion_nz_protein.out -cpu 8 -fasta $input/NZ/onion_nz_protein.fa \
-dir /home/sobczm/bin/hmmer-3.1b2/pfam30

#As this taking too long, going to parallelise in a crude way by splitting
#into small files with 100 protein sequences and running 20 at a time.
cd $input/KIM
file=onion_kim_protein.fa

#!/bin/bash
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $file>temp && mv temp $file
awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%100==0){file=sprintf("myseq%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < $file
for file in myseq*.fa
do
    Jobs=$(qstat | grep 'sub_hmmsca' | wc -l)
    while [ $Jobs -gt 20 ]
    do
        sleep 10
        printf "."
        Jobs=$(qstat | grep 'sub_hmmsca' | wc -l)
    done
qsub $scripts/sub_hmmscan.sh $file
done

cd $input/RAJ
file=onion_raj_protein.fa
#execute script lines 40-53.

cd $input/NZ
file=onion_nz_protein.fa
#execute script lines 40-53.

#concatenate the outputs into one file
for a in $input/KIM/*.out; do cat $a >> onion_kim_protein.out; done
for a in $input/RAJ/*.out; do cat $a >> onion_raj_protein.out; done
for a in $input/NZ/*.out; do cat $a >> onion_nz_protein.out; done

#Prepare input for domain analysis

mkdir -p $input/domains/kim
mkdir -p $input/domains/raj
mkdir -p $input/domains/nz

cp $input/KIM/onion_kim_protein* $input/domains/kim
cp $input/test2/db_descriptions.txt $input/domains/kim
mv $input/domains/kim/onion_kim_protein.fa onion_167_TAIR10.protein.fa
mv $input/domains/kim/onion_kim_protein.out onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.out

cp $input/RAJ/onion_raj_protein* $input/domains/raj
cp $input/test2/db_descriptions.txt $input/domains/raj
mv $input/domains/raj/onion_raj_protein.fa onion_167_TAIR10.protein.fa
mv $input/domains/raj/onion_raj_protein.out onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.out

cp $input/NZ/onion_nz_protein* $input/domains/nz
cp $input/test2/db_descriptions.txt $input/domains/nz
mv $input/domains/nz/onion_nz_protein.fa onion_167_TAIR10.protein.fa
mv $input/domains/nz/onion_nz_protein.out onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.out


perl $rgenes/processing_scripts/K-parse_Pfam_domains_v3.1.pl --pfam $input/domains/kim/onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.out \
--evalue 0.0001 --out $input/domains/kim/onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.parsed.verbose --verbose T
perl $rgenes/processing_scripts/K-parse_Pfam_domains_v3.1.pl --pfam $input/domains/raj/onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.out \
--evalue 0.0001 --out $input/domains/raj/onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.parsed.verbose --verbose T
perl $rgenes/processing_scripts/K-parse_Pfam_domains_v3.1.pl --pfam $input/domains/nz/onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.out \
--evalue 0.0001 --out $input/domains/nz/onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.parsed.verbose --verbose T

perl $rgenes/processing_scripts/K-parse_Pfam_domains_NLR-fusions-v2.3.pl \
--indir $input/domains/kim --evalue 0.0001 --outdir $input/domains/kim --db_description $input/domains/kim/db_descriptions.txt
perl $rgenes/processing_scripts/K-parse_Pfam_domains_NLR-fusions-v2.3.pl \
--indir $input/domains/raj --evalue 0.0001 --outdir $input/domains/raj --db_description $input/domains/raj/db_descriptions.txt
perl $rgenes/processing_scripts/K-parse_Pfam_domains_NLR-fusions-v2.3.pl \
--indir $input/domains/nz --evalue 0.0001 --outdir $input/domains/nz --db_description $input/domains/nz/db_descriptions.txt
