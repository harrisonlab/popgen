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

perl $rgenes/processing_scripts/K-parse_Pfam_domains_v3.1.pl --pfam $input/KIM/onion_kim_protein.out \
--evalue 0.0001 --out $input/KIM/onion_kim_protein.parsed --verbose T
perl $rgenes/processing_scripts/K-parse_Pfam_domains_v3.1.pl --pfam $input/RAJ/onion_raj_protein.out \
--evalue 0.0001 --out $input/RAJ/onion_raj_protein.parsed --verbose T
perl $rgenes/processing_scripts/K-parse_Pfam_domains_v3.1.pl --pfam $input/NZ/onion_nz_protein.out \
--evalue 0.0001 --out $input/NZ/onion_nz_protein.parsed --verbose T

perl $rgenes/processing_scripts/K-parse_Pfam_domains_NLR-fusions-v2.2.pl \
--indir $input/KIM --evalue 0.0001 --outdir $input/KIM --db_description $rgenes/processing_scripts/onion_ref.txt
perl $rgenes/processing_scripts/K-parse_Pfam_domains_NLR-fusions-v2.2.pl \
--indir $input/RAJ --evalue 0.0001 --outdir $input/RAJ --db_description $rgenes/processing_scripts/onion_ref.txt
perl $rgenes/processing_scripts/K-parse_Pfam_domains_NLR-fusions-v2.2.pl \
--indir $input/NZ --evalue 0.0001 --outdir $input/NZ --db_description $rgenes/processing_scripts/onion_ref.txt
