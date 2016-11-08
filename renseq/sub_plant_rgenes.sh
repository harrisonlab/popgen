#!/bin/bash


pfam_scan.pl -e_seq 1 -e_dom 1 -as -outfile test.out -cpu 8 -fasta onion_test_protein.fa \
-dir /home/sobczm/bin/hmmer-3.1b2/pfam30
perl $rgenes/processing_scripts/K-parse_Pfam_domains_v3.1.pl --pfam onion_test.out \
--evalue 0.0001 --out onion_test.parsed --verbose T
perl $rgenes/processing_scripts/K-parse_Pfam_domains_NLR-fusions-v2.2.pl \
--indir $input/test --evalue 0.0001 --outdir $input/test --db_description $rgenes/processing_scripts/onion_ref.txt
