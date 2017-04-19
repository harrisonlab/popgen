#!/bin/bash
input=/home/sobczm/popgen/clock/DNA_genomes
scripts=/home/sobczm/bin/popgen/clock/motif_discovery

#Get V. dahliae JR2 genome, extract promoters at 1 kbp, 2 kbp and 3 kbp using extract_promoter.sh.
cd $input/JR2
wget  ftp://ftp.ensemblgenomes.org/pub/fungi/release-34/fasta/verticillium_dahliaejr2/dna/Verticillium_dahliaejr2.GCA_000400815.2.dna_rm.toplevel.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-34/gff3/verticillium_dahliaejr2/Verticillium_dahliaejr2.GCA_000400815.2.34.gff3.gz

#Reverse complement the promoter sequences to get the orientation - START codon.
for a in *promoters*.fasta
do
/home/armita/prog/emboss/EMBOSS-4.0.0/bin/revseq $a ${a%.fasta}_revcom.fasta
done

#Look for the He (2005) motif with max. 30 bp gap between the two parts of the motif in the Vd JR2 genome.
cd $input/JR2
#For each promoter size (1000, 2000, 3000 bp)
python $scripts/find_gapped_motif.py Verticillium_dahliaejr2.GCA_000400815.2.dna_rm.toplevel_promoters_1000_revcom.fasta GATAC CGATA 0,30  Vd_He2005_A
python $scripts/find_gapped_motif.py Verticillium_dahliaejr2.GCA_000400815.2.dna_rm.toplevel_promoters_1000_revcom.fasta GATTC CGATT 0,30  Vd_He2005_T
python $scripts/find_gapped_motif.py Verticillium_dahliaejr2.GCA_000400815.2.dna_rm.toplevel_promoters_1000_revcom.fasta GATGC CGATG 0,30  Vd_He2005_G
python $scripts/find_gapped_motif.py Verticillium_dahliaejr2.GCA_000400815.2.dna_rm.toplevel_promoters_1000_revcom.fasta GATCC CGATC 0,30  Vd_He2005_C

python $scripts/total_no_motifs.py

#Now, in 2000bp of Neurospora promoters.
cd $input/promoters/extended/frq
python $scripts/find_gapped_motif.py Neurospora_crassa.NC12.dna_rm.toplevel_promoters_2000_filtered_revcom.fasta GATAC CGATA 0,30  Nc_He2005_A
python $scripts/find_gapped_motif.py Neurospora_crassa.NC12.dna_rm.toplevel_promoters_2000_filtered_revcom.fasta GATTC CGATT 0,30  Nc_He2005_T
python $scripts/find_gapped_motif.py Neurospora_crassa.NC12.dna_rm.toplevel_promoters_2000_filtered_revcom.fasta GATGC CGATG 0,30  Nc_He2005_G
python $scripts/find_gapped_motif.py Neurospora_crassa.NC12.dna_rm.toplevel_promoters_2000_filtered_revcom.fasta GATCC CGATC 0,30  Nc_He2005_C

#And get overall summary for all the different riffs on the motif
python $scripts/total_no_motifs.py