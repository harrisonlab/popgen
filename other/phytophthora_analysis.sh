#!/bin/bash
scripts=/home/sobczm/bin/popgen/other
input=/home/sobczm/popgen/other/phytophthora

#Get OrthoMCL results into a more readable tabular format
python $scripts/rearrange_orthomcl_results_with_rubi.py All_Strains_plus_rubi_orthogroups.txt
python $scripts/rearrange_orthomcl_results_without_rubi.py All_Strains_orthogroups.txt
#How many orthogroups present in all strains but BC16?
Orthogroups=All_Strains_orthogroups.txt
cat $Orthogroups | grep 'A4|' | grep 'Bc1|' | grep 'Bc23|' | grep 'Nov27|' | grep 'Nov5|' | grep 'Nov71|' | grep 'Nov77|' | grep 'Nov9|' | grep 'ONT3|' | grep 'SCRP245_v2|' | grep -v 'Bc16|' >no_Bc16_groups.txt
#There appears to be 230 such orthogroups, out of 19952.

scripts=/home/sobczm/bin/popgen/summary_stats
input=/home/sobczm/popgen/other/phytophthora
#Outgroup-based tests for selection in P. fragariae
#1
#P. rubi - used as outgroup (incorrectly, as derived from P. fragariae) in order
#to test the ancestral allele anonotation procedure using VCF polymorphism data stemming from alignment of genome outgroup reads to the focal species genome.
#2
#P. sojae, P. ramorum (P. ramorum is quite distantly related to the other two species)
#in order to test the ancestral allele anonotation procedure using Mauve-based whole-genome alignment.
#3
#Lastly, the results of ancenstral allele annotation using VCF-based annotation
#will be compared to the whole genome alignment annotation.

#Obtain masked P. fragariae, P. sojae, P. ramorum sequences
#P. fragariae
ProjDir=/home/groups/harrisonlab/project_files/phytophthora_fragariae
cd $input/genomes
cp $ProjDir/repeat_masked/P.fragariae/Bc16/filtered_contigs_repmask/95m_contigs_hardmasked.fa ./
#P. sojae
wget ftp://ftp.ensemblgenomes.org/pub/protists/release-34/fasta/phytophthora_sojae/dna/Phytophthora_sojae.P_sojae_V3_0.dna_rm.toplevel.fa.gz
#P. ramorum
wget ftp://ftp.ensemblgenomes.org/pub/protists/release-34/fasta/phytophthora_ramorum/dna/Phytophthora_ramorum.ASM14973v1.dna_rm.toplevel.fa.gz

##Use unannotated VCF files for analysis below to avoid conflict with the snpEff annotation
##1
python $scripts/annotate_vcf_aa.py $input/SNP_calling/95m_contigs_unmasked_filtered.vcf 2 SCRP249,,SCRP324,,SCRP333
##2
#mauveAligner
qsub $scripts/run_mauve_aligner.sh $input/genomes/aligner \
"95m_contigs_hardmasked.fa 95m_contigs_hardmasked.fa.sml Phytophthora_sojae.P_sojae_V3_0.dna_rm.toplevel.fa Phytophthora_sojae.P_sojae_V3_0.dna_rm.toplevel.fa.sml Phytophthora_ramorum.ASM14973v1.dna_rm.toplevel.fa Phytophthora_ramorum.ASM14973v1.dna_rm.toplevel.fa.sml"
#A C float error has occurred so running on a worker node instead
$mauve/mauveAligner --weight=300 --output=output.mauve --output-alignment=aligned_genomes.xmfa \
"95m_contigs_hardmasked.fa 95m_contigs_hardmasked.fa.sml Phytophthora_sojae.P_sojae_V3_0.dna_rm.toplevel.fa Phytophthora_sojae.P_sojae_V3_0.dna_rm.toplevel.fa.sml Phytophthora_ramorum.ASM14973v1.dna_rm.toplevel.fa Phytophthora_ramorum.ASM14973v1.dna_rm.toplevel.fa.sml"
#progressiveMauve
qsub $scripts/run_progressive_mauve.sh $input/genomes/progressive "95m_contigs_hardmasked.fa Phytophthora_sojae.P_sojae_V3_0.dna_rm.toplevel.fa Phytophthora_ramorum.ASM14973v1.dna_rm.toplevel.fa"
perl /home/sobczm/bin/popoolation_1.2.2/mauve-parser.pl --ref $input/genomes/95m_contigs_hardmasked.fa \
--input $input/genomes/progressive/aligned_genomes.xmfa --output $input/genomes/progressive/mel-guided-alignment.txt
#Correct output obtained, and as progressivemauve over best alignment overall, going to retain this part of the analysis.
#Option 'Y' specifies to print fake genotype into the VCF file encoding the identified ancestral alleles.
python $scripts/annotate_gen_aa.py $input/genomes/progressive/mel-guided-alignment.txt \
$input/SNP_calling/95m_contigs_unmasked_filtered.vcf 2 Y
#Carry out the analysis above without printing fake genotypes.
python $scripts/annotate_gen_aa.py $input/genomes/progressive/mel-guided-alignment.txt \
$input/SNP_calling/95m_contigs_unmasked_filtered.vcf 2 N
#3
#Compare the results of ancestral allele annotation obtained using VCF and genome alignment
#and print both AA field and fake genotype with the ancestral allele:
python $scripts/compare_outgroup_results.py $input/SNP_calling/95m_contigs_unmasked_filtered_aa.vcf \
$input/SNP_calling/95m_contigs_unmasked_filtered_gen_aa.vcf 2 N
