#!/bin/bash
scripts=/home/sobczm/bin/popgen/summary_stats
input=/home/sobczm/popgen/other/phytophthora

##### A)
################ Establish ancestral allele using genotype(s) from select outgroup species
#1
#P. rubi - used as outgroup (incorrectly, as derived from P. fragariae) in order
#to test the ancestral allele anonotation procedure using VCF polymorphism data 
# stemming from alignment of genome outgroup reads to the focal species genome.
#2
#Outgroups: P. sojae, P. ramorum (P. ramorum is quite distantly related 
#to the other two species) in order to test the ancestral allele anonotation procedure 
#using Mauve-based whole-genome alignment.
#3
#Lastly, the results of ancenstral allele annotation using VCF-based annotation (#1)
#will be compared to the whole genome alignment annotation (#2).

#Choice of the analysis to follow: 1, 2, 3 depends on the available resources and researcher preferences.
#Annotation with ancestral alleles can be used just to polarise the mutation status of SNPs of interest
#or can be used in the formal tests for selection (e.g. McDonald-Kreitman Test and Fay & Wu's H described below)

#Obtain masked P. fragariae, P. sojae, P. ramorum sequences
#P. fragariae
ProjDir=/home/groups/harrisonlab/project_files/phytophthora_fragariae
cd $input/genomes
cp $ProjDir/repeat_masked/P.fragariae/Bc16/filtered_contigs_repmask/95m_contigs_hardmasked.fa ./
#P. sojae
wget ftp://ftp.ensemblgenomes.org/pub/protists/release-34/fasta/phytophthora_sojae/dna/Phytophthora_sojae.P_sojae_V3_0.dna_rm.toplevel.fa.gz
#P. ramorum
wget ftp://ftp.ensemblgenomes.org/pub/protists/release-34/fasta/phytophthora_ramorum/dna/Phytophthora_ramorum.ASM14973v1.dna_rm.toplevel.fa.gz

##1
python $scripts/annotate_vcf_aa.py $input/SNP_calling/95m_contigs_unmasked_filtered.vcf 2 SCRP249,,SCRP324,,SCRP333

##2
#progressiveMauve
qsub $scripts/run_progressive_mauve.sh $input/genomes/progressive "95m_contigs_hardmasked.fa Phytophthora_sojae.P_sojae_V3_0.dna_rm.toplevel.fa Phytophthora_ramorum.ASM14973v1.dna_rm.toplevel.fa"
perl /home/sobczm/bin/popoolation_1.2.2/mauve-parser.pl --ref $input/genomes/95m_contigs_hardmasked.fa \
--input $input/genomes/progressive/aligned_genomes.xmfa --output $input/genomes/progressive/mel-guided-alignment.txt
#Option 'Y' specifies to print fake genotype into the VCF file encoding the identified ancestral alleles.
#Use this option when proceeding to use Popgenome in order to calculate outgroup-based
#statistics: Fay & Wu's H and McDonald-Kreitman test 
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

##### B)
################ Outgroup-based tests for selection 
##Fay & Wu's H (at least one outgroup genotype needed)
qsub $scripts/sub_calculate_faywu.sh

##McDonald-Kreitman test (a couple of genotypes from the outgroup species required)
qsub $scripts/sub_calculate_mkt.sh