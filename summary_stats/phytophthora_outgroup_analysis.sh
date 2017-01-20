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
##McDonald-Kreitman test (a couple of genotypes from the outgroup species required) calculated by PopGenome
##As an example, generate FASTA input using VCF created in A) 1
mkdir -p $input/mkt
cd $input/mkt
ref_genome=/home/sobczm/popgen/other/phytophthora/genomes/95m_contigs_hardmasked.fa
vcf_file=$input/SNP_calling/95m_contigs_unmasked_filtered_vcf_aa.vcf
python $scripts/vcf_to_fasta.py $vcf_file $ref_genome 2

#Prepare Popgenome input
function Popgenome {
mkdir contigs && mv *.fasta ./contigs
cd contigs
for f in *.fasta
do
folder=${f%.fasta}
mkdir $folder
mv $f $folder
done
#Gff files
cd ..
gff=/home/groups/harrisonlab/project_files/phytophthora_fragariae/gene_pred/codingquary/P.fragariae/Bc16/final/final_genes_appended.gff3
$scripts/split_gff_contig.sh $gff
mkdir gff && mv *.gff ./gff
#Check for orphan contigs with no matching gff file, which need to be removed prior to the run.
for a in $PWD/contigs/*/*.fasta
do
filename=$(basename "$a")
expected_gff="$PWD/gff/${filename%.fa*}.gff"
if [ ! -f "$expected_gff" ];
then
   rm -rf $(dirname $a)
fi
done
}
Popgenome

#Requires custom adjustment of the R script called below to include the samples being analysed.
qsub $scripts/sub_calculate_mkt.sh

################ Outgroup-based tests for selection 
##Fay & Wu's H (at least one outgroup genotype needed) calculated by PopGenome
##As an example, generate FASTA input using VCF created in A) 2
mkdir -p $input/faywuh
cd $input/faywuh
ref_genome=/home/sobczm/popgen/other/phytophthora/genomes/95m_contigs_hardmasked.fa
vcf_file=$input/SNP_calling/95m_contigs_unmasked_filtered_gen_aa.vcf
python $scripts/vcf_to_fasta.py $vcf_file $ref_genome 2
##Prepare Popgenome input
Popgenome
#Requires custom adjustment of the R script called below to include the samples being analysed.
qsub $scripts/sub_calculate_faywu.sh

