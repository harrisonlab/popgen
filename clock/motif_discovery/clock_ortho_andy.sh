#!/bin/bash
input=/home/sobczm/popgen/clock/pep_genomes
scripts=/home/armita/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
orthomcl=/home/armita/prog/ortho-mcl/orthomclSoftware-v2.0.9/bin

#Identification of orthologs amongst selected genomes using
#Andy's orthoMCL pipeline
#https://github.com/harrisonlab/fusarium/blob/master/pathogen/orthology/F.oxysporum_fsp.cepae_pathogen_vs_non-pathogen_orthology.md

#Set up a directory for the run, and copy all input fasta files there.
#Create the dir tree for analysis.
cd $input
mkdir orthoMCL
cp *.fa ./orthoMCL

mkdir orthoMCL/formatted
mkdir orthoMCL/goodProteins
mkdir orthoMCL/badProteins

cd orthoMCL
#Prepare orthoMCL input
Taxon_code=BoC
Fasta_file=Botrytis_cinerea.pep.fa
Id_field=1
$orthomcl/orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $input/orthoMCL/formatted/"$Taxon_code".fasta

Taxon_code=Fus2
Fasta_file=Fus2.pep.fa
Id_field=1
$orthomcl/orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $input/orthoMCL/formatted/"$Taxon_code".fasta

Taxon_code=FuG
Fasta_file=Fusarium_graminearum.pep.fa
Id_field=1
$orthomcl/orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $input/orthoMCL/formatted/"$Taxon_code".fasta

Taxon_code=FuS
Fasta_file=Fusarium_solani.pep.fa
Id_field=1
$orthomcl/orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $input/orthoMCL/formatted/"$Taxon_code".fasta

Taxon_code=MaO
Fasta_file=Magnaporthe_oryzae.pep.fa
Id_field=1
$orthomcl/orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $input/orthoMCL/formatted/"$Taxon_code".fasta

Taxon_code=NeD
Fasta_file=Neonectria_ditissima.pep.fa
Id_field=1
$orthomcl/orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $input/orthoMCL/formatted/"$Taxon_code".fasta

Taxon_code=NsC
Fasta_file=Neurospora_crassa.pep.fa
Id_field=1
$orthomcl/orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $input/orthoMCL/formatted/"$Taxon_code".fasta

Taxon_code=PoA
Fasta_file=Podospora_anserina.pep.fa
Id_field=1
$orthomcl/orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $input/orthoMCL/formatted/"$Taxon_code".fasta

Taxon_code=SoM
Fasta_file=Sordaria_macrospora.pep.fa
Id_field=1
$orthomcl/orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $input/orthoMCL/formatted/"$Taxon_code".fasta

Taxon_code=TrR
Fasta_file=Trichoderma_reesei.pep.fa
Id_field=1
$orthomcl/orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $input/orthoMCL/formatted/"$Taxon_code".fasta

Taxon_code=VeA
Fasta_file=Verticillium_alfalfae.pep.fa
Id_field=1
$orthomcl/orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $input/orthoMCL/formatted/"$Taxon_code".fasta

Taxon_code=VeD
Fasta_file=Verticillium_dahliae.pep.fa
Id_field=1
$orthomcl/orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $input/orthoMCL/formatted/"$Taxon_code".fasta

#Filter proteins into good and poor sets.
Input_dir=$input/orthoMCL/formatted/
Min_length=10
Max_percent_stops=20
Good_proteins_file=$input/orthoMCL/goodProteins/goodProteins.fasta
Poor_proteins_file=$input/orthoMCL/badProteins/poorProteins.fasta

$orthomcl/orthomclFilterFasta $Input_dir $Min_length $Max_percent_stops $Good_proteins_file $Poor_proteins_file

#Perform an all-vs-all blast of the proteins
mkdir $input/orthoMCL/blastall
BlastDB=$input/orthoMCL/blastall/ClockGenomes.db
makeblastdb -in $Good_proteins_file -dbtype prot -out $BlastDB
mkdir -p $input/orthoMCL/splitfiles
SplitDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
$SplitDir/splitfile_500.py --inp_fasta $Good_proteins_file --out_dir $input/orthoMCL/splitfiles --out_base goodProteins

ProgDir=/home/sobczm/bin/popgen/clock/motif_discovery
cd $input/orthoMCL/blastall
  for File in $input/orthoMCL/splitfiles/goodProteins*.fa; do
    Jobs=$(qstat | grep 'blast_500' | grep 'qw' | wc -l)
    while [ $Jobs -gt 1 ]; do
      sleep 3
      printf "."
      Jobs=$(qstat | grep 'blast_500' | grep 'qw' | wc -l)
    done
    printf "\n"
    echo $File
    BlastOut=$(echo $File | sed 's/.fa/.tab/g')
    qsub $ProgDir/blast_500.sh $BlastDB $File $BlastOut
  done

#Merge the all-vs-all blast results
MergeHits=clock_blast.tab
printf "" > $MergeHits
for Num in $(ls $input/orthoMCL/splitfiles/*.tab | rev | cut -f1 -d '_' | rev | sort -n); do
  File=$(ls $input/orthoMCL/splitfiles/*_$Num)
  cat $File
done > $MergeHits

#Perform ortholog identification
ProgDir=/home/sobczm/bin/popgen/clock/motif_discovery
  MergeHits=$input/orthoMCL/splitfiles/clock_blast.tab
  GoodProts=$input/orthoMCL/goodProteins/goodProteins.fasta
  qsub $ProgDir/sub_orthomcl.sh $MergeHits $GoodProts 1.5
