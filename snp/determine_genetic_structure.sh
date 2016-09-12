#!/bin/bash
input=/home/sobczm/popgen/snp/snp_calling/multithreaded
scripts=/home/sobczm/bin/popgen/snp

#Only retain biallelic high-quality SNPS with no missing data for genetic analyses.
cd $input
$scripts/vcf_parser_haploid.py --i Fus2_canu_contigs_unmasked.vcf
#General VCF stats (remember that vcftools needs to have the PERL library exported)
perl /home/sobczm/bin/vcftools/bin/vcf-stats \
Fus2_canu_contigs_unmasked.vcf >Fus2_canu_contigs_unmasked.stat
perl /home/sobczm/bin/vcftools/bin/vcf-stats \
Fus2_canu_contigs_unmasked_filtered.vcf >Fus2_canu_contigs_unmasked_filtered.stat
#Calculate the index for percentage of shared SNP alleles between the individs.
$scripts/similarity_percentage.py Fus2_canu_contigs_unmasked_filtered.vcf
#Visualise the output as heatmap and clustering dendrogram
Rscript --vanilla $scripts/distance_matrix.R Fus2_canu_contigs_unmasked_filtered_distance.log
#DAPC and AMOVA analysis
Rscript --vanilla $popgen/snp/amova_dapc.R
#Downsample SNPs for AMOVA analysis as enough information in 10% of the loci
#(and more not informative because of linkagee)
/home/sobczm/bin/vcflib/bin/vcfrandomsample \
--rate 0.1 Fus2_canu_contigs_unmasked_filtered.vcf > Fus2_canu_contigs_unmasked_filtered_subsampled.vcf
#Run STRUCTURE analysis for k 1 to 9, with 3 replicates for each k run consecutively.

#Prepare STRUCTURE input (PGDSpider does not work when wrapped up in a bash script, grrr)
#haploid (for diploid change the conversion script to vcf_to_structure_diploid.spid)
#!!!! Need to change the path to file with population definitions
pgdspid=/home/sobczm/bin/PGDSpider_2.1.0.3
input_file=Fus2_canu_contigs_unmasked_filtered_subsampled.vcf
filename=$(basename "$input_file")
outfile="${filename%.*}.struc"
java -jar $pgdspid/PGDSpider2-cli.jar -inputfile $input_file \
-inputformat VCF -outputfile $outfile -outputformat STRUCTURE -spid $pgdspid/vcf_to_structure_haploid_pop.spid
dos2unix $outfile

#Run replicate STRUCTURE runs, with K from 1 to 9
qsub $scripts/execute_structure.sh $input/Fus2_canu_contigs_unmasked_filtered_subsampled.struc 1 1 5
qsub $scripts/execute_structure.sh $input/Fus2_canu_contigs_unmasked_filtered_subsampled.struc 1 2 5
qsub $scripts/execute_structure.sh $input/Fus2_canu_contigs_unmasked_filtered_subsampled.struc 1 3 5
qsub $scripts/execute_structure.sh $input/Fus2_canu_contigs_unmasked_filtered_subsampled.struc 1 4 5
qsub $scripts/execute_structure.sh $input/Fus2_canu_contigs_unmasked_filtered_subsampled.struc 1 5 5
qsub $scripts/execute_structure.sh $input/Fus2_canu_contigs_unmasked_filtered_subsampled.struc 1 6 5
qsub $scripts/execute_structure.sh $input/Fus2_canu_contigs_unmasked_filtered_subsampled.struc 1 7 5
qsub $scripts/execute_structure.sh $input/Fus2_canu_contigs_unmasked_filtered_subsampled.struc 1 8 5
qsub $scripts/execute_structure.sh $input/Fus2_canu_contigs_unmasked_filtered_subsampled.struc 1 9 5
qsub $scripts/execute_structure.sh $input/Fus2_canu_contigs_unmasked_filtered_subsampled.struc 1 10 5
qsub $scripts/execute_structure.sh $input/Fus2_canu_contigs_unmasked_filtered_subsampled.struc 1 11 5

#Analyze STRUCTURE output
# Generate a folder containing all STRUCTURE output files for all K analyzed
mkdir structureHarvester
for d in $PWD/*
do
mv $d/*_f $PWD/structureHarvester
done
# structureHarvester
harvester=/home/sobczm/bin/structureHarvester/structureHarvester.py
$harvester --dir=$input/structureHarvester --out=$input/structureHarvester --evanno --clumpp
# CLUMPP
cd structureHarvester
clumpp=/home/sobczm/bin/CLUMPP_Linux64.1.1.2
cp $clumpp/paramfile_ind ./
mv paramfile_ind paramfile
#-i: indfile from StructureHarvester output
#-p: popfile from StructureHarvester output
#-o: output Q matrix for distruct input
#-k: K value
#-c: number of individuals (change according to STRUCTURE mainparam file)
#-r: number of replicate runs
for i in {1..3} #input range of K values tested
do
$clumpp/CLUMPP -i K$i.indfile -p K$i.popfile -o K$i.indivq -k $i -c 11 -r 2
done
cp $clumpp/paramfile_pop ./
mv paramfile_pop paramfile
for i in {1..3} #input range of K values tested
do
$clumpp/CLUMPP -i K$i.indfile -p K$i.popfile -o K$i.popq -k $i -c 11 -r 5
done
# DISTRUCT to visualise the results
#-K K value
#-M number of populations assigned in the Structure input file
#-N number of individuals
#-p input file (population q's)
#-i input file (individual q's)
#-a input file (labels atop figure)
#-b input file (labels below figure)
#-o output file
distruct=/home/sobczm/bin/distruct1.1
cp $distruct/drawparams ./
for i in {1..3} #input range of K values tested
do
$distruct/distructLinux1.1 -i K$i.indivq -p K$i.popq \
-a $input/fusarium.names -o K$i.ps -k $i -M 11 -N 11 -K $i
done
