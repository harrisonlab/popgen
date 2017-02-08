#!/bin/bash
input=/home/sobczm/popgen/snp/snp_calling/multithreaded
scripts=/home/sobczm/bin/popgen/snp
pgdspid=/home/sobczm/bin/PGDSpider_2.1.0.3

cd $input
#Downsample SNPs for Structure analysis as enough information in 10% of the loci
#(and more not likely informative because of linkage). In certain cases, when small number of markers detected,
#this step is unnecessary and all can be retained.
/home/sobczm/bin/vcflib/bin/vcfrandomsample \
--rate 0.1 Fus2_canu_contigs_unmasked_filtered.vcf > Fus2_canu_contigs_unmasked_filtered_subsampled.vcf
#Prepare STRUCTURE input (PGDSpider does not work when wrapped up in a bash script, grrr)
#haploid (for diploid change the conversion script to vcf_to_structure_diploid.spid)
#!!!! Need to change the path to file with population definitions !!!
input_file=$input/Fus2_canu_contigs_unmasked_filtered_subsampled.vcf
#Prepare population definition file. Each individual = new population
grep "#CHROM" $input_file | head -1 | awk '{for(i=10;i<=NF;++i)print $i " " $i "_pop"}' >"${input_file%.vcf}.lst"
#Copy the configuration file and change the path to the population definition file.
#For haploid organisms:
config=vcf_to_structure_haploid_pop.spid
#For diploid organisms:
#config=vcf_to_structure_diploid_pop.spid
cp $pgdspid/$config ./
dir=$PWD
list_file=$(echo "${input_file%.vcf}.lst")
sed -i 's,^\(VCF_PARSER_POP_FILE_QUESTION=\).*,\1'"$dir/$list_file"',' vcf_to_structure_haploid_pop.spid

#Also, create a label reference file to be used in the final step by distruct to label indidviduals in the output
names="${input_file%.vcf}.label"
grep "#CHROM" $input_file | head -1 | awk '{for(i=10;i<=NF;++i)print $i }' >temp
nl temp | sed 's/^ *//' | sed 's/\t/ /g' >$names
outfile="${input_file%.vcf}.struc"
#Execute VCF to .struc (input format for STRUCTURE) conversion
java -jar $pgdspid/PGDSpider2-cli.jar -inputfile $input_file \
-inputformat VCF -outputfile $outfile -outputformat STRUCTURE -spid vcf_to_structure_haploid_pop.spid
dos2unix $outfile

#Run STRUCTURE analysis to test for the presence of
#K (population clusters) 1 to 11, with 5 replicates for each K run consecutively. 
### Value of K needs to be changed in each analysis, depends on the number of likely clusters observed among individuals to be evaluated
#Run replicate STRUCTURE runs, with K from 1 to 11
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
# structureHarvester - summarise the results
harvester=/home/sobczm/bin/structureHarvester/structureHarvester.py
$harvester --dir=$input/structureHarvester --out=$input/structureHarvester --evanno --clumpp
# CLUMPP - permute the results
cd structureHarvester
clumpp=/home/sobczm/bin/CLUMPP_Linux64.1.1.2
cp $clumpp/paramfile_ind ./
mv paramfile_ind paramfile
#Options fed to CLUMPP
#-i: indfile from StructureHarvester output
#-p: popfile from StructureHarvester output
#-o: output Q matrix for distruct input
#-k: K value (number of clusters tested)

###!!! Options to be changed in each analysis manually
#c: number of individuals (change according to STRUCTURE mainparam file)
#r: number of replicate runs
#s: minimum number of population clusters (K) tested
#f: maximum number of population clusters (K) tested
c=11
r=5
s=1
f=11
for i in {1..11} #input range of K values tested
do
$clumpp/CLUMPP -i K$i.indfile -p K$i.popfile -o K$i.indivq -k $i -c $c -r $r
done
cp $clumpp/paramfile_pop ./
mv paramfile_pop paramfile
for i in {$s..$f} #input range of K values tested
do
$clumpp/CLUMPP -i K$i.indfile -p K$i.popfile -o K$i.popq -k $i -c $c -r $r
done
#Key options in the paramfile
# DISTRUCT to visualise the results
###!!!! Options to be changed in each analysis manually
#-M number of populations assigned in the Structure input file
#-N number of individuals
m=11
n=11
#-K K value
#-p input file (population q's)
#-i input file (individual q's)
#-a input file (labels atop figure)
#-b input file (labels below figure)
#-o output file
distruct=/home/sobczm/bin/distruct1.1
cp $distruct/drawparams ./
for i in {$s..$f} #input range of K values tested
do
$distruct/distructLinux1.1 -i K$i.indivq -p K$i.popq -a $input/$names -o K$i.ps -k $i -M $m -N $n -K $i
done

#Output is a number of PostScript files showing the average proportion of each
#individual's genome belonging to a given cluster and allowing some ancestry inference based on
#the most likely true number of population clusters as summarised by StructureHarvester.
