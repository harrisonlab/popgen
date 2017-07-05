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
#DO NOT CHANGE THE $start VALUE BELOW:
grep "#CHROM" $input_file | head -1 | awk '{for(start=10;start<=NF;++start)print $i " " $i "_pop"}' >"${input_file%.vcf}.lst"
#Copy the configuration file and change the path to the population definition file.
#For haploid organisms:
config=vcf_to_structure_haploid_pop.spid
#For diploid organisms:
#config=vcf_to_structure_diploid_pop.spid
cp $pgdspid/$config ./
dir=$PWD
list_file=$(echo "${input_file%.vcf}.lst")
sed -i 's,^\(VCF_PARSER_POP_FILE_QUESTION=\).*,\1'"$list_file"',' vcf_to_structure_haploid_pop.spid

#Also, create a label reference file to be used in the final step by distruct to label indidviduals in the output
#DO NOT CHANGE THE $start VALUE BELOW:
names="${input_file%.vcf}.label"
grep "#CHROM" $input_file | head -1 | awk '{for(start=10;start<=NF;++start)print $i }' >temp
nl temp | sed 's/^ *//' | sed 's/\t/ /g' >$names
outfile="${input_file%.vcf}.struc"
#Execute VCF to .struc (input format for STRUCTURE) conversion
java -jar $pgdspid/PGDSpider2-cli.jar -inputfile $input_file \
-inputformat VCF -outputfile $outfile -outputformat STRUCTURE -spid vcf_to_structure_haploid_pop.spid
dos2unix $outfile