#$ -S /bin/bash
#$ -cwd
#$ -pe smp 6
#$ -l h_vmem=1G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace12.blacklace

#### MAFFT alignment with the highest accuracy method L-INS-i (less than <200 sequences)
#### The scripts separately aligns all the FASTA sequences contained within each file in the directory.
#### OUTPUT: aligned FASTA file, with "aligned" suffix

mafft_path=/home/sobczm/bin/mafft-7.222/bin/mafft

for input in *.fasta
do
filename=$(basename "$input")
output="${filename%.*}_aligned.fasta"
$mafft_path --localpair --maxiterate 1000 $input >$output
done

#Convert to single line FASTA for easy parsing
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $output >temp && mv temp $output
