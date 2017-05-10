#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace11.blacklace|blacklace12.blacklace

#First argument: input STRUCTURE file
#Second argument: ploidy (value: 1|2)
#Third argument: K being tested (value: integer)
#Fourth argument: number of replicate runs (value: integer)
#Can change the number of burn-in and run replicates in the mainparams file
#Assumes that number of individuals = max number of populations (can change it below)
#Example usage: sh ./execute_structure.sh test.struc 1 1 3

input=$1
ploidy=$2
k=$3
n=$4
temp_dir="$TMPDIR"
cdir=$PWD

#Copies the structure folder to the current directory to allow parallel runs
struct=/home/sobczm/bin/structure

filename=$(basename "$input")
outfile="${filename%.*}"

### Prep
mkdir -p $temp_dir
cp -r $struct $temp_dir
mv $temp_dir/structure $temp_dir/structure_$k\_$n
cp $input $temp_dir/structure_$k\_$n

cd $temp_dir/structure_$k\_$n

# Substitute in the config file for the current run:
#input file name
sed -i 's,^\(#define INFILE \).*,\1'"$input"',' mainparams
#output
#sed -i 's,^\(#define OUTFILE \).*,\1'"$cdir/$outfile"',' mainparams
#ploidy
sed -i 's,^\(#define PLOIDY \).*,\1'"$ploidy"',' mainparams
#number of loci
a=`awk '{print $NF}' $input | head -1 | sed 's/SNP_//'`
sed -i 's,^\(#define NUMLOCI \).*,\1'"$a"',' mainparams
#number of inds
b=`wc -l <$input`
c=$(expr $b - 1)
sed -i 's,^\(#define NUMINDS \).*,\1'"$c"',' mainparams
#number of pops = number of inds
sed -i 's,^\(#define MAXPOPS \).*,\1'"$c"',' mainparams

for ((i=1;i<=$n;i++))
do
    $temp_dir/structure_$k\_$n/structure -K $k -o ${outfile}_k${k}_${i}
done

cp -r $temp_dir/structure_$k\_$n $cdir
