#$ -S /bin/bash
#$ -cwd
#$ -pe smp 16
#$ -l h_vmem=1G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace12.blacklace

### Setting variables
dir=$1
orthofinder=/home/sobczm/bin/OrthoFinder-1.0.7/orthofinder
anaconda=/home/sobczm/bin/anaconda2/bin/python
cpath=$PWD
temp_dir="$TMPDIR"
new_name="run_$(date +%F_%R)"

### Prep
mkdir -p $temp_dir
cp -r $cpath/$dir $temp_dir
cd $temp_dir

### Execute
$anaconda $orthofinder/orthofinder.py -f $dir -t 16
mv $dir/Results_* $dir/$new_name
$anaconda $orthofinder/trees_from_MSA.py $dir/$new_name -t 16

### Cleanup
cp -r $temp_dir $cpath
rm -rf $temp_dir
