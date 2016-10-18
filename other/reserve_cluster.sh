#$ -S /bin/bash
#$ -l h_vmem=1G
#$ -l mem_free=0.4G
#$ -l virtual_free=0.4G
#$ -l h_rt=9999:00:00
###$ -l h=blacklace01.blacklace
#$ -t 16

# hold a node for a qlogin session

#myqsub.sh ~/misc/sample_script.sh

while true
do
    sleep 10
done

```
The main resources you will be asking for are CPU(slots) and memory.

The memory can be requested with h_vmem, this is the maximum amount of memory a slot can take, if it exceeds this amount your job will be automatically killed. Another memory related option is mem_free

Now, we can also ask for more CPU using -pe option. This option stands for "parallel environment" and each cluster sets up a set of parallel environments for parallel processing. The syntax for the option is as follows:-pe [environment name] [number of slots].
```
