#$ -S /bin/bash
#$ -l h_vmem=5G
#$ -l mem_free=5G
#$ -l virtual_free=5G
#$ -pe smp 16
#$ -l h_rt=9999:00:00
#$ -l h=blacklace09.blacklace
#$ -l h=blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace11.blacklace

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
