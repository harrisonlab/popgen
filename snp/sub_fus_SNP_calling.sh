#$ -S /bin/bash
#$ -cwd 
#$ -pe smp 1
#$ -l h_vmem=16G 
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace11.blacklace|blacklace12.blacklace

input=/home/sobczm/popgen/input/mappings

reference=Fus2_canu_contigs_unmasked.fa
filename=$(basename "$reference")
output1="${filename%.*}.dict" 
output2="${filename%.*}.vcf" 

gatk=/home/sobczm/bin/GenomeAnalysisTK-3.6

java -jar $gatk/GenomeAnalysisTK.jar \
     -R $reference \
     -T HaplotypeCaller \
     -I $input/125_Fus2_canu_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/55_Fus2_canu_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/A1-2_Fus2_canu_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/A13_Fus2_canu_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/A23_Fus2_canu_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/A28_Fus2_canu_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/CB3_Fus2_canu_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/D2_Fus2_canu_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/Fus2_Fus2_canu_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/HB6_Fus2_canu_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/PG_Fus2_canu_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
	   -o $output2
	    
	    
#####################################
# Notes on GATK parallelisation
#####################################
# http://gatkforums.broadinstitute.org/gatk/discussion/1975/how-can-i-use-parallelism-to-make-gatk-tools-run-faster
