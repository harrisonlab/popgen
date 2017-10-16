#$ -S /bin/bash
#$ -cwd
#$ -pe smp 6
#$ -l h_vmem=4G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace12.blacklace

# Changes required in the script:
# VARIABLES
# Reference - the genome reference used in read mapping.
# INSIDE THE GATK command:
# To specify which BAM mapping files (output from pre_SNP_calling_cleanup.sh, filename ending with "_rg" -> that is, with 
# read group added) are to be used in SNP calling, use the -I argument with full path to each file following after that. 
# Each new BAM file has to be specified after a separate -I
#Since our reads were in FASTA, there are no base qualities. GATK needs them, so we set them to an arbitrary 30

input=/home/sobczm/popgen/renseq/strawberry/reads/albacore_emxfe
reference=$input/gatk/smartdenovo_barcode12_fenella_trimmed_all_racon_round_10.fasta
output=$input/gatk/renseq_canu_fenella_ref.vcf

filename=$(basename "$reference")

gatk=/home/sobczm/bin/GenomeAnalysisTK-3.6

java -jar $gatk/GenomeAnalysisTK.jar -R $reference -T HaplotypeCaller -ploidy 8 -nct 6 --allow_potentially_misencoded_quality_scores --defaultBaseQualities 30 -I $input/bwa_canu_emily_all_cross/reads.sorted_unique_sorted_rg.bam -I $input/bwa_canu_fenella_all_own/reads.sorted_unique_sorted_rg.bam -o $output
