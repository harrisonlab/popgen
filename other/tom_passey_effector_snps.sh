#!/bin/bash
input=/home/sobczm/popgen/other/passey
scripts=/home/sobczm/bin/popgen
cp -r /home/groups/harrisonlab/project_files/venturia/Maria $input

cd $input/Maria
#Annotate the file with fixed SNPs to filter out only non-synonymous variants
$scripts/summary_stats/annotate_snps_genome.sh Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed.vcf 172_pacbiov1.0

mkdir $input/Maria/vcf_files
mv *.vcf $input/Maria/vcf_files

#Prefilter the GFF files for different classes of genes to only keep lines showing gene annotation, and not individual exons etc.
awk '$3=="gene"' $input/Maria/all_genes/final_genes_appended.gff3 >$input/Maria/all_genes/final_genes_appended_gene.gff
awk '$3=="gene"' $input/Maria/cazy/172_pacbio_CAZY.gff >$input/Maria/cazy/172_pacbio_CAZY_gene.gff 
awk '$3=="gene"' $input/Maria/cazy_secreted/172_pacbio_CAZY_secreted.gff >$input/Maria/cazy_secreted/172_pacbio_CAZY_secreted_gene.gff
awk '$3=="gene"' $input/Maria/metabolite_cluster/metabolite_cluster_genes.gff >$input/Maria/metabolite_cluster/metabolite_cluster_genes_gene.gff 
#Simple intersect of different classes of genes with fixed SNPs and SVs.
#
for vcf_file in $input/Maria/vcf_files/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed_gene.vcf $input/Maria/vcf_files/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed_coding.vcf $input/Maria/vcf_files/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed_nonsyn.vcf $input/Maria/vcf_files/Ash_farm_struc_variants_fixed.vcf  
do
for gff_file in $input/Maria/all_genes/final_genes_appended_gene.gff $input/Maria/cazy/172_pacbio_CAZY_gene.gff $input/Maria/cazy_secreted/172_pacbio_CAZY_secreted_gene.gff $input/Maria/metabolite_cluster/metabolite_cluster_genes_gene.gff $input/Maria/smurf_cluster/Smurf_clusters.gff 
do
vcf=$(basename $vcf_file)
intersectBed -wb -a $vcf_file -b $gff_file > ${gff_file}_${vcf%.vcf}_overlap
done
done 