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
awk '$3=="gene"' $input/Maria/effectorp_secreted/v.inaequalis_172_pacbio_EffectorP_secreted.gff >$input/Maria/effectorp_secreted/v.inaequalis_172_pacbio_EffectorP_secreted_gene.gff 
#Simple intersect of different classes of genes with fixed SNPs and SVs.
#
for vcf_file in $input/Maria/vcf_files/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed_gene.vcf $input/Maria/vcf_files/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed_coding.vcf $input/Maria/vcf_files/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed_nonsyn.vcf $input/Maria/vcf_files/Ash_farm_struc_variants_fixed.vcf  
do
for gff_file in $input/Maria/all_genes/final_genes_appended_gene.gff $input/Maria/cazy/172_pacbio_CAZY_gene.gff $input/Maria/cazy_secreted/172_pacbio_CAZY_secreted_gene.gff $input/Maria/metabolite_cluster/metabolite_cluster_genes_gene.gff $input/Maria/smurf_cluster/Smurf_clusters.gff 
$input/Maria/effectorp_secreted/v.inaequalis_172_pacbio_EffectorP_secreted_gene.gff 
do
vcf=$(basename $vcf_file)
intersectBed -wb -a $vcf_file -b $gff_file > ${gff_file}_${vcf%.vcf}_overlap
done
done 

#In addition, get the transcript IDs of genes that have SSCP, and then filter the gene GFF file to output only those.
cd $input/Maria/sscp
cat 172_pacbio_sscp_headers.txt | cut -f1 >sscp_transcript_list
grep -f sscp_transcript_list $input/Maria/all_genes/final_genes_appended.gff3  | awk '$3=="mRNA"' >sscp_mrnas.gff
#Now, check for the overlap with fixed variants as for the other gene classes.
for vcf_file in $input/Maria/vcf_files/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed_gene.vcf $input/Maria/vcf_files/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed_coding.vcf $input/Maria/vcf_files/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed_nonsyn.vcf $input/Maria/vcf_files/Ash_farm_struc_variants_fixed.vcf  
do
for gff_file in sscp_mrnas.gff
do
vcf=$(basename $vcf_file)
intersectBed -wb -a $vcf_file -b $gff_file > ${gff_file}_${vcf%.vcf}_overlap
done
done

#Do the same for Swissprot annotation
cd $input/Maria/Swissprot 
cat swissprot_v2015_tophit_parsed.tbl | cut -f1 >swissprot_transcript_list
grep -f swissprot_transcript_list $input/Maria/all_genes/final_genes_appended.gff3  | awk '$3=="mRNA"' >swissprot_mrnas.gff

for vcf_file in $input/Maria/vcf_files/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed_gene.vcf $input/Maria/vcf_files/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed_coding.vcf $input/Maria/vcf_files/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed_nonsyn.vcf $input/Maria/vcf_files/Ash_farm_struc_variants_fixed.vcf  
do
for gff_file in swissprot_mrnas.gff
do
vcf=$(basename $vcf_file)
intersectBed -wb -a $vcf_file -b $gff_file > ${gff_file}_${vcf%.vcf}_overlap
done
done

#As a reference, for each file with overlaps in the swissprot directory, print a matching file listing trancript ids, swissprot hit species and swissprot hit gene function.
for overlap_file in $input/Maria/swissprot/*overlap
do
python $scripts/other/tom_passey_filter_swissprot_annotations.py $overlap_file swissprot_v2015_tophit_parsed.tbl 
done