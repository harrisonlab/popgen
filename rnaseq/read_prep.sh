#!/bin/bash
input=/home/sobczm/popgen/rnaseq
scripts=/home/sobczm/bin/popgen/rnaseq

#Read prep of RNA-Seq reads from the P.cactorum/F.ananassa experiment.

#Raw Data QC
for RawData in $(ls $input/*/*.fastq.gz); 
do
sh $scripts/run_fastqc.sh $RawData
done

#Quality/adaptor trimming
for reads in $(ls $input/*/*1.fastq.gz); 
do
    Jobs=$(qstat | grep 'sub_read_q' | wc -l)
    while [ $Jobs -gt 10 ]
    do
        sleep 10
        printf "."
        Jobs=$(qstat | grep 'sub_read_q' | wc -l)
    done
reads2=$(echo "$reads" | sed 's/1.fastq.gz/2.fastq.gz/')
qsub $scripts/sub_read_qc.sh $reads $reads2
done
## After all, simply copying files with the same QC carried out by Andy
cp -r /home/groups/harrisonlab/project_files/idris/qc_rna/paired/Transcriptome_Emily_Fenella_Pcactorum-2017-04-07 $input
#Copy over the vesca1.1 assembly and P.cactorum assemblies with GFF annotation.
cp -r /home/sobczm/popgen/renseq/strawberry/genome/fvesca_v1.1_all.fa $input
cp -r /home/sobczm/popgen/renseq/strawberry/genome/Fragaria_vesca_v1.1.a2.gff3 $input
cp -r /home/groups/harrisonlab/project_files/idris/assembly/merged_canu_spades/P.cactorum/414_v2/filtered_contigs/contigs_min_500bp_renamed.fasta $input
cp -r /home/groups/harrisonlab/project_files/idris/gene_pred/final_genes/P.cactorum/414_v2/final/final_genes_appended.gff3 $input
###Alignment of all reads with STAR - TEST with Andy's script
##To vesca1.1
cd $input
qsub $scripts/sub_star.sh fvesca_v1.1_all.fa Sample_2212_LIB26247_LDI23278/F/PRO1467_S1_totRNA_S17_L001_R1_trim.fq.gz Sample_2212_LIB26247_LDI23278/R/PRO1467_S1_totRNA_S17_L001_R2_trim.fq.gz test_vesca Fragaria_vesca_v1.1.a2.gff3
for fq_f in $(ls $PWD/*/F/*.fq.gz)
do
fq_fb=$(basename "$fq_f")
fq_fd=$(dirname "$fq_f")
fq_r="${fq_fb%_R1_trim.fq.gz}_R2_trim.fq.gz"
fq_p="${fq_fd%/F}/R"
qsub $scripts/sub_star_sensitive.sh $input/fvesca_v1.1_all.fa $fq_f $fq_p/$fq_r "vesca_${fq_fb%_R1_trim.fq.gz}" $input/Fragaria_vesca_v1.1.a2.gff3
done
##To P.cactorum
qsub $scripts/sub_star.sh contigs_min_500bp_renamed.fasta Sample_2212_LIB26247_LDI23278/F/PRO1467_S1_totRNA_S17_L001_R1_trim.fq.gz Sample_2212_LIB26247_LDI23278/R/PRO1467_S1_totRNA_S17_L001_R2_trim.fq.gz test_pcac final_genes_appended.gff3
for fq_f in $(ls $PWD/*/F/*.fq.gz)
do
fq_fb=$(basename "$fq_f")
fq_fd=$(dirname "$fq_f")
fq_r="${fq_fb%_R1_trim.fq.gz}_R2_trim.fq.gz"
fq_p="${fq_fd%/F}/R"
qsub $scripts/sub_star_sensitive.sh $input/contigs_min_500bp_renamed.fasta $fq_f $fq_p/$fq_r "pcac_${fq_fb%_R1_trim.fq.gz}" $input/final_genes_appended.gff3
done

#Remove not needed mapping files from the first round of mappind
for n in vesca_*/ pcac_*/
do
cd $n
rm *.bam
cd ../
done

###Cross-alignment of unmapped reads
##To vesca 1.1
for k in pcac_*/
do
cd $k
qsub $scripts/sub_star_sensitive.sh $input/fvesca_v1.1_all.fa star_aligmentUnmapped.out.mate1 star_aligmentUnmapped.out.mate2 vesca $input/Fragaria_vesca_v1.1.a2.gff3
cd ../
done

##To P.cactorum
for k in vesca_*/
do
cd $k
qsub $scripts/sub_star_sensitive.sh $input/contigs_min_500bp_renamed.fasta star_aligmentUnmapped.out.mate1  star_aligmentUnmapped.out.mate2 pcac $input/final_genes_appended.gff3
cd ../
done
