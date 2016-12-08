#!/bin/bash

input=/home/sobczm/popgen/input
scripts=/home/sobczm/bin/popgen/phylogenetics

##Process reference genomes
### Rename sequences in all FASTA files by prefixing with the species name:
cd $input/conglutians
for fasta in *.fa; do sed -i -e 's/>/>conglutinans_/' $fasta; done;
cd $input/cubense
for fasta in *.fa; do sed -i -e 's/>/>cubense_/' $fasta; done;
cd $input/fo47
for fasta in *.fa; do sed -i -e 's/>/>fo47_/' $fasta; done;
cd $input/fujikuroi
for fasta in *.fa; do sed -i -e 's/>/>fujikuroi_/' $fasta; done;
cd $input/lycopersici
for fasta in *.fa; do sed -i -e 's/>/>lycopersici_/' $fasta; done;
cd $input/melonis
for fasta in *.fa; do sed -i -e 's/>/>melonis_/' $fasta; done;
cd $input/pisi
for fasta in *.fa; do sed -i -e 's/>/>pisi_/' $fasta; done;
cd $input/radices-lycopersici
for fasta in *.fa; do sed -i -e 's/>/>radices-lycopersici_/' $fasta; done;
cd $input/raphani
for fasta in *.fa; do sed -i -e 's/>/>raphani_/' $fasta; done;
cd $input/vasinifectum
for fasta in *.fa; do sed -i -e 's/>/>vasinifectum_/' $fasta; done;
cd $input/verticillioides
for fasta in *.fa; do sed -i -e 's/>/>verticillioides_/' $fasta; done;

### Run BUSCO
cd /home/sobczm/popgen/phylogenetics/busco
qsub $scripts/sub_BUSCO2.sh $input/conglutians/Fusarium_oxysporum_f_sp_conglutinans_race_2_54008.FO_PHW808_V1.cds.all.fa Fungal
qsub $scripts/sub_BUSCO2.sh $input/cubense/Fusarium_oxysporum_f_sp_cubense_race_1.Foc1_1.0.cds.all.fa Fungal
qsub $scripts/sub_BUSCO2.sh $input/fo47/Fusarium_oxysporum_fo47.FO_Fo47_V1.cds.all.fa Fungal
qsub $scripts/sub_BUSCO2.sh $input/fujikuroi/Fusarium_fujikuroi.EF1.cds.all.fa Fungal
qsub $scripts/sub_BUSCO2.sh $input/lycopersici/Fusarium_oxysporum_f_sp_lycopersici_mn25.FO_MN25_V1.cds.all.fa Fungal
qsub $scripts/sub_BUSCO2.sh $input/melonis/Fusarium_oxysporum_f_sp_melonis_26406.FO_melonis_V1.cds.all.fa Fungal
qsub $scripts/sub_BUSCO2.sh $input/pisi/Fusarium_oxysporum_f_sp_pisi_hdv247.FO_HDV247_V1.cds.all.fa Fungal
qsub $scripts/sub_BUSCO2.sh $input/radices-lycopersici/Fusarium_oxysporum_f_sp_radicis_lycopersici_26381.FO_CL57_V1.cds.all.fa Fungal
qsub $scripts/sub_BUSCO2.sh $input/raphani/Fusarium_oxysporum_f_sp_raphani_54005.FO_PHW815_V1.cds.all.fa Fungal
qsub $scripts/sub_BUSCO2.sh $input/vasinifectum/Fusarium_oxysporum_f_sp_vasinfectum_25433.FO_Cotton_V1.cds.all.fa Fungal
qsub $scripts/sub_BUSCO2.sh $input/verticillioides/Fusarium_verticillioides.ASM14955v1.cds.all.fa Fungal
##Process in-house genomes
### Prepare input CDS assembly files (exon = CDS = cDNA identical in Braker and CodingQuary outputs)
& using the file "final_genes_combined.cdna.fasta" from each genome for phylogenetic analysis
### Rename input files in each folder by prefixing with the strain ID, and the same for sequences in all FASTA files:
cd $input/125/final
for filename in *; do mv "$filename" "125_$filename"; done;
for fasta in *.fasta; do sed -i -e 's/>/>125_/' $fasta; done;
cd $input/55/final
for filename in *; do mv "$filename" "55_$filename"; done;
for fasta in *.fasta; do sed -i -e 's/>/>55_/' $fasta; done;
cd $input/A1-2/final
for filename in *; do mv "$filename" "A1-2_$filename"; done;
for fasta in *.fasta; do sed -i -e 's/>/>A1-2_/' $fasta; done;
cd $input/A13/final
for filename in *; do mv "$filename" "A13_$filename"; done;
for fasta in *.fasta; do sed -i -e 's/>/>A13_/' $fasta; done;
cd $input/A23/final
for filename in *; do mv "$filename" "A23_$filename"; done;
for fasta in *.fasta; do sed -i -e 's/>/>A23_/' $fasta; done;
cd $input/A28/final
for filename in *; do mv "$filename" "A28_$filename"; done;
for fasta in *.fasta; do sed -i -e 's/>/>A28_/' $fasta; done;
cd $input/CB3/final
for filename in *; do mv "$filename" "CB3_$filename"; done;
for fasta in *.fasta; do sed -i -e 's/>/>CB3_/' $fasta; done;
cd $input/D2/final
for filename in *; do mv "$filename" "D2_$filename"; done;
for fasta in *.fasta; do sed -i -e 's/>/>D2_/' $fasta; done;
cd $input/Fus2_canu_new/final
for filename in *; do mv "$filename" "Fus2_$filename"; done;
for fasta in *.fasta; do sed -i -e 's/>/>Fus2_/' $fasta; done;
cd $input/HB6/final
for filename in *; do mv "$filename" "HB6_$filename"; done;
for fasta in *.fasta; do sed -i -e 's/>/>HB6_/' $fasta; done;
cd $input/PG/final
for filename in *; do mv "$filename" "PG_$filename"; done;
for fasta in *.fasta; do sed -i -e 's/>/>PG_/' $fasta; done;
cd $input/proliferatum/final
for filename in *; do mv "$filename" "proliferatum_$filename"; done;
for fasta in *.fasta; do sed -i -e 's/>/>proliferatum_/' $fasta; done;

## Run BUSCO
cd /home/sobczm/popgen/busco
qsub $scripts/sub_BUSCO2.sh $input/125/final/125_final_genes_combined.cdna.fasta Fungal
qsub $scripts/sub_BUSCO2.sh $input/55/final/55_final_genes_combined.cdna.fasta Fungal
qsub $scripts/sub_BUSCO2.sh $input/A1-2/final/A1-2_final_genes_combined.cdna.fasta Fungal
qsub $scripts/sub_BUSCO2.sh $input/A13/final/A13_final_genes_combined.cdna.fasta Fungal
qsub $scripts/sub_BUSCO2.sh $input/A23/final/A23_final_genes_combined.cdna.fasta Fungal
qsub $scripts/sub_BUSCO2.sh $input/A28/final/A28_final_genes_combined.cdna.fasta Fungal
qsub $scripts/sub_BUSCO2.sh $input/CB3/final/CB3_final_genes_combined.cdna.fasta Fungal
qsub $scripts/sub_BUSCO2.sh $input/D2/final/D2_final_genes_combined.cdna.fasta Fungal
qsub $scripts/sub_BUSCO2.sh $input/Fus2_canu_new/final/Fus2_final_genes_combined.cdna.fasta Fungal
qsub $scripts/sub_BUSCO2.sh $input/HB6/final/HB6_final_genes_combined.cdna.fasta Fungal
qsub $scripts/sub_BUSCO2.sh $input/PG/final/PG_final_genes_combined.cdna.fasta Fungal
qsub $scripts/sub_BUSCO2.sh $input/proliferatum/final/proliferatum_final_genes_combined.cdna.fasta Fungal

## Find the intersect of single-copy, complete genes across all the genomes being analysed.
# This command should generate a separate FASTA file for each BUSCO complete single copy gene 
# conserved across all your genomes analysed. You use those files as input to be aligned by MAFFT in the next step.

### Create a list of all BUSCO IDs

pushd /home/sobczm/bin/BUSCO_v1.22/fungi/hmms
ls -1 | sed -e 's/\..*$//' >../all_buscos_fungi
cp ../all_buscos_fungi /home/sobczm/popgen/phylogenetics/busco
popd

### Iteratively find the intersect of IDs of all 'complete' BUSCO genes present in the runs in the current directory

#!/bin/bash
cat all_buscos_fungi >temp_ref
for d in $PWD/run*
do
    if test -n "$(find $d -maxdepth 1 -name 'full_table*' -print -quit)"
    then
        awk '$2 == "Complete" { print $1}' $d/full_table* >temp
        grep -Fx -f temp temp_ref >final_list_ssc
        cat final_list_ssc >temp_ref
    else
        echo "There is no full_table result in $d."
    fi
done

### Pick random 50 IDs of all 'complete' BUSCO genes present in the intersect
sort -R final_list_ssc | head -n 50 >align_input_list.txt

#Note to self: having done this, the nucleotide diversity of 50 random genes is very consistent, with most of them
#in the range of 0.1-0.15. so will have to do the following steps for ALL shared BUSCOs. 

#Most likely this scenario will happen to you as well, so in most cases you will have to align all genes
#contained in final_list_ssc. In that case, just do this, instead of the above command:

cat final_list_ssc >align_input_list.txt

## Create FASTA files with separate alignment input for each of the selected genes.
#This command should generate a separate FASTA file for each BUSCO complete single copy gene conserved across all your genomes analysed. 
#You use those files as input to be aligned by MAFFT in the next step.
# Before executing the script, make sure you are in the directory with all the relevant BUSCO run results dirs.
perl $scripts/get_alignments.pl
