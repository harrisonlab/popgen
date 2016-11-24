#!/bin/bash
input=/home/sobczm/popgen/clock/coding_genomes
scripts=/home/sobczm/bin/popgen/clock/dn_ds

#Check for branch-site specific evolutionary pressure in core clock genes.
#freq, wc-1, wc-2, frh, fwd-1, ccg1, ccg4, ccg7, ccg14
#Velvet: VeA, VelB, VosA
#Photoreceptors: vvd, Phr, Phy-1

cd $input

#Download the CDS genome dataset from each species and rename the files
for a in *all.fa
do
new_name=$(echo $a | cut -d "_" -f1,2 | cut -d "." -f1)
mv $a $new_name
done

#Convert each FASTA file into single-line per sequences
for file in *
do
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $file>temp && mv temp $file
done
#Append an identifier with species id to each FASTA header
sed -i -e 's/>/>Bcit_/' Botrytis_cinerea
sed -i -e 's/>/>Cglo_/' Chaetomium_globosum
sed -i -e 's/>/>Cthe_/' Chaetomium_thermophilum
sed -i -e 's/>/>Ffuj_/' Fusarium_fujikuroi
sed -i -e 's/>/>Fgra_/' Fusarium_graminearum
sed -i -e 's/>/>Flan_/' Fusarium_langsethiae
sed -i -e 's/>/>Fpse_/' Fusarium_pseudograminearum
sed -i -e 's/>/>Fsol_/' Fusarium_solani
sed -i -e 's/>/>Fver_/' Fusarium_verticillioides
sed -i -e 's/>/>Ggra_/' Gaeumannomyces_graminis
sed -i -e 's/>/>Mory_/' Magnaporthe_oryzae
sed -i -e 's/>/>Mpoa_/' Magnaporthe_poae
sed -i -e 's/>/>Ndit_/' Neonectria_ditissima
sed -i -e 's/>/>Ncra_/' Neurospora_crassa
sed -i -e 's/>/>Ntet_/' Neurospora_tetrasperma
sed -i -e 's/>/>Pans_/' Podospora_anserina
sed -i -e 's/>/>Sbor_/' Sclerotinia_borealis
sed -i -e 's/>/>Sscl_/' Sclerotinia_sclerotiorum
sed -i -e 's/>/>Smac_/' Sordaria_macrospora
sed -i -e 's/>/>Tter_/' Thielavia_terrestris
sed -i -e 's/>/>Tatr_/' Trichoderma_atroviride
sed -i -e 's/>/>Tgam_/' Trichoderma_gamsii
sed -i -e 's/>/>Thar_/' Trichoderma_harzianum
sed -i -e 's/>/>Tree_/' Trichoderma_reesei
sed -i -e 's/>/>Tvir_/' Trichoderma_virens
sed -i -e 's/>/>Valf_/' Verticillium_alfalfae
sed -i -e 's/>/>Vdah_/' Verticillium_dahliae
sed -i -e 's/>/>Vlon_/' Verticillium_longisporum
#Extract sequences for each individual gene
python $scripts/prepare_dnds_sequences.py 
