
Copy gene models of sequenced *Fusarium oxysporum* f.sp. *cepae* and outgroup *F. proliferatum*
```
s_dir=/home/groups/harrisonlab/project_files/fusarium/gene_pred/codingquary/
o_dir=/home/sobczm/popgen/input
cd $s_dir/F.oxysporum_fsp_cepae
cp -r 125  55  A1-2  A13  A23  A28  CB3  D2 Fus2_canu_new HB6  PG $o_dir 
cd $s_dir/F.oxysporum_fsp_cepae
cp -r Fus2
cd $s_dir/F.proliferatum
cp -r A8 $o_dir
mv $o_dir/A8 $o_dir/proliferatum
```

Download remaining *F. oxysporum* genomes (CDS)
```
## Fusarium oxysporum fo47 - assembly FO_Fo47_V1
cd $o_dir
mkdir fo47
cd fo47
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-32/fasta/fungi_ascomycota1_collection/fusarium_oxysporum_fo47/cds/Fusarium_oxysporum_fo47.FO_Fo47_V1.cds.all.fa.gz
gunzip -d Fusarium_oxysporum_fo47.FO_Fo47_V1.cds.all.fa.gz

## Fusarium oxysporum f. sp. lycopersici - assembly FO_MN25_V1
cd $o_dir
mkdir lycopersici
cd lycopersici
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-32/fasta/fungi_ascomycota1_collection/fusarium_oxysporum_f_sp_lycopersici_mn25/cds/Fusarium_oxysporum_f_sp_lycopersici_mn25.FO_MN25_V1.cds.all.fa.gz
gunzip -d Fusarium_oxysporum_f_sp_lycopersici_mn25.FO_MN25_V1.cds.all.fa.gz

##F. oxysporum f. sp. raphani - assembly FO_PHW815_V1
cd $o_dir
mkdir raphani
cd raphani
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-32/fasta/fungi_ascomycota1_collection/fusarium_oxysporum_f_sp_raphani_54005/cds/Fusarium_oxysporum_f_sp_raphani_54005.FO_PHW815_V1.cds.all.fa.gz
gunzip -d Fusarium_oxysporum_f_sp_raphani_54005.FO_PHW815_V1.cds.all.fa.gz

##F. oxysporum f. sp. pisi (FoPi) - assembly FO_HDV247_V1
cd $o_dir
mkdir pisi
cd pisi
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-32/fasta/fungi_ascomycota1_collection/fusarium_oxysporum_f_sp_pisi_hdv247/cds/Fusarium_oxysporum_f_sp_pisi_hdv247.FO_HDV247_V1.cds.all.fa.gz
gunzip -d Fusarium_oxysporum_f_sp_pisi_hdv247.FO_HDV247_V1.cds.all.fa.gz

##F. oxysporum f. sp. radices-lycopersici (FoRL) - assembly FO_CL57_V1
cd $o_dir
mkdir radices-lycopersici
cd radices-lycopersici
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-32/fasta/fungi_ascomycota1_collection/fusarium_oxysporum_f_sp_radicis_lycopersici_26381/cds/Fusarium_oxysporum_f_sp_radicis_lycopersici_26381.FO_CL57_V1.cds.all.fa.gz
gunzip -d Fusarium_oxysporum_f_sp_radicis_lycopersici_26381.FO_CL57_V1.cds.all.fa.gz

##F. oxysporum f. sp. cubense (FoCub) - assembly Foc1_1.0
cd $o_dir
mkdir cubense
cd cubense
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-32/fasta/fungi_ascomycota1_collection/fusarium_oxysporum_f_sp_cubense_race_1/cds/Fusarium_oxysporum_f_sp_cubense_race_1.Foc1_1.0.cds.all.fa.gz
gunzip -d Fusarium_oxysporum_f_sp_cubense_race_1.Foc1_1.0.cds.all.fa.gz 

##F. oxysporum f. sp. vasinifectum (FoV) - assembly FO_Cotton_V1
cd $o_dir
mkdir vasinifectum
cd vasinifectum
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-32/fasta/fungi_ascomycota1_collection/fusarium_oxysporum_f_sp_vasinfectum_25433/cds/Fusarium_oxysporum_f_sp_vasinfectum_25433.FO_Cotton_V1.cds.all.fa.gz
gunzip -d Fusarium_oxysporum_f_sp_vasinfectum_25433.FO_Cotton_V1.cds.all.fa.gz

##F. oxysporum f. sp. melonis (FoM) - assembly FO_melonis_V1
cd $o_dir
mkdir melonis
cd melonis
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-32/fasta/fungi_ascomycota1_collection/fusarium_oxysporum_f_sp_melonis_26406/cds/Fusarium_oxysporum_f_sp_melonis_26406.FO_melonis_V1.cds.all.fa.gz
gunzip -d Fusarium_oxysporum_f_sp_melonis_26406.FO_melonis_V1.cds.all.fa.gz

##F. oxysporum f. sp. conglutinans (FoCo) - assembly FO_PHW808_V1 
cd $o_dir
mkdir conglutians 
cd conglutians
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-32/fasta/fungi_ascomycota1_collection/fusarium_oxysporum_f_sp_conglutinans_race_2_54008/cds/Fusarium_oxysporum_f_sp_conglutinans_race_2_54008.FO_PHW808_V1.cds.all.fa.gz
gunzip -d Fusarium_oxysporum_f_sp_conglutinans_race_2_54008.FO_PHW808_V1.cds.all.fa.gz

```

Download remaining outgroup genomes (CDS)
```
##F. fujikuroi  (assembly EF 1)
cd $o_dir
mkdir fujikuroi
cd fujikuroi
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-32/fasta/fusarium_fujikuroi/cds/Fusarium_fujikuroi.EF1.cds.all.fa.gz
gunzip -d Fusarium_fujikuroi.EF1.cds.all.fa.gz

##F. verticillioides (assembly ASM14955v1) 
cd $o_dir
mkdir verticillioides
cd verticillioides
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-32/fasta/fusarium_verticillioides/cds/Fusarium_verticillioides.ASM14955v1.cds.all.fa.gz
gunzip -d Fusarium_verticillioides.ASM14955v1.cds.all.fa.gz
```

Copy the Bowtie2 mappings and the reference (Fus2 unmasked genome)
```
mkdir mappings
cp /home/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs_repmask/Fus2_canu_contigs_unmasked.fa $o_dir/mappings
bowtie_mp=/home/groups/harrisonlab/project_files/fusarium/analysis/genome_alignment/bowtie/F.oxysporum_fsp_cepae/
cp -r $bowtie_mp/125/vs_Fus2_unmasked_max1200 $o_dir/mappings/125
cp -r $bowtie_mp/55/vs_Fus2_unmasked_max1200 $o_dir/mappings/55
cp -r $bowtie_mp/A1-2/vs_Fus2_unmasked_max1200 $o_dir/mappings/A1-2
cp -r $bowtie_mp/A13/vs_Fus2_unmasked_max1200 $o_dir/mappings/A13
cp -r $bowtie_mp/A23/vs_Fus2_unmasked_max1200 $o_dir/mappings/A23
cp -r $bowtie_mp/A28/vs_Fus2_unmasked_max1200 $o_dir/mappings/A28
cp -r $bowtie_mp/CB3/vs_Fus2_unmasked_max1200 $o_dir/mappings/CB3
cp -r $bowtie_mp/D2/vs_Fus2_unmasked_max1200 $o_dir/mappings/D2
cp -r $bowtie_mp/HB6/vs_Fus2_unmasked_max1200 $o_dir/mappings/HB6
cp -r $bowtie_mp/PG/vs_Fus2_unmasked_max1200 $o_dir/mappings/PG

cp -r $bowtie_mp/Fus2/vs_Fus2_unmasked_mixed $o_dir/mappings/Fus2

```
Copy Andy's gene annotations
```
mkdir annotations && annotations
cp -r $fus/gene_pred/annotations/F.oxysporum/fo47/fo47_gene_annotations.tab ./
cp -r $fus/gene_pred/annotations/F.oxysporum_fsp_cepae/125/125_gene_annotations.tab ./
cp -r $fus/gene_pred/annotations/F.oxysporum_fsp_cepae/A13/A13_gene_annotations.tab ./
cp -r $fus/gene_pred/annotations/F.oxysporum_fsp_cepae/A23/A23_gene_annotations.tab ./
cp -r $fus/gene_pred/annotations/F.oxysporum_fsp_cepae/A28/A28_gene_annotations.tab ./
cp -r $fus/gene_pred/annotations/F.oxysporum_fsp_cepae/CB3/CB3_gene_annotations.tab ./
cp -r $fus/gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab ./
cp -r $fus/gene_pred/annotations/F.oxysporum_fsp_cepae/PG/PG_gene_annotations.tab ./
cp -r $fus/gene_pred/annotations/F.oxysporum_fsp_lycopersici/4287/4287_gene_annotations.tab ./

```

