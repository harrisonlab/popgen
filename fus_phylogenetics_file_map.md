
Copy gene models of sequenced *Fusarium oxysporum* f.sp. *cepae* and outgroup *F. proliferatum*
```
s_dir=/home/groups/harrisonlab/project_files/fusarium/gene_pred/codingquary
o_dir=/home/sobczm/popgen/input
cd $s_dir
cd F.oxysporum_fsp_cepae
cp -r 125  55  A1-2  A13  A23  A28  CB3  D2  Fus2 HB6  PG $o_dir 
cd $s_dir
cd F.proliferatum
cp -r A8 $o_dir
```
Copy gene models from *F. oxysporum* reference sequences 
```
s_dir=/home/groups/harrisonlab/project_files/fusarium/assembly/external_group
## Fusarium oxysporum fo47
cd $s_dir
cp -r F.oxysporum/fo47 $o_dir
## F.oxysporum fsp lycopersici (FoL)
cd $s_dir
cd F.oxysporum_fsp_lycopersici
cp -r 4287/Fusox1 $o_dir
```
Download remaining *F. oxysporum* genomes (CDS)
```
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
##F. oxysporum f. sp. phaseoli (FoPh)

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
