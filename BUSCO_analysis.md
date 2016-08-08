```
input=/home/sobczm/popgen/input
cd /home/sobczm/popgen/busco

##Process reference genomes
sub_BUSCO_fungi.sh $input/conglutians/Fusarium_oxysporum_f_sp_conglutinans_race_2_54008.FO_PHW808_V1.cds.all.fa
sub_BUSCO_fungi.sh $input/cubense/Fusarium_oxysporum_f_sp_cubense_race_1.Foc1_1.0.cds.all.fa
sub_BUSCO_fungi.sh $input/fo47/Fusarium_oxysporum_fo47.FO_Fo47_V1.cds.all.fa
sub_BUSCO_fungi.sh $input/fujikuroi/Fusarium_fujikuroi.EF1.cds.all.fa
sub_BUSCO_fungi.sh $input/lycopersici/Fusarium_oxysporum_f_sp_lycopersici_mn25.FO_MN25_V1.cds.all.fa
sub_BUSCO_fungi.sh $input/melonis/Fusarium_oxysporum_f_sp_melonis_26406.FO_melonis_V1.cds.all.fa
sub_BUSCO_fungi.sh $input/pisi/Fusarium_oxysporum_f_sp_pisi_hdv247.FO_HDV247_V1.cds.all.fa
sub_BUSCO_fungi.sh $input/radices-lycopersici/Fusarium_oxysporum_f_sp_radicis_lycopersici_26381.FO_CL57_V1.cds.all.fa
sub_BUSCO_fungi.sh $input/raphani/Fusarium_oxysporum_f_sp_raphani_54005.FO_PHW815_V1.cds.all.fa
sub_BUSCO_fungi.sh $input/vasinifectum/Fusarium_oxysporum_f_sp_vasinfectum_25433.FO_Cotton_V1.cds.all.fa
sub_BUSCO_fungi.sh $input/vasinifectum/Fusarium_oxysporum_f_sp_vasinfectum_25433.FO_Cotton_V1.cds.all.fa

##Process in-house genomes
```

Find the intersect of single-copy, complete genes 
