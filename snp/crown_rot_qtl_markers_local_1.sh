#!/bin/bash
scripts=/home/sobczm/bin/popgen/snp
input=/home/sobczm/popgen/snp/snp_chip/crown_rot_local
cd $input

#Extract GFF files for SNPs within 0.5 Mbp interval of each focal SNP associated with a QTL peak.
for gff in Fa*1mbp.gff
do
intersectBed -wa -a istraw90_vesca_v1.1_snp_positions.gff3 -b $gff >${gff%.gff}_snps.gff
done
#Extract SNP names
for b in *_snps.gff
do
cut -f 18 $b | cut -d";" -f1 | sed 's/ID=//' >${b%.gff}.lst
done
#Extract genotypes for select samples. 
tr , '\n' < sample_ids_crown_rot.txt >sample_ids_crown_rot.lst
python $scripts/ananassa_genotypes_db.py sample_ids_crown_rot.lst sample_ids_crown_rot.out
#Output the genotypes in the VCF format with locations substituted  according to map positions relative to vesca 1.1. 
python $scripts/ananassa_genotypes_vcf.py sample_ids_crown_rot.out istraw90_vesca_v1.1_snp_positions.gff3

#Filter the VCF file with all markers down to only those around SNPs around QTL.
for gff in Fa*_snps.gff
do
head -2 sample_ids_crown_rot.out.vcf >${gff%.gff}.vcf
intersectBed -wa -a sample_ids_crown_rot.out.vcf -b $gff >>${gff%.gff}.vcf
done
#Concatenate gff files for SNPs of interest and then filter the main VCF file for those.
for gff in Fa*_snps.gff
do
cat $gff >>Fa_all.gff
done

head -2 sample_ids_crown_rot.out.vcf >Fa_all.vcf
intersectBed -wa -a sample_ids_crown_rot.out.vcf -b Fa_all.gff >>Fa_all.vcf

perl /home/sobczm/bin/vcftools/bin/vcf-stats Fa_all.vcf >Fa_all.stat
#Now need to select representative samples for a couple of clones with multiple samples. Criterion used: least number of missing genotypes.
#Use:
#Redgauntlet = 1113
#Hapil = 1114
#Flamenco = 859
#Holiday = 1255
#Korona = 1420  #NB this cultivar has particularly good calling rate.
#Fenella = 1203
#Emily = 1196
#Elsanta = 1341
#Honeoye = 1468
vcflib=/home/sobczm/bin/vcflib/bin
#Remove the following individuals:
for vcf in sample_ids_crown_rot.out.vcf
do
$vcflib/vcfremovesamples $vcf 1024 2079 1031 2031 1287 1130 1187 1231 1480 1415 2037 676 668 1989 841 852 >${vcf%.vcf}_nodup.vcf
done

#Classification based on resistance scores
#A)
#1-3: resistant
#3-5: intermediate
#>5: susceptible
#B)
#1-3: resistant
#>3: susceptible
#C)
#1-2: resistant
#2-4: intermediate
#>4: susceptible
#D)
#1-2: resistant
#>2: susceptible
vcftools=/home/sobczm/bin/vcftools/bin
#Calculate FST.
for vcf in *_nodup.vcf
do
$vcftools/vcftools --vcf $vcf --weir-fst-pop b_resistant.txt --weir-fst-pop b_susceptible.txt --out b_${vcf%.vcf}
done
for vcf in *_nodup.vcf
do
$vcftools/vcftools --vcf $vcf --weir-fst-pop d_resistant.txt --weir-fst-pop d_susceptible.txt --out d_${vcf%.vcf}
dones
for vcf in *_nodup.vcf
do
$vcftools/vcftools --vcf $vcf --weir-fst-pop a_resistant.txt --weir-fst-pop a_intermediate.txt --weir-fst-pop a_susceptible.txt --out a_${vcf%.vcf}
done
for vcf in *_nodup.vcf
do
$vcftools/vcftools --vcf $vcf --weir-fst-pop c_resistant.txt --weir-fst-pop c_intermediate.txt --weir-fst-pop c_susceptible.txt --out c_${vcf%.vcf}
done

#Now, remove the known individuals from EMxFE mapping population.
for vcf in *_nodup.vcf
do
$vcflib/vcfremovesamples $vcf 1595 551 554 587 601 577 568 644 675 677 >${vcf%.vcf}_noef.vcf
done

#And re-do the FST analysis above
for vcf in *_noef.vcf
do
$vcftools/vcftools --vcf $vcf --weir-fst-pop b_resistant.txt --weir-fst-pop b_susceptible.txt --out b_${vcf%.vcf}
done
for vcf in *_noef.vcf
do
$vcftools/vcftools --vcf $vcf --weir-fst-pop d_resistant.txt --weir-fst-pop d_susceptible.txt --out d_${vcf%.vcf}
done
for vcf in *_noef.vcf
do
$vcftools/vcftools --vcf $vcf --weir-fst-pop a_resistant.txt --weir-fst-pop a_intermediate.txt --weir-fst-pop a_susceptible.txt --out a_${vcf%.vcf}
done
for vcf in *_noef.vcf
do
$vcftools/vcftools --vcf $vcf --weir-fst-pop c_resistant.txt --weir-fst-pop c_intermediate.txt --weir-fst-pop c_susceptible.txt --out c_${vcf%.vcf}
done

##Repeat the same but FST for all loci across the genome.
cd $input/Genomewide
for vcf in *_noef.vcf
do
$vcftools/vcftools --vcf $vcf --weir-fst-pop b_resistant.txt --weir-fst-pop b_susceptible.txt --out b_${vcf%.vcf}
done
for vcf in *_noef.vcf
do
$vcftools/vcftools --vcf $vcf --weir-fst-pop d_resistant.txt --weir-fst-pop d_susceptible.txt --out d_${vcf%.vcf}
done
for vcf in *_noef.vcf
do
$vcftools/vcftools --vcf $vcf --weir-fst-pop a_resistant.txt --weir-fst-pop a_intermediate.txt --weir-fst-pop a_susceptible.txt --out a_${vcf%.vcf}
done
for vcf in *_noef.vcf
do
$vcftools/vcftools --vcf $vcf --weir-fst-pop c_resistant.txt --weir-fst-pop c_intermediate.txt --weir-fst-pop c_susceptible.txt --out c_${vcf%.vcf}
done

#Which outliers are witihin 1Mbp interval of QTL?
for gff in Fa*1mbp.gff
do
for bed in *.bed
do
intersectBed -wa -a $bed -b $gff >${bed%.bed}_vs_${gff%.gff}
done
done