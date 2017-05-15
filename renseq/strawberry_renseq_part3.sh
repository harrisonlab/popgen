#!/bin/bash
scripts=/home/sobczm/bin/popgen/renseq
input=/home/sobczm/popgen/renseq/strawberry

##Extract only lines matching the resistance loci in the 
#input GFF file
#First, filter the GFF file to retain only gene definitions
cd $input/rgaugury/vesca1.1/lists
awk '$3=="mRNA"' $input/genome/Fragaria_vesca_v1.1.a2.gff3 >Fragaria_vesca_v1.1.a2.mrna.gff3
for a in *.lst
do
grep -f $a Fragaria_vesca_v1.1.a2.mrna.gff3 > ${a%.lst}.gff3
done

#Check if the QTLs are within any gene?
intersectBed -wb -a qtls_helenc2.gff -b Fragaria_vesca_v1.1.a2.mrna.gff3 >qtls_helenc_overlap_mrna.gff
#All of them, apart from ID=89860579, are within a gene.

#What is the closest feature to this singleton QTL peak?
#first need to sort the Gff3 input
sort -k1,1 -k2,2n $input/genome/Fragaria_vesca_v1.1.a2.gff3 > Fragaria_vesca_v1.1.a2_sorted.gff3 
#And eliminate LG0
awk '$1!="LG0"' Fragaria_vesca_v1.1.a2_sorted.gff3  >Fragaria_vesca_v1.1.a2_sorted_elim.gff3

#Now check if any of those genes with a associated QTL is a resistance gene
for a in vesca*.gff3
do
intersectBed -wb -a qtls_helenc2.gff -b $a >${a%.gff3}_qtl_overlap
done

#Extract the CDS of the genes associated with QTLs
python $scripts/keep_list_genes2.py $input/genome/helen_qtl_associated_genes.txt $input/genome/fvesca_v1.1_all_annotated.fa No

#Now create expanded intervals around QTLs to check if they
#overlap with more resistance genes. Check 50/500 kbp upstream and downstream of each QTL peak.
awk '$4-=50000' qtls_helenc2.gff | awk '$5+=50000' | awk -v OFS="\t" '$1=$1' >qtls_helenc2_100kbp.gff

#Now check if any of those genes with a associated QTL is a resistance gene
for a in vesca*.gff3
do
intersectBed -wb -a qtls_helenc2_100kbp.gff -b $a >${a%.gff3}_qtl_overlap_100kbp
done

#Extract the CDS of the genes associated with extended intervals around QTLs
for b in *qtl_overlap_100kbp
do
cut -f 17 $b | cut -d";" -f1 | sed 's/ID=//' >${b}.lst
python $scripts/keep_list_genes2.py ${b}.lst $input/genome/fvesca_v1.1_all_annotated.fa No
done

####The same analysis, as above but for P. cactorum QTLs for Charlotte (preliminary analysis)
#Extract all gene ID within each QTLS.
intersectBed -wb -a charlotte_qtls.txt -b Fragaria_vesca_v1.1.a2.mrna.gff3 >qtls_charlotte_overlap_mrna.gff
#Only extract the gene names for matches
cat qtls_charlotte_overlap_mrna.gff | cut -f18 | cut -d ";" -f1 | sed 's/ID=//' >qtls_charlotte_overlap_mrna.lst
#Any resistance genes among them?
for a in vesca*.gff3
do
intersectBed -wb -a charlotte_qtls.txt -b $a >charlotte_${a%.gff3}_qtl_overlap
cat charlotte_${a%.gff3}_qtl_overlap | cut -f18 | cut -d ";" -f1 | sed 's/ID=//' >charlotte_${a%.gff3}.lst
done

#Extract all gene ID withins 1 Mbp of each QTLS.
intersectBed -wb -a charlotte_qtls_1mbp.txt -b Fragaria_vesca_v1.1.a2.mrna.gff3 >qtls_charlotte_1mbp_overlap_mrna.gff
#Only extract the gene names for matches
cat qtls_charlotte_1mbp_overlap_mrna.gff | cut -f18 | cut -d ";" -f1 | sed 's/ID=//' >qtls_charlotte_1mbp_overlap_mrna.lst

#Any resistance genes among them?
for a in vesca*.gff3
do
intersectBed -wb -a charlotte_qtls_1mbp.txt -b $a >charlotte_1mbp_${a%.gff3}_qtl_overlap
cat charlotte_1mbp_${a%.gff3}_qtl_overlap | cut -f18 | cut -d ";" -f1 | sed 's/ID=//' >charlotte_1mbp_${a%.gff3}.lst
done
