#! /usr/bin/env python
import os, sys, re, argparse
from collections import defaultdict as dd

#Filter out variants in VCF file to keep only those showing fixation or near-fixation
#(as per allele frequency threshold) in two different populations of samples.

######################## 
#   IMPORTANT REMARKS
########################

#Missing genotypes ignored in counting up allele frequencis, so advised to pre-filter the VCF file for variants with lots of #missing genotypes
 
#BUT script can handle >2 alleles per variant, so no prefiltering required for biallelic variants. 

ap = argparse.ArgumentParser()
ap.add_argument('--vcf',required=True,type=str,help='Input VCF file with variants')
ap.add_argument('--out',required=True,type=str,help='Name of the output VCF file')
ap.add_argument('--ply',required=True,type=int,help='Ploidy of the organism')
ap.add_argument('--pop1',required=True,type=str,help='Population A sample names, seperated by comma (,)')
ap.add_argument('--pop2',required=True,type=str,help='Population B sample names, seperated by comma (,)')
ap.add_argument('--thr',required=True,type=float,help='Min. allele frequency threshold', default=0.95)
args = ap.parse_args()

#population A sample names
pop1_names=args.pop1.split(",,")
#population B sample names
pop2_names=args.pop2.split(",,")
#Set ploidy
ploidy = args.ply
#Set allele frequency threshold
threshold = args.thr
#Input VCF
input_vcf = args.vcf
#Output VCF
output_vcf = args.out
file_h = open(output_vcf, 'w')

#Indexes of pop samples
pop1_indexes = list()
pop2_indexes = list()



def iterate_genotypes(indexes, fields):
    #Empty dictonaries to count up the genotypes
    pop_dict = dd(int)
    for a in indexes:
        #Ignore empty genotypes
        if fields[a] == ".":
            pass        
        else:
            n = fields[a].split(":")
       #Check for hidden empty genotypes:
            hidden_empty = re.match(r'\d', n[0])
            if hidden_empty:
                if ploidy == 1:
                    pop_dict[n[0]] += 1
                elif ploidy == 2:
                    gens = n[0].split("/")
                    for g in gens:
                        pop_dict[g] += 1
                else:
                    exit
            else:
                pass
    return pop_dict
def vcf_handling():
    vcf_h = open(input_vcf)
    for line in vcf_h:
        if line.startswith("##"):
            file_h.write(str(line))
        elif line.startswith("#"):
            file_h.write(str(line))
            fields = line.split()
            #Index all the individuals
            for f in fields[9:]:
                if f in pop1_names:
                    pop1_indexes.append(fields.index(f))
                elif f in pop2_names:
                    pop2_indexes.append(fields.index(f))
                else:
                    pass
        else:
            fields = line.split()
            #Iterate over population genotypes
            pop1_dict = iterate_genotypes(pop1_indexes, fields)
            #print(pop1_dict)
            pop2_dict = iterate_genotypes(pop2_indexes, fields)
            #Sum the number of genotypes available in each population
            pop1_total = 0
            pop2_total = 0
            for v in pop1_dict.values():
                pop1_total += v
            for v in pop2_dict.values():
                pop2_total += v
            #Check if number of genotypes >0 in both populations
            if (pop1_total > 0 and pop2_total > 0):
                #Get the most common genotype
                pop1_top = max(pop1_dict, key=pop1_dict.get)
                pop2_top = max(pop2_dict, key=pop2_dict.get)
                #Check that different top alleles present in each population
                if pop1_top != pop2_top:
                    pop1_freq = pop1_dict[pop1_top]/pop1_total 
                    pop2_freq = pop2_dict[pop2_top]/pop2_total
                    #Check that top allele frequency higher than threshold in both populations
                    if (pop1_freq >= threshold and pop2_freq >= threshold):
                        #Print the selected fixed variants to file
                        #print ("Writing fixed variant to file.")
                        file_h.write(str(line))    
    vcf_h.close()
    file_h.close()

vcf_handling()
