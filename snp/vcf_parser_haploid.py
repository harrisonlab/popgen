#! /usr/bin/env python
import os, sys, re, argparse
from sys import argv

## Filter VCF calls (ver. 4.2) from GATK to obtain high confidence variants
## in a haploid organism (will tweak it later to include a diploid option).
## Prints out how many SNPs filtered out at each step in the log file.

##### Default values (global per SNP call):
QUAL=40
## QUAL s is the phred-scaled probability of a SNP occurring at this site. NOT AT ALL informative in GATK!
MQ=30
## MQ gives phred-scaled probability that the read is mapped to the correct location (a low map score will occur if reads are not mapped uniquely at that site)

##### Default values (per individual genotype call):
DP=10
## Minimum depth at the site per sample
GQ=30
## GQ is the phred-scaled probability that the sample genotype being called is correct, given that there is a SNP at that site.


#### Other filters:
### As this is a haploid organism, no reads should map to the alternative allele
### (field AD). Default: on
### Only include variants with no missing data- from any sample. Default: on

script, vcf = argv
vcf_h = open(vcf)

bare = r"(\w+)(.vcf)"
out_sub = r"\1_filtered.vcf"
log_sub = r"\1_filtered.log"
out = re.sub(bare, out_sub, vcf)
log = re.sub(bare, log_sub, vcf)
vcf_out = open(out, 'w')
vcf_log = open(log, 'w')

#### Counters
# Counter for eliminated variants based on QUAL criterion
qual_c = 0
# Counter for eliminated variants based on MQ criterion
mq_c = 0
# Counter for eliminated variants based on DP criterion
dp_c = 0
# Counter for eliminated variants based on GQ criterion
gq_c = 0
# Counter for eliminated variants based on presence of alternative alleles criterion
aa_c = 0
# Counter for eliminated variants based on missing genotypes criterion
na_c = 0

def gen(fields):
    #Reset val flag
    val = 1
    mq_p = r"(MQ=)(\d+)"
    #Field with general variant stats for all samples
    general = fields[7].split(";")
    #Match up the MQ value
    m = re.search(mq_p, general[9])
    mqs = float(m.group(2))
    #Check for the QUAL score
    if float(fields[5]) < QUAL:
        global qual_c
        qual_c += 1
        val = 0
    else:
        pass
    #Check for the MQ score
    if mqs < MQ:
        global mq_c
        mq_c += 1
        val = 0
    else:
        pass

    return val

def inds(fields):
    #Reset val flag
    val = 1
    called = False
    called_dp = False
    called_gq = False
    called_aa = False
    for f in fields[9:-1]:
        #Check for the presence of missing genotypes
        if f == ".":
            global na_c
            if not called:
                na_c += 1
                called = True
            val = 0
        else:
            n = f.split(":")
            #First check for 'hidden' missing genotypes
            if n[2] == ".":
                if not called:
                    na_c += 1
                    called = True
                val = 0
            else:
                #Check for the read depth
                if float(n[2]) < DP:
                    global dp_c
                    if not called_dp:
                        dp_c += 1
                        called_dp = True
                        val = 0
                    else:
                        pass
                #Check for genotype quality
                if float(n[3]) < GQ:
                    global gq_c
                    if not called_gq:
                        gq_c += 1
                        called_gq = True
                        val = 0
                    else:
                        pass
                #Check for the presence of alternative alleles
                a = n[1].split(",")
                if (int(a[0]) and int(a[1])):
                    global aa_c
                    if not called_aa:
                        aa_c += 1
                        called_aa = True
                        val = 0
                    else:
                        pass

    return val


for line in vcf_h:
    if line.startswith("#"):
        vcf_out.write(str(line))
    else:
        fields = line.split("\t")
        if (gen(fields) == 1 and inds(fields) == 1):
            vcf_out.write(str(line))
        else:
            pass

vcf_log.write("The number of variants filtered out due to variant quality lower than %s is %s.\n" % (QUAL, qual_c))
vcf_log.write("The number of variants filtered out due to mapping qualiy lower than %s is %s.\n" % (MQ, mq_c))
vcf_log.write("The number of variants filtered out due to genotype depth lower than %s is %s.\n" % (DP, dp_c))
vcf_log.write("The number of variants filtered out due to genotype quality lower than %s is %s.\n" % (GQ, gq_c))
vcf_log.write("The number of variants filtered out due to presence of alternative allele is %s.\n" % aa_c)
vcf_log.write("The number of variants filtered out due to genotype missing altogether is %s.\n" % na_c)

vcf_h.close
vcf_out.close
vcf_log.close
