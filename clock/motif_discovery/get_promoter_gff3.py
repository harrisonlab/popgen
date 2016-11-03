#! /usr/bin/env python
import os, sys, re
from sys import argv

#Write out the CDS regions corresponding to the beginning
#of the coding region in Ensembl gff files
#And modify them so that they correspond to a given chunk of the
#putative promoter upstream of the start codon (bp range given as a 3rd argument)

#Write out plus and minus strands to separate files.

#Also, write a conversion table with output coordinates
#and protein name to cross-reference with the
#the results from GFF to FASTA script from bedtools

script, gff_file, interval = argv
gff_h = open(gff_file)

## File handles
bare = r"(\w+)(.gff3)"
plus = r"\1_plus_strand_" + re.escape(interval) + r".gff3"
minus = r"\1_minus_strand_" + re.escape(interval) + r".gff3"
coords = r"\1." + re.escape(interval)
plus_h = re.sub(bare, plus, gff_file)
minus_h = re.sub(bare, minus, gff_file)
coords_h = re.sub(bare, coords, gff_file)

out_plus = open(plus_h, 'w')
out_minus = open(minus_h, 'w')
out_coords = open(coords_h, 'w')

out_plus.write("##gff-version   3" + "\n")
out_minus.write("##gff-version   3" + "\n")

## Conversion
seen_proteins = dict()
previous_line = list()

for line in gff_h:
    if not line.startswith("#"):
        fields = line.split()
        if fields[2] == "CDS":
            ids = fields[8].split(";")
            protein_id = ids[2].split("=")
            if protein_id[1] not in seen_proteins:
                seen_proteins[protein_id[1]] = 1
                if fields[6] == "+":
                    new_start = int(fields[3]) - int(interval)
                    new_end = int(fields[3]) - 1
                    if new_start < 0:
                        new_start = 1
                    out_plus.write(fields[0] + "\t" + fields[1]  + "\t" + fields[2]  + "\t")
                    out_plus.write(str(new_start) + "\t" + str(new_end)  + "\t")
                    out_plus.write(fields[5] + "\t" + fields[6]  + "\t" + fields[7]  + "\t")
                    out_plus.write(fields[8] + "\n")
                    out_coords.write(fields[0] + ":" + str(new_start - 1) + "-" + str(new_end))
                    out_coords.write("\t" + protein_id[1] + "\n")
                if (fields[6] == "-" and len(previous_line) > 1):
                    fields2 = previous_line
                    new_start = int(fields2[4]) + 1
                    new_end = int(fields2[4]) + int(interval)
                    out_minus.write(fields2[0] + "\t" + fields2[1]  + "\t" + fields2[2]  + "\t")
                    out_minus.write(str(new_start) + "\t" + str(new_end)  + "\t")
                    out_minus.write(fields2[5] + "\t" + fields2[6]  + "\t" + fields2[7]  + "\t")
                    out_minus.write(fields2[8] + "\n")
                    out_coords.write(fields2[0] + ":" + str(new_start - 1) + "-" + str(new_end))
                    out_coords.write("\t" + protein_id[1] + "\n")
                    previous_line = fields
            else:
                if fields[6] == "-":
                    previous_line = fields


out_plus.close()
out_minus.close()
out_coords.close()
