#! /usr/bin/env python
from __future__ import division
from matplotlib import pyplot as plt
from collections import defaultdict as dd
import os, sys, re, argparse, numpy
"""
Takes in a) self-on-self pre-filtered blast results after DAGChainer analysis
(filter_blast.py -> blast_to_dagchainer.py -> dag_chainer_analysis.sh)
and b) output from GFF file processed by cds_to_chromose_coords.py script
Outputs three files:
1) summary table listing gene names, number of all duplicates,
number of tandem duplicates, number of segmental duplicates
and their respective gene names
2) density plots of tandem genome duplications across
contigs, calculated in 100 kbp partitions (x axis - contig coordinate, y axis - number
of genes tandemly duplicated in the region as well as the total number of genes
in the region), with underlying raw data output as well (3rd file).
#Works only on PacBio genomes - not fragmented genomes

################### Options
--b BLAST input
--g Gene table input
--o specify if tandem duplications defined by max. distance in basepairs (argument: distance) or max. number of intermittent genes (argument: gene)
--t specify integer specifying the max. number of genes or distance in basepairs (option --o) allowed for tandem duplications
"""

###Check if proper filtering of BLAST results in terms of gene coverage and id has been carried out!!

ap = argparse.ArgumentParser()
ap.add_argument('--b',required=True,type=str,help='BLAST input')
ap.add_argument('--g',required=True,type=str,help='Gene table input')
ap.add_argument('--o',required=True,type=str,help='How tandem duplication defined: gene or distance', choices=('gene', 'distance'))
ap.add_argument('--t',required=True,type=int,help='Max. number of genes or basepairs allowed for tandem duplications')

args = ap.parse_args()

metric = args.o
dist = args.t

blast = open(args.b)
table = open(args.g)

bare = r"(\w+)(.aligncoordsf)"
out_sub = r"\1_summary"
raw_sub = r"\1_raw"
out = re.sub(bare, out_sub, args.b)
out2 = re.sub(bare, raw_sub, args.b)
sum_out = open(out, 'w')
raw_out = open(out2, 'w')

coo = dd(lambda: dd(int))
ordering = dd(lambda: dd(int))
interval = 10000

all_duplications_list = []
seg_coords = []
tan_coords = []

stats_tandem = {}
stats_segmental = {}

def determine_gene_order():
    #dictionary of dictonaries(contigs) with keys: gene start and values: gene_id
    #when populated sorted by key value and created a dict of contig names and
    #gene names with values giving gene order on the contig
    table.seek(0)
    first_pass = dd(lambda: dd(int))
    for line in table:
        fields = line.split()
        first_pass[fields[1]][int(fields[2])] = fields[0]
        #Iterate over each chromosome
        for x in sorted(first_pass.keys()):
            count = 0
            #Iterate over each gene
            for v in sorted(first_pass[x].keys()):
                count += 1
                ordering[x][first_pass[x][v]] = count
                print (first_pass[x][v], count)


def save_tandem(gene_name1, gene_name2):
    stats_tandem[gene_name1][0] += 1
    stats_tandem[gene_name1][1].append(gene_name2)

def save_segmental(gene_name1, gene_name2):
    stats_segmental[gene_name1][0] += 1
    stats_segmental[gene_name1][1].append(gene_name2)

def classify_duplications_g():
    for line in blast:
        fields = line.split()
        #Check if genes on the same contig
        if fields[0] == fields[4]:
            gene_pos_1 = ordering[field[0]][field[1]]
            gene_pos_2 = ordering[field[4]][field[5]]
            diff = abs(gene_pos_2 - gene_pos_1)
            if diff <= dist:
                all_duplications_list.append(fields[1])
                seg_coords.append[fields[6]]
                tan_coords.append[fields[6]]
                save_tandem(fields[1], fields[5])
            else:
                all_duplications_list.append(fields[1])
                save_segmental(fields[1], fields[5])
        #If not, classify as segemental duplications
        else:
            all_duplications_list.append(fields[1])
            seg_coords.append[fields[6]]
            save_segmental(fields[1], fields[5])

def classify_duplications_d():
    for line in blast:
        fields = line.split()
        #Check if genes on the same contig
        if fields[0] == fields[4]:
        #Find out the distance between the contigs
            diff = abs(fields[6] - fields[2])
            if diff <= dist:
                all_duplications_list.append(fields[1])
                tan_coords.append[fields[6]]
                save_tandem(fields[1], fields[5])
            else:
                all_duplications_list.append(fields[1])
                save_segmental(fields[1], fields[5])
                seg_coords.append[fields[6]]

        #If not, classify as segemental duplications
        else:
            all_duplications_list.append(fields[1])
            seg_coords.append[fields[6]]
            save_segmental(fields[1], fields[5])

def print_summary():
    sum_out.write("gene name" + "\t" + "tandem duplicate no"
    + "\t" + "tandem duplicate names" + "\t" + "segmental duplicate number" + "\t"
    + "segmental duplicate names" + "\n")
    for k in sorted(all_duplications_list):
    sum_out.write(k + "\t" + stats_tandem[k][0] "\t" + stats_tandem[k][1]
    + stats_segmental[k][0] + stats_segmental[k][1])


def number_genes():
    return "number_genes"

def distance():
    return "distance"

def terminate():
    sys.exit(0)

def parse_gene_table():
    #Gene to be counted in an interval where it starts
    table.seek(0)
    for line in table:
        fields = line.split()
        #Dictionary (key contig) of dictonaries (key interval) storing values (number of genes)
        #in the interval
        i = numpy.ceil(int(fields[2]) / interval)
        coo[fields[1]][i] += 1
    for x in sorted(coo.keys()):
        for v in sorted(coo[x].keys()):
            #print(coo[x][v])

def print_figure():
    #Get rid of repeated values - each gene should appear just once
    seg = set(seg_coords)
    tan = set(tan_coords)
    for x in sorted(coo.keys()):
        coordinates = []
        all_genes = []
        for v in sorted(coo[x].keys()):
            coordinates.append(v)
            all_genes.append(coo[x][v])
    #number of elements in both x and y lists have to match
    #coordinates = ["0", "10", "20", "30", "40", "50", "60", "70", "80", "90"]
    #all_genes = ["10", "5", "16", "80", "20", "4", "14", "50", "11", "2"]
    #tandem_duplications = ["0", "3", "0", "10", "1", "2", "7", "1", "3", "0"]
        #plt.scatter(coordinates, all_genes, color='black', marker='o', label='All genes')
        plt.plot(coordinates, all_genes, color='black', marker='o', label='All genes')
    #plt.plot(coordinates, all_duplications, color='blue', marker='o', linestyle='dashed', label='All duplicated genes')
    #plt.plot(coordinates, tandem_duplications, color='red', marker='o', linestyle='solid', label='Tandemly duplicated genes')
        plt.title("Duplication level across " + str(x))
        plt.ylabel("Number of genes in the interval")
        plt.xlabel("Basepair *" + str(interval))
    #Location of the legen, upper right-hand corner
        plt.legend(loc=1)
        output_name = x + '.png'
        plt.savefig(output_name, bbox_inches='tight')
        #Clear figure
        plt.clf()


determine_gene_order()
if metric == gene:
    classify_duplications_g()
elif metric == distance:
    classify_duplications_d
parse_gene_table()
print_summary()
#print_figure()


blast.close()
table.close()
sum_out.close()
raw_out.close()
