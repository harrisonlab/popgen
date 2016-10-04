#! /usr/bin/env python
from __future__ import division
from matplotlib import pyplot as plt
from collections import defaultdict as dd
from pylab import rcParams
import os, sys, re, argparse, numpy
"""
Takes in a) self-on-self pre-filtered blast results after DAGChainer analysis:
###Check if proper filtering of BLAST results in terms of gene coverage and identity has been carried out!!
(filter_blast.py -> blast_to_dagchainer.py -> dag_chainer_analysis.sh)
and b) output from GFF file processed by cds_to_chromosome_coords.py script
Outputs two types of files:
1) summary table listing gene names, number of all duplicates,
number of tandem duplicates, number of segmental duplicates
and their respective gene names
2) density plots of tandem genome duplications across
contigs, calculated in 10 kbp partitions (x axis - contig coordinate, y axis - number
of genes tandemly duplicated in the region as well as the total number of genes
in the region).
#Works only on PacBio genomes - not fragmented genomes

################### Options
--b BLAST input derieved from DAG chainer pipeline
--g Gene table input (output from cds_to_chromosome_coords.py)
--o specify if tandem duplications defined by max. distance in basepairs (argument: distance) or max. number of intermittent genes (argument: gene)
--t specify integer giving the max. number of genes or distance in basepairs (option --o) allowed for tandem duplications
"""

ap = argparse.ArgumentParser()
ap.add_argument('--b',required=True,type=str,help='BLAST input derieved from DAG chainer pipeline')
ap.add_argument('--g',required=True,type=str,help='Gene table input (output from cds_to_chromosome_coords.py)')
ap.add_argument('--o',required=True,type=str,help='How tandem duplication defined: gene or distance', choices=('gene', 'distance'))
ap.add_argument('--t',required=True,type=int,help='Max. number of genes or basepairs allowed to define as tandem duplications')

args = ap.parse_args()

metric = args.o
dist = args.t

blast = open(args.b)
table = open(args.g)

bare = r"(\w+)(.aligncoordsf)"
out_sub = r"\1_summary"
out = re.sub(bare, out_sub, args.b)
sum_out = open(out, 'w')


coo = dd(lambda: dd(int))
seg_coords = dd(lambda: dd(int))
tan_coords = dd(lambda: dd(int))

ordering = dd(lambda: dd(int))
interval = 10000

all_duplications_list = set()

stats_tandem_names =  dd(lambda: dd(str))
stats_segmental_names =  dd(lambda: dd(str))

rcParams['figure.figsize'] = 35, 15

def parse_gene_table():
    #Genes to be counted in an interval where it starts (input for duplication levels figure)
    table.seek(0)
    for line in table:
        fields = line.split()
        #Dictionary (key: contig) of dictonaries (key: interval) storing values (number of genes)
        #in the interval
        i = numpy.ceil(int(fields[2]) / interval)
        coo[fields[1]][i] += 1

def determine_gene_order():
    #dictionary of dictonaries (contigs) with keys: gene start, and values: gene_id
    #when populated sorted by key value (gene start) and created a dict of contig names and
    #gene names with values giving gene order on the contig (input for classify_duplications_g)
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

def save_tandem(gene_name1, gene_name2):
    stats_tandem_names[gene_name1][gene_name2] = "1"

def save_segmental(gene_name1, gene_name2):
    stats_segmental_names[gene_name1][gene_name2] = "2"

def classify_duplications_g():
    gene_name_full = r"(\w+)(.\w+$)"
    for line in blast:
        fields = line.split()
        i = numpy.ceil(int(fields[2]) / interval)
        #Check if genes on the same contig
        if fields[0] == fields[4]:
            gene_name_short1 = re.search(gene_name_full, fields[1])
            gene_name_short2 = re.search(gene_name_full, fields[5])
            if gene_name_short1:
                gene_pos_1 = ordering[fields[0]][gene_name_short1.group(1)]
                gene_pos_2 = ordering[fields[4]][gene_name_short2.group(1)]
            else:
                gene_pos_1 = ordering[fields[0]][fields[1]]
                gene_pos_2 = ordering[fields[4]][fields[5]]
            diff = abs(gene_pos_2 - gene_pos_1)
            if diff < dist:
                if fields[1] not in stats_tandem_names:
                    tan_coords[fields[0]][i] += 1
                save_tandem(fields[1], fields[5])
                all_duplications_list.add(fields[1])
            else:
                if fields[1] not in stats_segmental_names:
                    seg_coords[fields[0]][i] += 1
                save_segmental(fields[1], fields[5])
                all_duplications_list.add(fields[1])
        #If not, classify as segemental duplications
        else:
            if fields[1] not in stats_segmental_names:
                seg_coords[fields[0]][i] += 1
            save_segmental(fields[1], fields[5])
            all_duplications_list.add(fields[1])

def classify_duplications_d():
    for line in blast:
        fields = line.split()
        i = numpy.ceil(int(fields[2]) / interval)
        #Check if genes on the same contig
        if fields[0] == fields[4]:
        #Find out the distance between the contigs
            diff = abs(int(fields[6]) - int(fields[2]))
            if diff < dist:
                if fields[1] not in stats_tandem_names:
                    tan_coords[fields[0]][i] += 1
                save_tandem(fields[1], fields[5])
                all_duplications_list.add(fields[1])
            else:
                if fields[1] not in stats_segmental_names:
                    seg_coords[fields[0]][i] += 1
                save_segmental(fields[1], fields[5])
                all_duplications_list.add(fields[1])

        #If not, classify as segemental duplications
        else:
            if fields[1] not in stats_segmental_names:
                seg_coords[fields[0]][i] += 1
            save_segmental(fields[1], fields[5])
            all_duplications_list.add(fields[1])

def print_summary():
    sum_out.write("gene name" + "\t" + "tandem duplicate no"
    + "\t" + "tandem duplicate names" + "\t" + "segmental duplicate number" + "\t"
    + "segmental duplicate names" + "\n")
    for k in sorted(all_duplications_list):
        sum_out.write(k + "\t")
        sum_out.write(str(len(stats_tandem_names[k])) + "\t")
        for z in stats_tandem_names[k]:
            sum_out.write(z + ",")
        sum_out.write("\t")
        sum_out.write(str(len(stats_segmental_names[k])) + "\t")
        for z in stats_segmental_names[k]:
            sum_out.write(z + ",")
        sum_out.write("\n")

def print_figure():
    #Get a list of interval coordinates, and a list of all genes falling within each interval
    #Over each chromosome
    for x in sorted(coo.keys()):
        coordinates = []
        all_genes = []
        tandem_genes = []
        segmental_genes = []
    #Over each interval
        for v in sorted(coo[x].keys()):
            coordinates.append(v)
            all_genes.append(coo[x][v])
            tandem_genes.append(tan_coords[x][v])
            segmental_genes.append(seg_coords[x][v])
        plt.plot(coordinates, all_genes, color='black', marker='o', label='All genes')
        plt.plot(coordinates, segmental_genes, color='blue', marker='o', linestyle='dashed', label='Segmentally duplicated genes')
        plt.plot(coordinates, tandem_genes, color='red', marker='o', linestyle='solid', label='Tandemly duplicated genes')
        plt.title("Duplication level across " + str(x))
        plt.ylabel("Number of genes in the interval")
        plt.xlabel("Basepair *" + str(interval))
    #Location of the legen, upper right-hand corner
        plt.legend(loc=1)
        output_name = x + '.png'
        plt.savefig(output_name, bbox_inches='tight', dpi=300)
        #Clear figure
        plt.clf()

parse_gene_table()
if metric == "gene":
    determine_gene_order()
    classify_duplications_g()
elif metric == "distance":
    classify_duplications_d()
print_summary()
print_figure()


blast.close()
table.close()
sum_out.close()
