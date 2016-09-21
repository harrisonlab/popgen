#! /usr/bin/env python
from __future__ import division
import os, sys, re, argparse
"""
Takes in self-on-self pre-filtered blast results (filter_blast.py)
Outputs three files: 1) filtred blast results with self-hits (reciprocal) eliminated
2) summary table listing gene names, number of all duplicates,
number of tandem duplicates, number of segmental duplicates
and their respective gene names 3) density plots of tandem genome duplications across
contigs, calculated in 10 kbp partitions (x axis - contig coordinate, y axis - fraction
of genes tandemly duplicated in the region as well as the total number of genes
in the region), with underlying raw data output as well.
#Works only on PacBio genomes
"""
# definition of tandem duplication: no more than 5-10 genes (changeable) seperating
# the genes or a given distance
