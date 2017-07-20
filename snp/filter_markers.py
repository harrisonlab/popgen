#!/usr/bin/python

#
# filter out markers where marker is flagged as:
# multiallelic, sexerror or crosslg
# usage: zcat genotypes_input.tsv.gz | filter_markers.py | gzip > genotypes_output.tsv.gz
#

import sys
import MySQLdb

#main database
db = MySQLdb.connect(host="mongo",   #192.168.1.100
                     user="vicker",
                     passwd=open("/home/vicker/passwords/mysql_mongo_vicker").read(),
                     db="strawberry_samples")

cc = db.cursor()

#get list of markers tagged as multiallelic, crosslg or sexerror
cc.execute('select marker_id from marker_tag where tag in ("sexerror","multiallelic","crosslg")')

drop_markers = { int(x[0]):None for x in cc }

inpdata = sys.stdin

#filter genotype data
header_line = inpdata.readline().strip()
header = header_line.split('\t')
assert header == ['sample_id','file','clone_id','name','alias_id','marker_id','snp_id','pipeline_id','genotype']
print header_line
for line in inpdata:
    tok = line.strip().split('\t')

    #filter out flagged markers
    marker_id = int(tok[5])
    if marker_id in drop_markers: continue

    print line.strip()

db.close()
