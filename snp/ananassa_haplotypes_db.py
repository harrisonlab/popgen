#!/usr/bin/python
import sys
import csv
from sys import argv
import MySQLdb 

script, id_list, output = argv

with open (id_list) as id_listh:
    ids = id_listh.read().splitlines()

with open ("/home/sobczm/bin/mysql_sample_database/login", "r") as myfile:
    password = myfile.read().strip()

db = MySQLdb.connect(host="mongo",   #192.168.1.100
                     user="marias",
                     passwd=password,
                     db="strawberry_copy")

cur = db.cursor()
cur.execute("SELECT g.id, g.pipeline_id, g.sample_id, a.platform, a.probe_id, m.name, a.marker_id, h.phased_genotype, m.alleles FROM genotype g JOiN alias a ON g.alias_id = a.id JOIN marker m ON a.marker_id = m.id JOIN haplotype h on h.genotype_id = g.id WHERE g.sample_id IN (" + ",".join(("%s",) * len(ids)) + ")", tuple(ids))

num_fields = len(cur.description)
field_names = [i[0] for i in cur.description]
rows = cur.fetchall()
out = open (output, 'w')
myFile = csv.writer(out, delimiter="\t")
myFile.writerow(field_names)
myFile.writerows(rows)

out.close
    
cur.close()
db.close()
