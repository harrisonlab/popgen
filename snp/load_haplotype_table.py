#!/usr/bin/python
import sys
from sys import argv
import MySQLdb 

script, haplotype = argv

with open ("/home/sobczm/bin/mysql_sample_database/login", "r") as myfile:
    password = myfile.read().strip()

db = MySQLdb.connect(host="mongo",   #192.168.1.100
                     user="marias",
                     passwd=password,
                     db="strawberry_copy")

cur = db.cursor()
with open (haplotype) as haplotype_h:
    for line in haplotype_h:
        data = line.strip().split("\t")
        cur.execute("INSERT INTO haplotype (genotype_id, pipeline_id, chromosome, subgenome, phased_genotype) VALUES (%s, %s, %s, %s, %s)", data)
    db.commit()

cur.close()
db.close()

