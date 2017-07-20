#!/usr/bin/python
import MySQLdb as mdb

with open ("/home/sobczm/bin/mysql_sample_database/login", "r") as myfile:
    password = myfile.read().strip()

con = mdb.connect(host="mongo",   #192.168.1.100
                     user="marias",
                     passwd=password,
                     db="strawberry_copy")

with con:
    cur = con.cursor()
    cur.execute("INSERT INTO pipeline_run (name, script_path, data_path, date, notes) VALUES('haplotyping3_vesca4x', '/home/vicker/git_repos/axiom_strawberry/haplotyping3/000_affycall_and_qc_vescax4.sh', '/home/vicker/octoploid_mapping/haplotyping3/vescax4', '2017-05-23', 'the newer vescax4 ordering (the consensus map + vesca genome order combined)')")
