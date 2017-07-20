#!/usr/bin/python

#
# filter out progeny of all the main crosses
# except for emxfe, leave some in
# usage: zcat genotypes_input.tsv.gz | filter_progeny.py | gzip > genotypes_output.tsv.gz
#

import sys
import MySQLdb

#main database
db = MySQLdb.connect(host="mongo",   #192.168.1.100
                     user="vicker",
                     passwd=open("/home/vicker/passwords/mysql_mongo_vicker").read(),
                     db="strawberry_samples")

cc = db.cursor()

#retain some progeny of these crosses
retain = 4
retain_some =\
[
    [3,4], #emily x fenella
]

#exclude all, including parents
excl_all =\
[
    [18,19], #p1 x p2
]

#exclude all progeny with these parents
excl_progeny=\
[
    [1,2],   #redgauntlet x hapil
    [5,6],   #flamenco x chandler
    [7,8],   #bsp14 x elvira
    [14,13], #holiday x korona
    [11,12], #darselect x monterey
    [17,17], #camarosa x dover (self)
    [10,15], #capitola x cf1116
]

drop_samples = {}

#exclude progeny
for parents in excl_progeny+excl_all:
    #get list of clones matching these parents
    cc.execute('select id from clone where mat_id=%d and pat_id=%d'%(parents[0],parents[1]))
    res = cc.fetchall()

    #get matching sample ids
    for row in res:
        cc.execute('select id from sample where clone_id=%d'%int(row[0]))
        for item in cc: drop_samples[int(item[0])] = None

#retain some progeny
for parents in retain_some:
    #get list of clones matching these parents
    cc.execute('select id from clone where mat_id=%d and pat_id=%d'%(parents[0],parents[1]))
    res = cc.fetchall()

    #get matching sample ids
    for row in res[retain:]:
        cc.execute('select id from sample where clone_id=%d'%int(row[0]))
        for item in cc: drop_samples[int(item[0])] = None

#exclude parents
for parents in excl_all:
    #get the parents
    cc.execute('select id from clone where id in (%d,%d)'%(parents[0],parents[1]))
    res = cc.fetchall()

    #get matching sample ids
    for row in res:
        cc.execute('select id from sample where clone_id=%d'%int(row[0]))
        for item in cc: drop_samples[int(item[0])] = None

inpdata = sys.stdin

#filter genotype data
header_line = inpdata.readline().strip()
header = header_line.split('\t')
assert header == ['sample_id','file','clone_id','name','alias_id','marker_id','snp_id','pipeline_id','genotype']
print header_line
for line in inpdata:
    tok = line.strip().split('\t')

    #filter out flagged sampleids
    sample_id = int(tok[0])
    if sample_id in drop_samples: continue

    print line.strip()

db.close()
