USE strawberry_copy;

LOAD DATA INFILE '/home/sobczm/popgen/snp/snp_chip/haplotypes/haplotype_table_keys.txt' 
INTO TABLE haplotype
FIELDS TERMINATED BY '\t' 
LINES TERMINATED BY '\n'
