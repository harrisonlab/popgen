USE strawberry_copy;

DROP TABLE IF EXISTS haplotype;
CREATE TABLE haplotype
(
    haplo_id BIGINT(20) NOT NULL AUTO_INCREMENT,
    genotype_id BIGINT(20) NOT NULL,
    pipeline_id INT NOT NULL,
    chromosome VARCHAR(45) NULL,
    subgenome VARCHAR(45) NULL,
    phased_genotype VARCHAR(45) NOT NULL,
    PRIMARY KEY (haplo_id)
);
