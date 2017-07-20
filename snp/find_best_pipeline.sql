#
# choose the best pipeline per sample
#

#where more than one pipeline is available for a sample
#simply pick the one with the smallest pipeline_id as this is sufficient
#to select the consensusmap pipeline in preference to the haplotype pipeline
create or replace view
    rjvview1 as
select
    sample_id,
    min(pipeline_id) as best_id,
    count(distinct pipeline_id) as npipelines,
    group_concat(distinct pipeline_id order by pipeline_id separator ',') as pipeline_list
from
    genotype
group by
    sample_id;

#dump results to stdout
select
    r.sample_id,
    r.best_id
from
    rjvview1 r;
