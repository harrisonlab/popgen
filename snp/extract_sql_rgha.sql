#
# dump genotypes of all samples as a list to csv
#

select
   g.sample_id,
   s.file,
   s.clone_id,
   c.name,
   g.alias_id,
   a.marker_id,
   m.snp_id,
   g.pipeline_id,
   g.genotype
from genotype g
     join alias a on g.alias_id = a.id
     join marker m on a.marker_id = m.id
     join sample s on g.sample_id = s.id
     join clone c on s.clone_id = c.id
where s.file like '%_RGxHA%.CEL';