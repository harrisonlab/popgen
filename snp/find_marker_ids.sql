#
# list marker ids in order
#

#find all distinct marker ids
create or replace view
    rjvview2 as
select
    distinct a.marker_id as marker_id
from
    genotype g join
    alias a on g.alias_id = a.id;

#dump to stdout
select
    r.marker_id
from
    rjvview2 r
order by
    marker_id;
