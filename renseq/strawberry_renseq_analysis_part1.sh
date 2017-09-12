#!/bin/bash
scripts=/home/sobczm/bin/popgen/renseq
input=/home/sobczm/popgen/renseq/strawberry/reads/Helen_Bates_EMR.RH.ENQ-1704.A.01/analysis
cd $input
#Following the running of BLAST on RG and Hapil HQ CCS and HGAP assembly against F. vesca NBS-LRRs CDS, we do a reciprocal blast, and carry out a RBB analysis. 
for db in S1_ccs_3_99_nucl.db S2_ccs_3_99_nucl.db S1_HGAP_polished_assembly_nucl.db S2_HGAP_polished_assembly_nucl.db
do
blastn -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand"  -num_threads 1 -max_target_seqs 1000000000 -evalue 0.0000000001 -query vesca_v1.1_nblrrs_augustus_cds.fasta -db $db >> vesca_v1.1_nblrrs_augustus_cds.fasta_vs_$db
done

#Identify RBB matches.
for sample in S1_ccs_3_99 S2_ccs_3_99 S1_HGAP_polished_assembly S2_HGAP_polished_assembly
do
python /home/sobczm/bin/popgen/other/rbb.py ${sample}.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db vesca_v1.1_nblrrs_augustus_cds.fasta_vs_${sample}_nucl.db >${sample}_vesca_v1.1_nblrrs_augustus_cds_rbb
done

#Count up the length of CDS, 5' UTR, 3' UTR using the original renseq vs vesca CDS BLAST output
for sample in S1_ccs_3_99 S2_ccs_3_99 S1_HGAP_polished_assembly S2_HGAP_polished_assembly
do
python $scripts/blast_matches_summary2.py vesca_v1.1_nblrrs_augustus_cds_lengths.txt ${sample}_lengths.txt ${sample}.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db > ${sample}.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_cds
done
#Filter those results to only retain the RBB hits.
for sample in S1_ccs_3_99 S2_ccs_3_99 S1_HGAP_polished_assembly S2_HGAP_polished_assembly
do
python $scripts/filter_by_two_columns.py ${sample}_vesca_v1.1_nblrrs_augustus_cds_rbb ${sample}.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_cds > ${sample}.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_cds_rbb_filt
python $scripts/filter_by_two_columns.py ${sample}_vesca_v1.1_nblrrs_augustus_cds_rbb ${sample}.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_5utr > ${sample}.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_5utr_rbb_filt 
python $scripts/filter_by_two_columns.py ${sample}_vesca_v1.1_nblrrs_augustus_cds_rbb ${sample}.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_3utr > ${sample}.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_3utr_rbb_filt 
done