#!/bin/bash
scripts=/home/sobczm/bin/popgen/renseq
input=/home/sobczm/popgen/renseq/strawberry/reads/Helen_Bates_EMR.RH.ENQ-1704.A.01
cd $input/assembly
#Final functional analysis  of HGAP, Falcon assemblies and CCS (HQ: 99% precision) reads from strawberry renseq

#NBS Parser
for a in S1_ccs_3_99.fasta S2_ccs_3_99.fasta S2_HGAP_polished_assembly.fasta S1_HGAP_polished_assembly.fasta
do
cp $a ../analysis/nbs-parser
cd ../analysis/nbs-parser
sh $scripts/sub_nlrparser.sh $(basename $a)
cd ../../assembly
done

#Create blast databases out of the assemblies
for assembly in S1_ccs_3_99.fasta S2_ccs_3_99.fasta S2_HGAP_polished_assembly.fasta S1_HGAP_polished_assembly.fasta
do
makeblastdb -in $assembly -input_type fasta -dbtype nucl \
-title "${assembly%.*}"_nucl.db -parse_seqids -out "${assembly%.*}"_nucl.db
done

#Search the assemblies for Rpf2
for db in S1_ccs_3_99_nucl.db S2_ccs_3_99_nucl.db S1_HGAP_polished_assembly_nucl.db S2_HGAP_polished_assembly_nucl.db
do
blastn -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand"  -num_threads 1 -max_target_seqs 100 -evalue 0.0000000001 -query F.vesca_Rpf2.fasta -db $db >> F.vesca_Rpf2_vs_$db
done

#Run blast search against the CDS database of chosen genes as well as the entire strawberry genome CDS database as well as baits database.
cd $input/analysis
for db in probes-R4-final.nucl.db Fragaria_vesca_v1.1.a2_cds_removed.fasta_nucl.db vesca_v1.1_nblrrs_augustus_cds_nucl.db
do
for ass in S1_ccs_3_99.fasta S2_ccs_3_99.fasta S2_HGAP_polished_assembly.fasta S1_HGAP_polished_assembly.fasta
do
    blastn -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand"  -num_threads 1 -max_target_seqs 1000000000 -evalue 0.0000000001 -query $ass -db $db >> $(basename $ass)_vs_$db
done
done
#Search the baits sequences against the databases of sequences of Ren-Seq reads/assemblies.
for db in F06.assembly_nucl.db G06.assembly_nucl.db S1_ccs_3_99_nucl.db S2_ccs_3_99_nucl.db S1_HGAP_polished_assembly_nucl.db S2_HGAP_polished_assembly_nucl.db
do
blastn -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand"  -num_threads 1 -max_target_seqs 1000000000 -evalue 0.0000000001 -query probes-R4-final.fas -db $db >> probes-R4-final.fas_vs_$db
done


#Generate a file with sequence lengths for CDS input files.
#Note for future reference: BLAST can automatically output these values if fields slen and qlen specified.
for a in S1_ccs_3_99.fasta S2_ccs_3_99.fasta S2_HGAP_polished_assembly.fasta S1_HGAP_polished_assembly.fasta
do
python $scripts/write_seq_length.py $a
done

