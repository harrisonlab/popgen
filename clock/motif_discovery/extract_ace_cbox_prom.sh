#!/bin/bash
input=/home/sobczm/popgen/clock/DNA_genomes/promoters/extended
scripts=/home/sobczm/bin/popgen/clock/motif_discovery

cd $input
python $scripts/keep_list_genes.py $input/ace/Botrytis_cinerea.pep.fa.ace \
Botrytis_cinerea.ASM83294v1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/ace/Chaetomium_globosum.pep.fa.ace \
Chaetomium_globosum_cbs_148_51.ASM14336v1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/ace/Chaetomium_thermophilum.pep.fa.ace \
Chaetomium_thermophilum_var_thermophilum_dsm_1495.CTHT_3.0.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/ace/Fus2.pep.fa.ace \
Fus2_canu_contigs_hardmasked_upstream2000.fa
python $scripts/keep_list_genes.py $input/ace/Fusarium_fujikuroi.pep.fa.ace \
Fusarium_fujikuroi.EF1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/ace/Fusarium_graminearum.pep.fa.ace \
Fusarium_graminearum.RR.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/ace/Fusarium_langsethiae.pep.fa.ace \
Fusarium_langsethiae.ASM129263v1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/ace/Fusarium_pseudograminearum.pep.fa.ace \
Fusarium_pseudograminearum.GCA_000303195.1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/ace/Fusarium_solani.pep.fa.ace \
Fusarium_solani.v2.0.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/ace/Fusarium_verticillioides.pep.fa.ace \
Fusarium_verticillioides.ASM14955v1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/ace/Gaeumannomyces_graminis.pep.fa.ace \
Gaeumannomyces_graminis.Gae_graminis_V2.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/ace/Magnaporthe_oryzae.pep.fa.ace \
Magnaporthe_oryzae.MG8.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/ace/Magnaporthe_poae.pep.fa.ace \
Magnaporthe_poae.Mag_poae_ATCC_64411_V1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/ace/Neonectria_ditissima.pep.fa.ace \
Neonectria_ditissima.R0905_v2.0.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/ace/Neurospora_crassa.pep.fa.ace \
Neurospora_crassa.NC12.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/ace/Neurospora_tetrasperma.pep.fa.ace \
Neurospora_tetrasperma_fgsc_2508.v2.0.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/ace/Podospora_anserina.pep.fa.ace \
Podospora_anserina_s_mat_.ASM22654v1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/ace/Sclerotinia_borealis.pep.fa.ace \
Sclerotinia_borealis_f_4157.SBOR_1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/ace/Sclerotinia_sclerotiorum.pep.fa.ace \
Sclerotinia_sclerotiorum.ASM14694v1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/ace/Sordaria_macrospora.pep.fa.ace \
Sordaria_macrospora.ASM18280v2.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/ace/Thielavia_terrestris.pep.fa.ace \
Thielavia_terrestris_nrrl_8126.ASM22611v1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/ace/Trichoderma_atroviride.pep.fa.ace \
Trichoderma_atroviride_imi_206040.TRIAT_v2_0.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/ace/Trichoderma_gamsii.pep.fa.ace \
Trichoderma_gamsii.ASM148177v1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/ace/Trichoderma_harzianum.pep.fa.ace \
Trichoderma_harzianum.ASM98886v1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/ace/Trichoderma_reesei.pep.fa.ace \
Trichoderma_reesei.GCA_000167675.2.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/ace/Trichoderma_virens.pep.fa.ace \
Trichoderma_virens.ASM17099v1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/ace/Verticillium_alfalfae.pep.fa.ace \
Verticillium_alfalfae_vams_102.ASM15082v1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/ace/Verticillium_dahliae.pep.fa.ace \
Verticillium_dahliae.ASM15067v2.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/ace/Verticillium_longisporum.pep.fa.ace \
Verticillium_longisporum_gca_001268165.vl2.denovo.v1.dna_rm.toplevel_promoters_2000.fasta
mv $input/*filtered* $input/ace

python $scripts/keep_list_genes.py $input/cbox/Botrytis_cinerea.pep.fa.cbox \
Botrytis_cinerea.ASM83294v1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/cbox/Chaetomium_globosum.pep.fa.cbox \
Chaetomium_globosum_cbs_148_51.ASM14336v1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/cbox/Chaetomium_thermophilum.pep.fa.cbox \
Chaetomium_thermophilum_var_thermophilum_dsm_1495.CTHT_3.0.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/cbox/Fus2.pep.fa.cbox \
Fus2_canu_contigs_hardmasked_upstream2000.fa
python $scripts/keep_list_genes.py $input/cbox/Fusarium_fujikuroi.pep.fa.cbox \
Fusarium_fujikuroi.EF1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/cbox/Fusarium_graminearum.pep.fa.cbox \
Fusarium_graminearum.RR.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/cbox/Fusarium_langsethiae.pep.fa.cbox \
Fusarium_langsethiae.ASM129263v1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/cbox/Fusarium_pseudograminearum.pep.fa.cbox \
Fusarium_pseudograminearum.GCA_000303195.1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/cbox/Fusarium_solani.pep.fa.cbox \
Fusarium_solani.v2.0.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/cbox/Fusarium_verticillioides.pep.fa.cbox \
Fusarium_verticillioides.ASM14955v1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/cbox/Gaeumannomyces_graminis.pep.fa.cbox \
Gaeumannomyces_graminis.Gae_graminis_V2.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/cbox/Magnaporthe_oryzae.pep.fa.cbox \
Magnaporthe_oryzae.MG8.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/cbox/Magnaporthe_poae.pep.fa.cbox \
Magnaporthe_poae.Mag_poae_ATCC_64411_V1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/cbox/Neonectria_ditissima.pep.fa.cbox \
Neonectria_ditissima.R0905_v2.0.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/cbox/Neurospora_crassa.pep.fa.cbox \
Neurospora_crassa.NC12.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/cbox/Neurospora_tetrasperma.pep.fa.cbox \
Neurospora_tetrasperma_fgsc_2508.v2.0.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/cbox/Podospora_anserina.pep.fa.cbox \
Podospora_anserina_s_mat_.ASM22654v1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/cbox/Sclerotinia_borealis.pep.fa.cbox \
Sclerotinia_borealis_f_4157.SBOR_1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/cbox/Sclerotinia_sclerotiorum.pep.fa.cbox \
Sclerotinia_sclerotiorum.ASM14694v1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/cbox/Sordaria_macrospora.pep.fa.cbox \
Sordaria_macrospora.ASM18280v2.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/cbox/Thielavia_terrestris.pep.fa.cbox \
Thielavia_terrestris_nrrl_8126.ASM22611v1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/cbox/Trichoderma_atroviride.pep.fa.cbox \
Trichoderma_atroviride_imi_206040.TRIAT_v2_0.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/cbox/Trichoderma_gamsii.pep.fa.cbox \
Trichoderma_gamsii.ASM148177v1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/cbox/Trichoderma_harzianum.pep.fa.cbox \
Trichoderma_harzianum.ASM98886v1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/cbox/Trichoderma_reesei.pep.fa.cbox \
Trichoderma_reesei.GCA_000167675.2.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/cbox/Trichoderma_virens.pep.fa.cbox \
Trichoderma_virens.ASM17099v1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/cbox/Verticillium_alfalfae.pep.fa.cbox \
Verticillium_alfalfae_vams_102.ASM15082v1.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/cbox/Verticillium_dahliae.pep.fa.cbox \
Verticillium_dahliae.ASM15067v2.dna_rm.toplevel_promoters_2000.fasta
python $scripts/keep_list_genes.py $input/cbox/Verticillium_longisporum.pep.fa.cbox \
Verticillium_longisporum_gca_001268165.vl2.denovo.v1.dna_rm.toplevel_promoters_2000.fasta
mv $input/*filtered* $input/cbox
