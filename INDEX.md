##Index to "popgen" analyses available to use
###**Please make sure that the program you are trying to use has its dependencies in the PATH. See:**
###https://github.com/harrisonlab/popgen/blob/master/programs.md
My current bash profile
`/home/sobczm/bin/maria_bash_profile`

*What directories are there?*

1. Phylogenetics
2. Snp
3. Summary stats
4. Clock
5. Codon

###*Each directory contains a README.md file listing shell scripts which contain a model (example) analysis using scripts in a given directory.
###*Read the header of each individual script you are trying to exectute to find out about the options, input and output file.

*What is in each directory? Only listing re-usable analyses*

1. **Phylogenetics:** Making Bayesian phylogenetic trees from gene sequences
2. **Snp:** Mapping reads to reference, cleaning the mapped files and calling SNPs. Filtering SNPs and their basic summary stats. Various methods of establishing population structure, including SNP-based NJ trees.
3. **Summary stats:** More complex analyses involving SNP sampling across populations identified in 2. Divergence, natural selection, recombination.
4. **Clock:** De novo motif discovery and motif scanning. Orthology analysis and automated generation of ortholog trees. Pairwise, branch-site and branch dN/dS codeml models.
5. **Codon:** Analysis of gene codon usage and gene duplication levels.

##Phylogenetics
###Scripts to make a Bayesian phylogenetic tree based on gene sequences.
**Model analysis file:** [BUSCO_analysis.sh] (https://github.com/harrisonlab/popgen/blob/master/phylogenetics/BUSCO_analysis.sh)

[sub_BUSCO2.sh] (https://github.com/harrisonlab/popgen/blob/master/phylogenetics/sub_BUSCO2.sh) 
Establish the number of single copy Fungal/Eukaryotic/Plant/Prokaryotic/ conserved genes in a genome or transcriptome. 

[get_alignments.pl] (https://github.com/harrisonlab/popgen/blob/master/phylogenetics/get_alignments.pl) 
This is the script to use if you want to follow the steps in the model analysis file and create individual alignments for each BUSCO gene found to be complete in each genome analysed. These gene alignments (only CDS recommended) can then be used to make trees in the next step.

**Model analysis file:** [pre_BEAST_prep.sh] (https://github.com/harrisonlab/popgen/blob/master/phylogenetics/pre_BEAST_prep.sh)

[sub_mafft_alignment.sh] (https://github.com/harrisonlab/popgen/blob/master/phylogenetics/sub_mafft_alignment.sh) 
Quickly align the sequences generated in 1. to be used for tree construction. The script loops over all the FASTA files in the directory and tries to align sequences inside each one.

[calculate_nucleotide_diversity.py] (https://github.com/harrisonlab/popgen/blob/master/phylogenetics/calculate_nucleotide_diversity.py)
calculates basic sequence diversity stats for each alignment which could then guide selection of individual genes towards making the
phylogenetic tree (see the model analysis file for a example of an analysis)

*The remainder of the model analysis file describes establishing the correct model of evolution for each gene using PartitionFinder that can be then implemented within BEAST. With some prior information on that, this step can be potentially abandoned.*

**Model analysis file:** [BEAST_run.sh] (https://github.com/harrisonlab/popgen/blob/master/phylogenetics/BEAST_run.sh)
The BEAST analysis has to be so far set-up by hand. A guide to do so and obtain a final tree is given in the model analysis file above.

##SNP
###Scripts to call SNPs on multiple individuals using a single genome/transcriptome reference, filter (and downsample) them, and establish the basic population structure.

**Model analysis file:** [pre_SNP_calling_cleanup.sh] (https://github.com/harrisonlab/popgen/blob/master/snp/pre_SNP_calling_cleanup.sh)

[sub_pre_snp_calling.sh] (https://github.com/harrisonlab/popgen/blob/master/snp/sub_pre_snp_calling.sh)
Script accepts SAM mappings output by Bowtie2 along with sample ID, and outputs filtered, indexed and ID-tagged BAM files to be used for variant calling

**Model analysis file:** [fus_SNP_calling_multithreaded.sh] (https://github.com/harrisonlab/popgen/blob/master/snp/fus_SNP_calling_multithreaded.sh)
This script allows variant calling with GATK. It needs to be modified for each GATK run. First, it prepares genome reference indexes required by GATK and then calls the variant with the GATK package in the custom script similar to [sub_fus_SNP_calling_multithreaded.sh]
(https://github.com/harrisonlab/popgen/blob/master/snp/sub_fus_SNP_calling_multithreaded.sh). 

**Model analysis file:** [determine_genetic_structure.sh] (https://github.com/harrisonlab/popgen/blob/master/snp/determine_genetic_structure.sh)
The script contains a collection of scripts carrying out the following SNP-based population structure analyses:

1. Variant stats: 
`/home/sobczm/bin/vcftools/bin/vcf-stats` 
2. Default variant filtering on the input VCF file recommended before structure analysis and required by some tools:
[sub_vcf_parser.sh] (https://github.com/harrisonlab/popgen/blob/master/snp/sub_vcf_parser.sh) 
3. Subsampling of SNPs down to 1 SNP per x bp on a contig with `/home/sobczm/bin/vcftools/bin/vcftools`
4. Basic handle on the partitioning of SNPs between individuals. Calculate percentage of shared SNP alleles between each sample in
a pairwise manner and k-mean cluster the samples, and generate a heatmap and a dendrogram as the output:
[similarity_percentage.py] (https://github.com/harrisonlab/popgen/blob/master/snp/similarity_percentage.py)
[distance_matrix.R] (https://github.com/harrisonlab/popgen/blob/master/snp/distance_matrix.R) 
5. Another way to detect any population genetic structure between samples: PCA
[pca.R] (https://github.com/harrisonlab/popgen/blob/master/snp/pca.R)
6. Show relationships between samples using a neighbour-joining tree:
[nj_tree.sh] (https://github.com/harrisonlab/popgen/blob/master/snp/nj_tree.sh)
7. Carry out a custom AMOVA analysis to try to partition genetic variation between a hypothetical factor, such as virulence level or geographic origin (requires considerable adaptation to individual analysis):
[amova_dapc.R] (https://github.com/harrisonlab/popgen/blob/master/snp/amova_dapc.R)

**Model analysis file:** [structure_analysis.sh]
(https://github.com/harrisonlab/popgen/blob/master/snp/structure_analysis.sh)

Downsample the VCF file with SNPs prior to analysis with the STRUCTURE program.
`/home/sobczm/bin/vcflib/bin/vcfrandomsample --rate 0.1 $input_vcf >$output_vcf`

Run the STRUCTURE analysis to test for thelikely number of population clusters (K) (can be in the range of: K=1 up to K=number of individuals tested), summarise the results with StructureHarvester and CLUMPP, visualise with DISTRUCT. 

##Summary stats 
###Scripts for functional annotation of SNPs, and calculation of general population genetics parameters (Fst, nuclotide diversity, Tajima's D) which can be informative about demographic and selection processes operating on a given gene(s) in tested populations. Analyses available include both haplotype- and nucleotide- based.
 **Model analysis file:** [fus_variant_annotation.sh] (https://github.com/harrisonlab/popgen/blob/master/summary_stats/fus_variant_annotation.sh)

1. Create VCF file subsets only with certain individuals retained for downstream analysis with `/home/sobczm/bin/vcflib/bin/vcfremovesamples` and remove monomorphic (idential) positions in the output with `/home/sobczm/bin/vcftools/bin/vcftools`

2. Create custom SNPEff annotation for each new genome which allows classification of variants into different functional categories
First, build genome database with [build_genome_database.sh] (https://github.com/harrisonlab/popgen/blob/master/summary_stats/build_genome_database.sh) Secondly, annotate the variants in select VCF file and create subsets of SNPs with different effect (genic, coding, synonymous, non-synonymous, 4-fold degenerate) with [annotate_snps_genome.sh] (https://github.com/harrisonlab/popgen/blob/master/summary_stats/annotate_snps_genome.sh)

3. Convert a VCF file to independent FASTA sequences for each individual (sample), incorporating the polymorphisms stored in the VCF file - that is, generate a custom version of each genome based on the common reference and SNPs detected during variant calling [vcf_to_fasta.py] (https://github.com/harrisonlab/popgen/blob/master/summary_stats/vcf_to_fasta.py)

4. Split a GFF file with gene annotation for the reference genome into independent file per contig. Required by PopGenome used in the next script. [split_gff_contigs.sh] (https://github.com/harrisonlab/popgen/blob/master/summary_stats/split_gff_contig.sh)

**Model analysis file:** [fus_popgenome_analysis.sh] (https://github.com/harrisonlab/popgen/blob/master/summary_stats/fus_popgenome_analysis.sh)
The scripts used in this part of the analysis described in the model analysis file above require a special FASTA and GFF input generated with scripts [vcf_to_fasta.py] (https://github.com/harrisonlab/popgen/blob/master/summary_stats/vcf_to_fasta.py) and [split_gff_contigs.sh] (https://github.com/harrisonlab/popgen/blob/master/summary_stats/split_gff_contig.sh). The input needs to be also arranged in a particular way described in the model analysis file above. Lastly, the submitted R scripts themselves require customisation for each analysis regarding sample names and their population assingment.

1. [sub_calculate_nucleotide_diversity.sh] (https://github.com/harrisonlab/popgen/blob/master/summary_stats/sub_calculate_nucleotide_diversity.sh) Calculate, plot (histogram and line plots) and output in tabular format (available per contig and per entire genome) the following statistics calculated over a gene and sliding window-based intervals (two types of plots: individual populations and multiple specified populations on one plot): nucleotide diversity within populations - Pi (Nei, 1987), ratio of Pi (nonsynonymous/synonymous), pairwise nucleotide diversity between populations - Dxy. 

2. [sub_calculate_neutrality_stats.sh] (https://github.com/harrisonlab/popgen/blob/master/summary_stats/sub_calculate_neutrality_stats.sh)
Calculate, plot (histogram and line plots) and output in tabular format (available per contig and per entire genome) the following statistics calculated over a gene and sliding window-based intervals (two types of plots: individual populations and multiple specified populations on one plot): Tajima's D, and related (Fu & Li's F, Fu & Li's D), Watterson's Theta per site, Average rate of segregating sites.

3. [sub_calculate_fst.sh] (https://github.com/harrisonlab/popgen/blob/master/summary_stats/sub_calculate_fst.sh)
Calculate, plot (histogram and line plots) and output in tabular format (available per contig and per entire genome) the following statistics calculated over a gene and sliding window-based intervals (two types of plots: individual populations and multiple specified populations on one plot): Weir and Cockerham's Fst (across all populations), Weir and Cockerham's Fst (pairwise between populations), Hudson's Kst.

4. [sub_calculate_haplotype_based_stats.sh] (https://github.com/harrisonlab/popgen/blob/master/summary_stats/sub_calculate_haplotype_based_stats.sh)
Calculate EXACTLY the same diversity statistiscs, as in 2. [sub_calculate_neutrality_stats.sh] (https://github.com/harrisonlab/popgen/blob/master/summary_stats/sub_calculate_neutrality_stats.sh) and [sub_calculate_fst.sh] (https://github.com/harrisonlab/popgen/blob/master/summary_stats/sub_calculate_fst.sh) BUT using HAPLOTYPE not NUCLEOTIDE sequences as input. In addition, carry out a four gamete test on each population to check for presence of recombination.
 
**Model analysis file:** [fus_linkage_disequilibrum.sh] (https://github.com/harrisonlab/popgen/blob/master/summary_stats/fus_linkage_disequilibrum.sh)
Use `/home/sobczm/bin/vcftools/bin/vcftools` to calculate D, D' and r2 for SNPs seperated by a specific range of intervals to estimate recombination rates and subsequently visualise the results (D' and r2 versus physical distance, histogram of D' values) using [sub_plot_ld.sh] (https://github.com/harrisonlab/popgen/blob/master/summary_stats/sub_plot_ld.sh). 

