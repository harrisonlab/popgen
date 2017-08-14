The order in which main bash scripts in the folder are executed:

1) [pre_SNP_calling_cleanup.sh](https://github.com/harrisonlab/popgen/blob/master/snp/pre_SNP_calling_cleanup.sh)

2) [fus_SNP_calling_multithreaded](https://github.com/harrisonlab/popgen/blob/master/snp/fus_SNP_calling_multithreaded.sh) (with fus samples hardcoded)

3) [structural_variants.sh](https://github.com/harrisonlab/popgen/blob/master/snp/structural_variants.sh)

4) [determine_genetic_structure.sh](https://github.com/harrisonlab/popgen/blob/master/snp/determine_genetic_structure.sh)

5) [structure_analysis.sh](https://github.com/harrisonlab/popgen/blob/master/snp/structure_analysis.sh)

6) [fast_structure_analysis.sh](https://github.com/harrisonlab/popgen/blob/master/snp/fast_structure_analysis.sh)

