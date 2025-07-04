# Code used to clump eQTLs across many different conditions

### Installation

If using Sanger farm. No downloads or installs are required. 
Otherwise, singularity image can be downloaded from dockerhub:
```
export SINGULARITY_CACHEDIR=$PWD/.singularity_cache
mkdir -p "$SINGULARITY_CACHEDIR"
singularity pull docker://bh18/ld_clump:4
```
- Then adjust the singularity paths in `scripts/run_ldclump.smk` to the location of the downloaded image.

## 1. Preperation of tests
- Before running the analysis, need to prep tests and populate the input files. 
- These include a list of genes to test (as eQTLs will be grouped within these), and variant-gene maps
- An example script for doing this from various TensorQTL/coloc outputs is detailed in `scripts/001_prep_tests.r`

## 2. Run pipeline
- To execute (on LSF at Sanger):
```bsub -M 10000 -a "memlimit=True" -R "select[mem>10000] rusage[mem=10000] span[hosts=1]" -o sm_logs/snakemake_master-%J-output.log -e sm_logs/snakemake_master-%J-error.log -q oversubscribed -J "snakemake_master_CLUMP" < submit_snakemake_BH.sh ```
- If non-Sanger, this will need adjusting for specific singularity/snakemake versions (we use 7)
