# Bradley April 2025
mkdir -p ld_between_conditional_temp
eqtldir=/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_qtlight/2025_03_20-multi_tissue_base_results/TensorQTL_eQTLS/

# Make log dir:
mkdir -p logs/check_ld_conditional

# For each condition, get the conditional results and look at the LD between the conditional variants
for f in $eqtldir*; do 
    condition=$(basename $f)
    echo $condition
    bsub -M 1000 -a "memlimit=True" -R "select[mem>1000] rusage[mem=1000] span[hosts=1]" \
        -o logs/check_ld_conditional/${condition}-%J.out \
        -e logs/check_ld_conditional/${condition}-%J.err \
        -q normal \
        -J "check_LD_${condition}" \
        "bash scripts/000_check_ld_conditional_submit.sh ${eqtldir} ${condition} ld_between_conditional_temp"
done
