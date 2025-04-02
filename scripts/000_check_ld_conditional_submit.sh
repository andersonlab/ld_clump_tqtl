# Bradley April 2025
#BSUB -o logs/ld_conditional-%J-output.log
#BSUB -e logs/ld_conditional-%J-error.log 
#BSUB -q oversubscribed
#BSUB -G team152
#BSUB -n 1
#BSUB -M 1000
#BSUB -a "memlimit=True"
#BSUB -R "select[mem>1000] rusage[mem=1000] span[hosts=1]"
#BSUB -J 1

# get options
basedir="$1"
echo "basedir: $basedir" # basedir=/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_qtlight/2025_03_20-multi_tissue_base_results/TensorQTL_eQTLS/
condition="$2"
echo "condition: $condition" # condition=dMean__unannotated_ti_all
pathout="$3"
echo "pathout: $pathout" # pathout=ld_between_conditional_temp

# Get a list of unique, significant genes from the qvalue results
qval_f="${basedir}${condition}/OPTIM_pcs/base_output/base/Cis_eqtls_qval.tsv"
if [ -f "$qval_f" ]; then # Only run rest if this exists
    echo "..Getting sig egenes"
    awk -F'\t' '$18 < 0.05' ${qval_f} | awk '{print $1}' | uniq > ${pathout}/${condition}_sigegenes.txt

    # For each gene, make a file for the unique variants
    echo "..Looping through conditional effects"
    conditional_f="${basedir}${condition}/OPTIM_pcs/base_output/base/Cis_eqtls_independent.tsv"
    while read g; do
        if [ $(grep "$g" "$conditional_f" | wc -l) -gt 1 ]; then # Only run if >1 conditional variants
            echo $g
            grep "$g" "$conditional_f" | awk '{print $7}' > "${pathout}/${condition}_${g}_cond_vars.txt"
            num_chr=$(awk -F':' 'NR==1 {{gsub("chr", "", $1); print $1}}' "${pathout}/${condition}_${g}_cond_vars.txt")
            plink --bfile /lustre/scratch126/humgen/projects/sc-eqtl-ibd/data/genotypes/eQTL_genotypes_march_2024/plink_march_2025/plink_geno_chr${num_chr}_plink --extract "${pathout}/${condition}_${g}_cond_vars.txt" --make-bed --out "${pathout}/${condition}_${g}"
            # Compute LD of these
            plink --bfile ${pathout}/${condition}_${g} --ld-snp-list "${pathout}/${condition}_${g}_cond_vars.txt" --ld-window-kb 50000000 --ld-window-r2 0 --r2 --out "${pathout}/${condition}_${g}"
            # Add to out file
            awk -v cond="$condition" -v gene="$g" 'NR==1 {print $0, "condition", "gene"; next} {print $0, cond, gene}' \
                "${pathout}/${condition}_${g}.ld" | tail -n +2 >> "${pathout}/LD_${condition}_independent.txt"
            # Remove intermediate files
            rm ${pathout}/${condition}_${g}*
        fi
    done <${pathout}/${condition}_sigegenes.txt

    # If there is any output, compress:
    if [ -f "${pathout}/LD_${condition}_independent.txt" ]; then 
        # Compress the output
        echo "..zipping"
        gzip ${pathout}/LD_${condition}_independent.txt
    else
        echo "..No secondary or higher effects!"
    fi

    # Remove the sigegene list
    rm ${pathout}/${condition}_sigegenes.txt
else
    echo "..No eQTLs!"
fi