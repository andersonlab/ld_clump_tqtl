# Bradley March 2025
library(ggrepel)
library(shadowtext)
library(patchwork)
library(tidyverse)
library(forcats)
library(ggpubr)

# Get paths
repo.dir <- '/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/code/IBDVerse-sc-eQTL-code/'
sumstats.all.basedir <- '/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_qtlight/2024_12_31-multi_tissue_base_results/TensorQTL_eQTLS/'
out.dir <- 'input/per_gene'
data.dir <- paste0(repo.dir,'/data/')
coloc.dir <- "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/snakemake_colocalisation/results/2025_03_IBDverse_coloc_all_gwas/collapsed/"
source(paste0(repo.dir,'qtl_plot/helper_functions.R'))
source(paste0('qtl_plot/pleitropy_helpers.r'))

# Make outdirs
if(!file.exists(out.dir)){
    dir.create(out.dir, recursive=T)
}

##################
# Read in eQTLs and prep
##################
sumstat.df <- read_eqtls(sumstats.all.basedir)
sig.sumstat.df <- sumstat.df %>% 
  filter(qval < 0.05) 

# get a unique list of phenotype IDs and variants
qtl_gene_leads = sig.sumstat.df %>% 
    group_by(variant_id, phenotype_id) %>% 
    summarise(
        P=min(pval_nominal) # Not interested in whether the same gene-snp pair is more/less significant in different conditions
    ) %>% 
    mutate(
        type="eQTL"
    )

##################
# Read in colocs and prep
##################

colocs_all = get_colocs(coloc.dir, known_ibd_only = F) %>% 
    filter(
        !(gwas_trait %in% c("ibd", "cd", "uc")),
        PP.H4.abf > 0.75
    )

# get a unique list of phenotype IDs and variants
coloc_gene_leads = colocs_all %>% 
    mutate(variant_id = gsub("\\_", "\\:", snp_id)) %>% 
    group_by(variant_id, phenotype_id) %>% 
    summarise(
        P=min(qtl_pval) # Not interested in whether the same gene-snp pair is more/less significant in different conditions
    ) %>% 
    mutate(
        type="coloc"
    )


##################
# Combine these together, save a per gene file as well as a summary of the genes
##################
coqtl = coloc_gene_leads %>% 
    bind_rows(qtl_gene_leads) %>% 
    arrange(phenotype_id, variant_id, P) %>% 
    rename(SNP = variant_id) %>% 
    select(phenotype_id, SNP, P)

# summary
write.table(coqtl %>% pull(phenotype_id) %>% unique(), paste0(out.dir, "/gene_list.txt"), col.names=F, row.names=F, quote=F)

for(g in unique(coqtl$phenotype_id)){
    write.table(coqtl %>% 
        filter(phenotype_id == g)
    , paste0(out.dir, "/", g, "_variants_to_clump.txt"), row.names=F, quote=F, sep = "\t")
}


