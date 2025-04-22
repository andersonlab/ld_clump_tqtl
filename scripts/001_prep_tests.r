# Bradley March 2025
library(tidyverse)

# Get paths
repo.dir <- '/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/code/IBDVerse-sc-eQTL-code/'
sumstats.all.basedir <- '/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_qtlight/2025_03_20-multi_tissue_base_results/TensorQTL_eQTLS/'
sumstats.interaction.basedir <- '/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_qtlight/2025_03_26-multi_tissue_interaction_results/TensorQTL_eQTLS/'
out.dir <- 'input'
data.dir <- paste0(repo.dir,'/data/')
coloc.dir <- "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/snakemake_colocalisation/results/2025_03_IBDverse_coloc_all_gwas/collapsed/"
int.coloc.dir <- "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/results/coloc/2025_03_IBDverse_coloc_interaction_all_gwas/collapsed"
source(paste0(repo.dir,'qtl_plot/helper_functions.R'))

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

# Get sig gene x condition pairs
sig.gene.condition = sig.sumstat.df %>% 
    mutate(pheno_annotation = paste0(phenotype_id, "-", annotation)) %>% 
    pull(pheno_annotation) %>% 
    unique()

length(sig.gene.condition) # 330504

##################
# Read in conditional eQTLs and prep
##################
cond.sumstat.df <- read_conditional_eqtls(sumstats.all.basedir)
cond.sumstat.df <- cond.sumstat.df %>% 
    mutate(pheno_annotation = paste0(phenotype_id, "-", annotation))

# get a unique list of phenotype IDs and variants
qtl_gene_leads = cond.sumstat.df %>% 
    filter(pheno_annotation %in% sig.gene.condition) %>% 
    group_by(variant_id, phenotype_id) %>% 
    summarise(
        P=min(pval_nominal) # Not interested in whether the same gene-snp pair is more/less significant in different conditions
    ) %>% 
    mutate(
        type="eQTL"
    ) # 196657

# Load in interaction eQTLs
interactions = read_ieqtls(sumstats.interaction.basedir) %>% 
    filter(
        pval_adj_bh < 0.05
    ) %>% 
    mutate(
        variant_id = factor(variant_id, levels=unique(variant_id)),
        interaction_new = interaction.mapping[interaction]
    )

int_gene_leads = interactions %>%
    group_by(variant_id, phenotype_id) %>%
    filter(pval_adj_bh < 0.05) %>%  
    summarise(
        P=min(pval_adj_bh)
    ) %>% 
    mutate(
        type="interaction_eQTL"
    )

# Add this
qtl_gene_leads = qtl_gene_leads %>% bind_rows(int_gene_leads) %>% 
    group_by(variant_id, phenotype_id) %>% 
    slice_min(P) %>% 
    distinct() # 197,770


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
    ) # 10,051

# Interactions: 
int_colocs_all = get_colocs(int.coloc.dir, known_ibd_only = F) %>% 
    filter(
        !(gwas_trait %in% c("ibd", "cd", "uc")),
        PP.H4.abf > 0.75
    )

int_coloc_gene_leads = int_colocs_all %>% 
    mutate(variant_id = gsub("\\_", "\\:", snp_id)) %>% 
    group_by(variant_id, phenotype_id) %>% 
    summarise(
        P=min(qtl_pval) # Not interested in whether the same gene-snp pair is more/less significant in different conditions
    ) %>% 
    mutate(
        type="interaction_coloc"
    ) # 16

# Combine
coloc_gene_leads = coloc_gene_leads %>% bind_rows(int_coloc_gene_leads) %>% 
    group_by(variant_id, phenotype_id) %>% 
    slice_min(P) %>% 
    distinct() # 10,067

##################
# Combine these together, saveas well as a summary of the genes
##################
coqtl = coloc_gene_leads %>% 
    bind_rows(qtl_gene_leads) %>% 
    group_by(variant_id, phenotype_id) %>% 
    slice_min(P) %>% 
    distinct() %>% 
    arrange(phenotype_id, variant_id, P) %>% 
    rename(SNP = variant_id) %>% 
    select(phenotype_id, SNP, P, type)

nrow(coqtl) # 198045

#            coloc              eQTL interaction_coloc  interaction_eQTL 
#             4784            192148                16              1097

# summary_genes
write.table(coqtl %>% pull(phenotype_id) %>% unique(), paste0(out.dir, "/gene_list.txt"), col.names=F, row.names=F, quote=F)
# variants per gene
write.table(coqtl, paste0(out.dir, "/variants_per_gene_list.txt"), col.names=T, row.names=F, quote=F, sep = "\t")