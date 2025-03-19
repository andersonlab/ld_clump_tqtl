# Bradley March 2025
library(dplyr)
library(tidyr)
args = commandArgs(trailingOnly=TRUE)
infile = args[1] # infile="temp/ENSG00000000457.clumped"
outfile = args[2] # outfile="results/per_gene/ENSG00000000457_clumped_format.txt"

# Load in the infile
res = read.delim(infile, sep="")

# Extract gene
gene = gsub(".clumped", "", unlist(strsplit(infile, "\\/"))[2])

# Transform
res_transformed <- res %>%
  mutate(SP2 = gsub("\\(1\\)", "", SP2)) %>%  # Remove "(1)"
  mutate(SP2 = ifelse(SP2 == "NONE", "", SP2)) %>%  # Replace "NONE" with empty string
  separate_rows(SP2, sep = ",") %>%  # Split multiple values into separate rows
  rename(qtl_clump_index = SNP, variant_id = SP2) %>%
  select(qtl_clump_index, variant_id) %>%
  bind_rows(res %>% select(qtl_clump_index = SNP) %>% mutate(variant_id = qtl_clump_index)) %>%
  arrange(qtl_clump_index, variant_id) %>% 
  filter(
    variant_id != ""
  ) %>% 
  mutate(phenotype_id = gene) %>% 
  select(phenotype_id, variant_id, qtl_clump_index)

# Save
write.table(res_transformed, outfile, sep = "\t", quote=F, row.names=F, col.names=F)
