# Bradley April 2025
# Checking the LD between indepdent variants

library(ggplot2)
library(dplyr)

data.dir <- 'ld_between_conditional_temp'
out.dir <- 'ld_between_conditional_temp/summary'
if(!file.exists(out.dir)){
    dir.create(out.dir)
}

# Load data
res = do.call(rbind, lapply(list.files(data.dir, full.names=TRUE, pattern="*.txt.gz"), function(x){
    read.delim(x, sep = "", header=F)
    })) %>% 
    rename(SNPA=V3, SNPB=V6, R2=V7, condition=V8, gene=V9) %>% 
    filter(SNPA != SNPB) %>% 
    select(SNPA, SNPB, R2, condition, gene) %>%
    mutate(
        pair_id = paste(pmin(SNPA, SNPB), pmax(SNPA, SNPB), sep="_")
    ) %>% 
    distinct(
        pair_id, condition, gene, R2
    ) %>% arrange(R2)

# Plot the distribution of the whole set
p = ggplot(res, aes(x = R2)) + 
    geom_density(fill = "lightgrey", alpha = 0.4) +
    geom_vline(xintercept = median(res$R2, na.rm = TRUE), 
                linetype = "dashed", color = "red", size = 1) +
    geom_vline(xintercept = 0.5, 
                linetype = "dashed", color = "blue", size = 1) +
    labs(x = "LD between independent effects for same gene x condition", 
        y = "Density") +
    theme_classic()

ggsave(paste0(out.dir,"/LD_leads_distribution.png"), p, width = 4.8, height = 5)

# Where is the 95th percentile:
quantile(res$R2, 0.95)