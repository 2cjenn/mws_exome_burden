
library(data.table)
library(R.utils)
library(dplyr)
library(tidyr)
library(stringr)
library(yaml)
library(here)

# Load the project config file for filepaths etc
config = read_yaml(here::here("./config.yml"))

# Source shared functions
source(config$scripts$shared)

pheno <-data.table::fread("/shared/MWS_Regeneron_Data/MWS_regeneron/Jennifer/exwas/phenofiles/bc_famhist_proxy_pheno.txt") %>% 
  select(IID, bc_famhist_proxy)

covars <- data.table::fread("/shared/MWS_Regeneron_Data/MWS_regeneron/Jennifer/exwas/covarfile/covars.txt") %>%
  select(IID, starts_with("pc"))

print_model <- function(data, formula, phenotype){
  warned <- 0
  withCallingHandlers(
    model <- glm(formula, family=binomial(link='logit'),data=data),
    warning = function(w) {warned <<- 1}
  )
  
  summ <- summary(model)
  coeff <- c(summ$coefficients[phenotype,], warned)
  return(coeff)
}

wilcox_method <- function(gene_vars, gene_list, outfile, pheno, covars, unrelated=FALSE) {
  write(x=paste0(c("SYMBOL", "Estimate", "SE", "z_score", "P", "Warning", "Unrelated"), collapse="\t"), file = outfile)
  phenotype <- "bc_famhist_proxy"
  if(unrelated) {
    pheno_unrel <- unrelated_max_subset(phenotype_file = pheno,
                                        pheno_want = phenotype)
  }
  
  for (i in seq(1, length(gene_list), by=1)) {
    gene <- gene_list[i]
    gene_variants <- gene_vars %>% filter(SYMBOL==gene)
    
    raw_vars <- extract_variants_bfile(
      gene_variants %>% pull(ID),
      bfile = "/shared/MWS_Regeneron_Data/MWS_regeneron/pVCF/exome_PLINK/OXFORD-MWS_Freeze_One.norm",
      plink_bin = config$software$plink2) 
    
    sum_burden <- as.data.frame(raw_vars) %>%
      select(-FID, -PAT, -MAT, -SEX, -PHENOTYPE) %>%
      mutate(across(starts_with("chr"), ~2-as.numeric(.x))) %>%
      pivot_longer(cols=starts_with("chr"), names_to="ID", values_to="genotype") %>%
      mutate(ID = str_remove(ID, "chr"),
             ID = str_remove(ID, "_[^.]*$")) %>%
      left_join(gene_variants, by="ID") %>%
      group_by(IID, SYMBOL) %>%
      summarise(burden = as.numeric(sum(!is.na(genotype) & genotype != 0) == 1), .groups="keep") %>%
      pivot_wider(names_from="SYMBOL", values_from="burden")
    
    print(paste0(gene, ": ", phenotype))
    formula <- paste0("`", gene, "` ~ ", phenotype, " + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10")
    
    data <- sum_burden %>% 
      inner_join(pheno %>% select(IID, bc_famhist_proxy), by="IID") %>%
      inner_join(covars, by="IID")
    
    coeff <- print_model(data=data, formula=formula, phenotype=phenotype)
    print(coeff)
    
    write(paste0(c(gene, coeff, 0), collapse="\t"), file=outfile, append=TRUE)
    
    if(unrelated) {
      data_unrel <- sum_burden %>% 
        inner_join(pheno_unrel %>% select(IID, bc_famhist_proxy), by="IID") %>%
        inner_join(covars, by="IID")
      
      coeff <- print_model(data=data_unrel, formula=formula, phenotype=phenotype)
      print(coeff)
      
      write(paste0(c(gene, coeff, 1), collapse="\t"), file=outfile, append=TRUE)
    }
  }
}


#######
# PTV # --------------------------------------------------------------------------------------------------------------------------
# Takes ages and might crash so only run if ready for that!
#######

all_genes <- data.table::fread("/shared/MWS_Regeneron_Data/MWS_regeneron/Jennifer/gene_panel/regenie/BCAC_masks/BCAC_rare_missense_anno.txt",
                               col.names=c("ID", "SYMBOL", "MASK"))
gene_vars <- all_genes %>%
  filter(MASK == "PTV")

gene_list <- gene_vars %>% pull(SYMBOL) %>% unique()
outfile <- "/shared/MWS_Regeneron_Data/MWS_regeneron/Jennifer/meta/FaH_Proxy_PTV.txt"

wilcox_method(gene_vars=gene_vars, gene_list=gene_list, outfile=outfile,
              pheno=pheno, covars=covars, unrelated=FALSE)



#################
# Rare missense # --------------------------------------------------------------------------------------------------------------------------
#################

wilcox_method(
  gene_vars=data.table::fread("/shared/MWS_Regeneron_Data/MWS_regeneron/Jennifer/gene_panel/regenie/BCAC_masks/BCAC_rare_missense_anno.txt",
                              col.names=c("ID", "SYMBOL", "MASK")) %>%
    filter(MASK == "rare_missense"), 
  gene_list=c("CHEK2", "SAMHD1"), 
  outfile="/shared/MWS_Regeneron_Data/MWS_regeneron/Jennifer/meta/FaH_Proxy_Rare.txt",
  pheno=pheno, covars=covars, unrelated=TRUE
)

#################
# CADD missense # --------------------------------------------------------------------------------------------------------------------------
#################

wilcox_method(
  gene_vars=data.table::fread("/shared/MWS_Regeneron_Data/MWS_regeneron/Jennifer/gene_panel/regenie/BCAC_masks/BCAC_cadd_missense_anno.txt",
                              col.names=c("ID", "SYMBOL", "MASK")), 
  # Don't filter MASK because the mask we want is PTV & CADD missense
  gene_list=c("CHEK2", "BRCA2", "PALB2", "BRCA1", "ATM"), 
  outfile="/shared/MWS_Regeneron_Data/MWS_regeneron/Jennifer/meta/FaH_Proxy_CADD.txt",
  pheno=pheno, covars=covars, unrelated=TRUE
)

##################
# Helix missense # --------------------------------------------------------------------------------------------------------------------------
##################


wilcox_method(
  gene_vars=data.table::fread("/shared/MWS_Regeneron_Data/MWS_regeneron/Jennifer/gene_panel/regenie/BCAC_masks/BCAC_helix_missense_anno.txt",
                              col.names=c("ID", "SYMBOL", "MASK")), 
  # Don't filter MASK because the mask we want is PTV & Helix missense
  gene_list=c("BRCA2", "BRCA1", "CHEK2", "PALB2", "ATM", "MAP3K1", "LZTR1"), 
  outfile="/shared/MWS_Regeneron_Data/MWS_regeneron/Jennifer/meta/FaH_Proxy_Helix.txt",
  pheno=pheno, covars=covars, unrelated=TRUE
)