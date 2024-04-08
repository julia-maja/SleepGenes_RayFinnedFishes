##### Libraries  ---------------------------------
rm(list=ls())

set.seed(2712)


library("caper")
library("patchwork")
library("ape")
library(dplyr)
library(ggplot2)
library(tidyverse)
library("lattice")
library(reshape2)
library(adephylo)
library(phylobase)
library(data.table)
library(phytools)
library("corrplot")
library("geiger")
library("ggpubr")
library(phylolm)
library(janitor)
library(RColorBrewer)
library(ggtreeExtra)
library(rcartocolor)
library(ggtree)
library(ggcorrplot)

split_tibble <- function(tibble, column = 'col') {
  tibble %>% split(., .[,column]) %>% lapply(., function(x) x[,setdiff(names(x),column)])
}



#### My functions + Colors palettes  ---------------------------------

PGLS_pvalue <- function(pgls_rslt) {
  sum_cor <- summary(pgls_rslt)
  pvalue =  formatC(sum_cor$coefficients[8], digits = 3)
  if (pvalue == "   0"){ pvalue = 2.2e-16}
  return(pvalue)
}

PGLS_R2 <- function(pgls_rslt) {
  sum_cor <- summary(pgls_rslt)
  r2 = formatC(sum_cor$r.squared, digits = 2)
  return(r2)
}

PGLS_lambda <- function(pgls_rslt) {
  sum_cor <- summary(pgls_rslt)
  PGLS_lambda = sum_cor$param[2]
  return(PGLS_lambda)
}


GLS_function <- function(gls_rslt) {
  GLS_cc <- coef(gls_rslt)
  GLS_f <- function(x) GLS_cc[1] + GLS_cc[2]*x
  return(GLS_f)
}

library(RColorBrewer)
n <- 30
colrs <- brewer.pal.info[brewer.pal.info$colorblind == TRUE, ]
col_vec = unlist(mapply(brewer.pal, colrs$maxcolors, rownames(colrs)))
col <- sample(col_vec, n)
area <- rep(1,n)
pie(area, col = col)

opsins_clades_colors <- 
  c("exorh"="#E6F5D0",
    "lws"="#67001F",
    "opn3"="#F4A582",
    "opn4m1_3"="#35978F",
    "opn4m2"="#238443",
    "opn4x"="#525252",
    "opn5"="#9970AB",
    "opn6"="#B15928",
    "opn7a"="#08519C",
    "opn7b"="#EF6548",
    "opn8a"="#B3B3B3",
    "opn8b"="#2D004B",
    "opn8c"="goldenrod",
    "opn9"="#88419D",
    "parapinopsin"="#A8DDB5",
    "parietopsin"="#FDAE61",
    "pinopsin"="#3690C0",
    "rgr"="#C2A5CF",
    "rh1"="black",
    "rh2"="#276419",
    "rrh"="#ABD9E9",
    "sws1"="#FEE391",
    "sws2"="#41B6C4",
    "tmt1"="#FC4E2A",
    "tmt2"="#045A8D",
    "tmt3"="#D1E5F0",
    "va"="#8C510A"
  )

##### Data load - Species info table -----------------

species_table <- 
  read.table("Actino_species_table_info_v3.tsv",
             sep="\t",
             header=TRUE)
species_table <- 
  species_table %>% 
  filter(genome_name != "allo-octoploid_hybrid_of_Carassius_auratus_x_Cyprinus_carpio")
species_table[(species_table$genome_name == "Taenioides_sp._WSHXM2023"),"genome_name"] <- "Taenioides_sp_WSHXM2023"
species_table <- 
  species_table %>%
  mutate(perc_BUSCO = BUSCO_complete/BUSCO_searched)



## Add information about recent WGD 

species_table <- 
  species_table %>%
  mutate(WGD = if_else(
    family %in% c("Salmonidae", "Cyprinidae", "Acipenseridae", "Polyodontidae", "Catostomidae"),
    "Yes", 
    "No"
  ))

species_table[(species_table$genome_name == "Labeo_rohita"),"WGD"] <- "No"
species_table[(species_table$genome_name == "Labeo_catla"),"WGD"] <- "No"
species_table[(species_table$genome_name == "Acrossocheilus_wenchowensis"),"WGD"] <- "No"
species_table[(species_table$genome_name == "Paedocypris_micromegethes"),"WGD"] <- "No"
species_table[(species_table$genome_name == "Gymnocypris_eckloni_scoliostomus"),"WGD"] <- "No"
species_table[(species_table$genome_name == "Schizopygopsis_malacanthus"),"WGD"] <- "No"
species_table[(species_table$genome_name == "Schizopygopsis_pylzovi"),"WGD"] <- "No"
species_table[(species_table$genome_name == "Schizothorax_lantsangensis"),"WGD"] <- "No"
species_table[(species_table$genome_name == "Aspiorhynchus_laticeps"),"WGD"] <- "No"
species_table[(species_table$genome_name == "Puntigrus_tetrazona"),"WGD"] <- "No"
species_table[(species_table$genome_name == "Onychostoma_macrolepis"),"WGD"] <- "No"




##### Data load - Opsins table Original  ---------------------------------

opsins_table <- 
  read.table("ALL_opsins_classifications.csv",
             sep=",",
             header=FALSE)

colnames(opsins_table) <- c("gene_state", "gene_name", "clade", "species")


species_table_sel <- species_table %>% dplyr::select(genome_name, assembly_accession)
colnames(species_table_sel) <- c("species", "accession")

opsins_table <- 
  left_join(opsins_table, species_table_sel, by="species")

opsins_table_complete <- 
  opsins_table %>% filter(gene_state == "Complete")

## count the number of complete genes

summary_opsins_count <- 
  as.data.frame(
  opsins_table_complete %>%
  group_by(species, accession, clade) %>%
  summarise(count = n())
  )

opsins_table_complete_wide <- 
  as.data.frame(
    summary_opsins_count %>%
      pivot_wider(names_from = clade, values_from = count)
  )

opsins_table_complete_wide <- mutate_all(opsins_table_complete_wide, ~replace_na(.,0))



##### Data load - BUSCO alternative assemblies  ---------------------------------

BUSCO_alternate <- 
  read.table("BUSCO_results_alternative_assemblies.tsv",
             sep="\t",
             header=TRUE)

BUSCO_alternate <- 
  BUSCO_alternate %>%
  mutate(perc_BUSCO = BUSCO_complete/BUSCO_searched)


BUSCO_erosion_df <- c()
for (i in seq(0, 1, 0.01)){
  
  current_df <- 
    BUSCO_alternate %>% 
    dplyr::filter(perc_BUSCO >= i) %>% 
    summarise(number = n()) %>%
    mutate(Busco_percentage = i)
  

  BUSCO_erosion_df <- rbind(BUSCO_erosion_df, current_df)
  
  
}



#Define the number of species filtered 
BUSCO_90_nb <- BUSCO_erosion_df %>% filter(Busco_percentage == 0.9) %>% pull(number)


BUSCO_erosion_df %>%
  ggplot(., aes(x=Busco_percentage, y=number)) +
  geom_point() +
  theme_classic() +
  geom_line() +
  ylab("Number of genomes") +
  xlab("Proportion of complete BUSCO genes")+
  geom_segment(aes(x = 0.9, y = 0, xend = 0.9, yend = BUSCO_90_nb),color="#DC267F", linetype="dashed") +
  geom_segment(aes(x = 0, y = BUSCO_90_nb, xend = 0.9, yend = BUSCO_90_nb),color="#DC267F", linetype="dashed") +
  #geom_segment(aes(x = 0.8, y = 0, xend = 0.8, yend = BUSCO_80_nb_nonteleost),color="#FFB000", linetype="dashed") +
  #geom_segment(aes(x = 0, y = BUSCO_80_nb_nonteleost, xend = 0.8, yend = BUSCO_80_nb_nonteleost),color="#FFB000", linetype="dashed") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none")


##### Data load - Opsins alternative assemblies   ---------------------------------


opsins_table_alternate <- 
  read.table("Robustness_opsins_classifications.csv",
             sep=",",
             header=FALSE)

colnames(opsins_table_alternate) <- c("gene_state", "gene_name", "clade", "species", "accession")


## count the number of complete genes

summary_opsins_count_alternate <- 
  as.data.frame(
    opsins_table_alternate %>%
      group_by(species, accession, clade) %>%
      summarise(count = n())
  )

opsins_table_complete_wide_alternate <- 
  as.data.frame(
    summary_opsins_count_alternate %>%
      pivot_wider(names_from = clade, values_from = count)
  )

opsins_table_complete_wide_alternate <- mutate_all(opsins_table_complete_wide_alternate, ~replace_na(.,0))




##### Correlation between counts on the same species assemblies  ---------------------------------

opsins_table_complete_wide <- 
  opsins_table_complete_wide %>%
  mutate(VO = rh1 + rh2 + lws + sws1 + sws2) %>%
  mutate(NVO = exorh + opn3 + opn4m1_3 + opn4m2 + opn4x + opn5 + opn6 + opn7a + opn7b + opn8a +
           opn8b + opn8c + opn9 + parapinopsin + parietopsin + rgr + rrh + tmt1 + tmt2 + tmt3 + 
           va + pinopsin)


opsins_table_complete_wide_alternate <- 
  opsins_table_complete_wide_alternate %>%
  mutate(VO = rh1 + rh2 + lws + sws1 + sws2) %>%
  mutate(NVO = exorh + opn3 + opn4m1_3 + opn4m2 + opn4x + opn5 + opn6 + opn7a + opn7b + opn8a +
           opn8b + opn8c + opn9 + parapinopsin + parietopsin + rgr + rrh + tmt1 + tmt2 + tmt3 + 
           va + pinopsin)


list_uniq_species <- opsins_table_complete_wide_alternate %>% pull(species) %>% unique()



NVO_comparison_df <- as.data.frame(NULL)
VO_comparison_df <- as.data.frame(NULL)
for(curr_species in list_uniq_species){
  
  study_NVO_count <- 
    opsins_table_complete_wide %>%
    filter(species == curr_species) %>%
    pull(NVO)
     
  study_VO_count <- 
    opsins_table_complete_wide %>%
    filter(species == curr_species) %>%
    pull(VO)   
  
  
  alternative_counts_NVO <- 
    opsins_table_complete_wide_alternate %>%
    filter(species == curr_species) %>%
    pull(NVO)
  
  alternative_counts_VO <- 
    opsins_table_complete_wide_alternate %>%
    filter(species == curr_species) %>%
    pull(VO)
  
  
  
  for(NVO_c in alternative_counts_NVO){
    curr_NVO_df <- as.data.frame(cbind(curr_species, study_NVO_count, NVO_c))
    colnames(curr_NVO_df) <- c("species", "study_NVO", "alternative_NVO")
    NVO_comparison_df <- rbind(NVO_comparison_df, curr_NVO_df)
  }
  
  for(VO_c in alternative_counts_VO){
    curr_VO_df <- as.data.frame(cbind(curr_species, study_VO_count, VO_c))
    colnames(curr_VO_df) <- c("species", "study_VO", "alternative_VO")
    VO_comparison_df <- rbind(VO_comparison_df, curr_VO_df)
  }

}


VO_comparison_df$study_VO <- as.numeric(VO_comparison_df$study_VO)
VO_comparison_df$alternative_VO <- as.numeric(VO_comparison_df$alternative_VO)

NVO_comparison_df$study_NVO <- as.numeric(NVO_comparison_df$study_NVO)
NVO_comparison_df$alternative_NVO <- as.numeric(NVO_comparison_df$alternative_NVO)







VO_comparison_df <- 
  VO_comparison_df %>%
  mutate(Category_total = case_when(
    alternative_VO > study_VO ~ "alternative_more",
    alternative_VO < study_VO ~ "study_more",
    alternative_VO == study_VO ~ "Equal"
  )) 


NVO_comparison_df <- 
  NVO_comparison_df %>%
  mutate(Category_total = case_when(
    alternative_NVO > study_NVO ~ "alternative_more",
    alternative_NVO < study_NVO ~ "study_more",
    alternative_NVO == study_NVO ~ "Equal"
  )) 


category_annot_color <- 
  c("alternative_more"="#DC267F",
    "study_more"="#648FFF",
    "Equal"="#FFB000")


NVO_comparison_df %>% mutate(abs_diff = abs(study_NVO - alternative_NVO)) %>% arrange(abs_diff)


cor_test_result <- 
  cor.test(
    NVO_comparison_df %>% pull(study_NVO),
    NVO_comparison_df %>% pull(alternative_NVO),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)
lm_NVO <- 
  lm(data = NVO_comparison_df, 
     formula = alternative_NVO ~ study_NVO)
summary(lm_NVO)

NVO_comparison_df %>%
  ggplot(., aes(x=study_NVO, y=alternative_NVO, color=Category_total)) +
  geom_jitter(size=2) +
  theme_classic() +
  scale_color_manual(values = category_annot_color) + 
  geom_smooth(method = "lm", color = "black") +
  geom_abline(slope=1, linetype="dashed", color = "gray") +
  xlab("Number of non-visual opsins - Primary assembly") +
  ylab("Number of non-visual opsins - Alternative assembly") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none")

cor_test_result <- 
  cor.test(
    VO_comparison_df %>% pull(study_VO),
    VO_comparison_df %>% pull(alternative_VO),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)

lm_VO <- 
  lm(data = VO_comparison_df, 
   formula = alternative_VO ~ study_VO)
summary(lm_VO)


VO_comparison_df %>%
  ggplot(., aes(x=study_VO, y=alternative_VO, color=Category_total)) +
  geom_jitter(size=2) +
  scale_color_manual(values = category_annot_color) + 
  theme_classic() +
  geom_smooth(method = "lm", color = "black") +
  geom_abline(slope=1, linetype="dashed", color = "gray") +
  xlab("Number of visual opsins - Primary assembly") +
  ylab("Number of visual opsins - Alternative assembly") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none")




VO_comparison_df %>%
  group_by(Category_total) %>%
  summarise(n())
NVO_comparison_df %>%
  group_by(Category_total) %>%
  summarise(n())

##### Absolute difference VO and NVO  ---------------------------------

NVO_abs_diff_df <- 
  NVO_comparison_df %>%
  mutate(abs_diff = abs(study_NVO - alternative_NVO)) %>%
  dplyr::select(abs_diff) %>%
  mutate(NVO_VO = "NVO")

VO_abs_diff_df <- 
  VO_comparison_df %>%
  mutate(abs_diff = abs(study_VO - alternative_VO)) %>%
  dplyr::select(abs_diff) %>%
  mutate(NVO_VO = "VO")


All_diff_counts_df <- 
  rbind(VO_abs_diff_df, NVO_abs_diff_df)



All_diff_counts_df %>%
  ggplot(., aes(x=NVO_VO, y=abs_diff)) +
  geom_boxplot() +
  theme_classic() +
  xlab("") +
  ylab("Absolute difference") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none")



All_diff_counts_df %>%
  group_by(NVO_VO) %>%
  summarise(mean_diff = mean(abs_diff),
            min_diff = min(abs_diff),
            max_diff = max(abs_diff))





NVO_comparison_df %>%
  mutate(abs_diff = abs(study_NVO - alternative_NVO)) %>%
  filter(abs_diff == 7)
  
VO_comparison_df %>%
  mutate(abs_diff = abs(study_VO - alternative_VO)) %>%
  filter(abs_diff == 1)

