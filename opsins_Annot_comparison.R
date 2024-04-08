##### Libraries  ---------------------------------

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

##### NCBI - Annotation number comparison  ---------------------------------

#Import annotation comparison table
Annotation_comp_df <- 
  read.table("Comparison_annotations_BUSCO.csv",
             sep=",",
             header=FALSE)
colnames(Annotation_comp_df) <- c("species", "NCBI_complete", "NCBI_incomplete", "study_complete", "study_incomplete")
category_annot_color <- 
  c("NCBI_more"="#DC267F",
    "study_more"="#648FFF",
    "Equal"="#FFB000")

Annotation_comp_df <- 
  Annotation_comp_df %>%
  mutate(NCBI_total = NCBI_complete + NCBI_incomplete) %>%
  mutate(study_total = study_complete + study_incomplete)


Annotation_comp_df <- 
  Annotation_comp_df %>%
  mutate(Category_total = case_when(
    NCBI_total > study_total ~ "NCBI_more",
    NCBI_total < study_total ~ "study_more",
    NCBI_total == study_total ~ "Equal"
  )) %>%
  mutate(Category_complete = case_when(
    NCBI_complete > study_complete ~ "NCBI_more",
    NCBI_complete < study_complete ~ "study_more",
    NCBI_complete == study_complete ~ "Equal"
  )) %>%
  mutate(Category_incomplete = case_when(
    NCBI_incomplete > study_incomplete ~ "NCBI_more",
    NCBI_incomplete < study_incomplete ~ "study_more",
    NCBI_incomplete == study_incomplete ~ "Equal"
  )) 



#Make a correlation graph between number of opsins retrieved in NCBI vs this study


cor_test_result <- 
  cor.test(
    Annotation_comp_df %>% pull(NCBI_total),
    Annotation_comp_df %>% pull(study_total),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)


Annotation_comp_df %>%
  ggplot(., aes(x=study_total,y=NCBI_total, color=Category_total)) +
  geom_jitter(size=2) +
  theme_classic() +
  geom_smooth(method = "lm", color = "black") +
  geom_abline(slope=1, linetype="dashed", color = "gray") +
  xlab("Number of opsins - Present study") +
  ylab("Number of opsins - NCBI") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  scale_color_manual(values = category_annot_color) +
  theme(legend.position="none")




#Extract species which have more in NCBI than in my study and check why ?

Annotation_comp_df %>%
  group_by(Category_total) %>%
  summarise(count = n())



cor_test_result <- 
  cor.test(
    Annotation_comp_df %>% pull(NCBI_complete),
    Annotation_comp_df %>% pull(study_complete),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)


Annotation_comp_df %>%
  ggplot(., aes(x=study_complete,y=NCBI_complete, color=Category_complete)) +
  geom_point(size=2) +
  theme_classic() +
  geom_smooth(method = "lm", color = "black") +
  geom_abline(slope=1, linetype="dashed", color = "gray") +
  xlab("Number of opsins - Present study") +
  ylab("Number of opsins - NCBI") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  scale_color_manual(values = category_annot_color)




cor_test_result <- 
  cor.test(
    Annotation_comp_df %>% pull(NCBI_incomplete),
    Annotation_comp_df %>% pull(study_incomplete),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)


Annotation_comp_df %>%
  ggplot(., aes(x=study_incomplete,y=NCBI_incomplete, color=Category_incomplete)) +
  geom_point(size=2) +
  theme_classic() +
  geom_smooth(method = "lm", color = "black") +
  geom_abline(slope=1, linetype="dashed", color = "gray") +
  xlab("Number of opsins - Present study") +
  ylab("Number of opsins - NCBI") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  scale_color_manual(values = category_annot_color)


##### NCBI - Annotation sequences comparison -- Boxplot  ---------------------------------

seq_identity_df <- 
  read.table("ALL_opsins.identity.csv",
             sep=",",
             header=FALSE)

colnames(seq_identity_df) <- 
  c("species", "subfamily", "gene_state", "annotated_gene_name","mygene_name","DNA_identity")


seq_identity_df$subfamily <-
  factor(seq_identity_df$subfamily ,
         levels=c("rh1", "exorh", "rh2", "sws2", "sws1", "lws", "pinopsin", "va", "parapinopsin",
                  "parietopsin", "tmt1", "tmt2", "tmt3","opn3", "opn8b", "opn8a", 
                  "opn8c", "opn7b", "opn7a", "opn5", "opn9", "opn6", "rgr", "rrh", "opn4m1_3", 
                  "opn4m2", "opn4x"))


seq_identity_df %>%
  ggplot(., aes(x=subfamily, y=DNA_identity, fill=subfamily)) +
  geom_boxplot() +
  scale_fill_manual(values = opsins_clades_colors) +
  xlab("Subfamily") +
  ylab("DNA identity") +
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  ylim(0, 1)



nrow(seq_identity_df %>%
  filter(DNA_identity == 1.000))

nrow(seq_identity_df %>%
       filter(DNA_identity > 0.9))

nrow(seq_identity_df)


##### NCBI - Annotation sequences comparison -- Erosion  ---------------------------------


subfamily_df <- 
  as.data.frame(enframe(c(seq_identity_df %>% pull(subfamily) %>% unique())) %>% dplyr::select(value))
misc <- as.data.frame("Total")
colnames(misc) <- "value"
subfamily_df <- 
  rbind(subfamily_df, misc)
colnames(subfamily_df) <- "subfamily"
seqid_erosion_df <- c()
for (i in seq(0, 1, 0.01)){
  
  
  current_sequences <- 
    seq_identity_df %>% 
    dplyr::filter(DNA_identity >= i) %>% 
    pull(mygene_name)
  
  
  current_df <- 
    seq_identity_df %>% 
    dplyr::filter(mygene_name %in% current_sequences) %>% 
    group_by(subfamily) %>% 
    summarise(number = n()) %>%
    adorn_totals("row") %>%
    mutate(DNA_identity_value = i)
  
  
  merged_df <- left_join(subfamily_df, current_df, by="subfamily")
  merged_df$number[is.na(merged_df$number)] <- 0
  merged_df$DNA_identity_value[is.na(merged_df$DNA_identity_value)] <- i
  
  seqid_erosion_df <- rbind(seqid_erosion_df, merged_df)
  
  
}




seqid_erosion_df %>%
  filter(subfamily != "Total") %>%
  ggplot(., aes(x=DNA_identity_value, y=number, color=subfamily)) +
  geom_point() +
  theme_classic() +
  geom_line() +
  ylab("Number of sequences") +
  xlab("DNA identity")+
  scale_color_manual(values = opsins_clades_colors) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none")

seqid_erosion_df %>%
  filter(subfamily != "Total") %>%
  ggplot(., aes(x=DNA_identity_value, y=number, color=subfamily)) +
  theme_classic() +
  geom_line() +
  ylab("Number of sequences") +
  xlab("DNA identity")+
  scale_color_manual(values = opsins_clades_colors) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none")


seqid_erosion_df %>%
  filter(subfamily != "Total") %>%
  ggplot(., aes(x=DNA_identity_value, y=number, color=subfamily)) +
  geom_point() +
  theme_classic() +
  geom_line() +
  ylab("Number of sequences") +
  xlab("DNA identity")+
  scale_color_manual(values = opsins_clades_colors) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) 


#More than 80% (81.9%) of sequences have 100% sequence identity with the NCBI annotation. 
#More than 98% (98.5%) of sequence have >90% sequence identity with the NCBI annotation.
#More than 95% (95.8%) of sequence have >95% sequence identity with the NCBI annotation.


nrow(seq_identity_df %>%
       filter(DNA_identity == 1.000))

nrow(seq_identity_df %>%
       filter(DNA_identity > 0.95))

nrow(seq_identity_df %>%
       filter(DNA_identity > 0.9))

nrow(seq_identity_df)


##### Musilova and Coresti comparison  ---------------------------------

#Import annotation comparison table
Musilova_comp_df <- 
  read.table("Musilova_and_Coresti.tsv",
             sep="\t",
             header=TRUE)
colnames(Musilova_comp_df) <- c("species", "order", "musilova_genome", "m_sws1", "m_sws2", "m_lws", "m_rh2", "total")

#import my table and merge

visual_opsins_table_count_wide <- 
  read.table("visual_opsins_table_count_wide.tsv", sep="\t", header=TRUE)

species_table <- 
  read.table("species_table.tsv", sep="\t", header=TRUE)



visual_opsins_table_count_wide_plus_acc <- 
  left_join(visual_opsins_table_count_wide, species_table, by="species")


My_study_cones <- 
  visual_opsins_table_count_wide_plus_acc %>%
  dplyr::select(species, sws1, sws2, lws, rh2, assembly_accession) %>%
  mutate(total_my_study = sws1 + sws2 + lws + rh2)

My_study_cones <- 
  left_join(My_study_cones, Musilova_comp_df, by="species") %>%
  filter(! is.na(total)) 


My_study_cones$musilova_genome <- sub("GCF_", "GCA_", My_study_cones$musilova_genome) 

My_study_cones <- 
  My_study_cones %>%
  mutate(accession_diff = if_else(
    assembly_accession == musilova_genome,
    "Same_assembly",
    "Different_assembly"
  ))


category_annot_color <- 
  c("Musilova_more"="#DC267F",
    "study_more"="#648FFF",
    "Equal"="#FFB000")


My_study_cones <- 
  My_study_cones %>%
  mutate(Category_total = case_when(
    total > total_my_study ~ "Musilova_more",
    total < total_my_study ~ "study_more",
    total == total_my_study ~ "Equal"
  )) 


My_study_cones %>%
  filter(Category_total == "Musilova_more") 

My_study_cones %>%
  filter(Category_total == "Equal")

My_study_cones %>%
  filter(Category_total == "study_more")

My_study_cones %>%
  filter(Category_total == "Musilova_more") %>%
  dplyr::select(m_sws1, sws1)
#Make a correlation graph between number of opsins retrieved in NCBI vs this study


cor_test_result <- 
  cor.test(
    My_study_cones %>% pull(total),
    My_study_cones %>% pull(total_my_study),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)


My_study_cones %>%
  ggplot(., aes(x=total_my_study, y=total, color=Category_total)) +
  geom_jitter(size=2, aes(shape=accession_diff)) +
  theme_classic() +
  geom_smooth(method = "lm", color = "black") +
  geom_abline(slope=1, linetype="dashed", color = "gray") +
  xlab("Number of cone opsins - Present study") +
  ylab("Number of cone opsins - Musilova and Cortesi 2023") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  scale_color_manual(values = category_annot_color) +
  theme(legend.position="none")





cor_test_result_only_same <- 
  cor.test(
    My_study_cones %>% filter(accession_diff == "Same_assembly") %>% pull(total),
    My_study_cones %>% filter(accession_diff == "Same_assembly") %>% pull(total_my_study),
    method="pearson"
  )
estimate_p <- round(cor_test_result_only_same$estimate, 3)
pvalue_p <- cor_test_result_only_same$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)




##### Emile et al comparison  ---------------------------------

#Import annotation comparison table
Emile_comp_df <- 
  read.table("Emile_et_al_2017.tsv",
             sep="\t",
             header=FALSE)
#the table was modified according to the legend of the table in emile et al 2017. 
colnames(Emile_comp_df) <- c("species",
                             "exorh", "pinopsin", "va", "parietopsin", "parapinopsin",
                             "opn3", "tmt2", "tmt3", "tmt1", "opn4x", "opn4m", "opn5_9",
                             "opn6a", "opn6b", "opn7a", "opn7","opn8",
                             "rrh", "rgr", "total_emile")

non_visual_opsins_table_count_wide <- 
  read.table("non_visual_opsins_table_count_wide.tsv", sep="\t", header=TRUE)


non_visual_opsins_table_count_wide_plus_acc <- 
  left_join(non_visual_opsins_table_count_wide, species_table, by="species")



colnames(non_visual_opsins_table_count_wide_plus_acc)  <- 
  sub("-", "_", colnames(non_visual_opsins_table_count_wide_plus_acc))

My_study_non_visual <- 
  non_visual_opsins_table_count_wide_plus_acc %>%
  mutate(total_my_study = exorh + opn3 + opn4m1_3 + opn4m2 + opn4x + opn5 + opn6 + opn7a + 
           opn7b + opn8a + opn8b + opn8c + opn9 + parapinopsin + 
           parietopsin + rgr + rrh + tmt1 + tmt2 + tmt3 + va + pinopsin)

My_study_non_visual <- 
  left_join(My_study_non_visual, Emile_comp_df, by="species") %>%
  filter(! is.na(total_emile)) 



category_annot_color <- 
  c("Emile_more"="#DC267F",
    "study_more"="#648FFF",
    "Equal"="#FFB000")

My_study_non_visual <- 
  My_study_non_visual %>%
  mutate(Category_total = case_when(
    total_emile > total_my_study ~ "Emile_more",
    total_emile < total_my_study ~ "study_more",
    total_emile == total_my_study ~ "Equal"
  )) 





#Make a correlation graph between number of opsins retrieved in NCBI vs this study


cor_test_result <- 
  cor.test(
    My_study_non_visual %>% pull(total_emile),
    My_study_non_visual %>% pull(total_my_study),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)


My_study_non_visual %>%
  ggplot(., aes(x=total_my_study,y=total_emile, color=Category_total)) +
  geom_point(size=2) +
  theme_classic() +
  geom_smooth(method = "lm", color = "black") +
  geom_abline(slope=1, linetype="dashed", color = "gray") +
  xlab("Number of cone opsins - Present study") +
  ylab("Number of cone opsins - Beaudry et al. 2017") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  scale_color_manual(values = category_annot_color) +
  theme(legend.position="none")











