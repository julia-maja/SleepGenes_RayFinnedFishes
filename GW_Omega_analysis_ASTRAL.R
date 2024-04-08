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
library(gprofiler2)
library(topGO)
library(gt)
library(orthogene)

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


slope_sign_color <- 
  c("positive"="#648FFF",
    "negative"="#FFB000")

#### Data load - Bad OGG to remove ------

opsins_OGG <- scan("list_opsin_OGG.txt", what="character") 

#### Data load - OGG count per species ------



OGG_table <- 
  read.table("Orthogroups.GeneCount.tsv",
             header=TRUE,
             sep="\t")


#Remove noVar and opsins OGG from the OGG table

OGG_table_table_filtered <- 
  OGG_table %>%
  filter(! Orthogroup %in% c(opsins_OGG))

#Total number of genes in each OGG


OGG_table_table_count <- 
  OGG_table_table_filtered %>%
  dplyr::select(Orthogroup, Total)

colnames(OGG_table_table_count) <- c("OGG", "Total")

OGG_table_table_count <- 
  OGG_table_table_count %>%
  mutate(half_total = floor(Total/2))


#### Data load - OGG omega pGLS ------

opsins_omega_OGG_pgls_df <- 
  read.table(
    "pGLS_opsins_omega_OGG_teleost.tsv",
    header=FALSE,
    sep="\t")

colnames(opsins_omega_OGG_pgls_df) <-
  c("Response", "R2", "pvalue", "lambda", "slope", "predictor")


opsins_omega_OGG_pgls_df <- 
  opsins_omega_OGG_pgls_df %>%
  filter(! predictor %in% c(opsins_OGG))

opsins_omega_OGG_pgls_df %>%
  group_by(Response) %>%
  summarise(n())




#### Data load - Danio gene names ------


OGG_DanioName_df <- 
  read.table("OGG_DanioName_table.csv",
             header=FALSE,
             sep=",")
colnames(OGG_DanioName_df) <- c("OGG", "gene_name")

#remove opsins OGG
OGG_DanioName_df <- 
  OGG_DanioName_df %>%
  filter(! OGG %in% c(opsins_OGG))

#convert all gene name to GOterm

#keep only OGG present in the dataset 
datasetOGG <- opsins_omega_OGG_pgls_df %>% pull(predictor) %>% unique()
OGG_DanioName_df <- 
  OGG_DanioName_df %>%
  filter(OGG %in% datasetOGG)




OGG_DanioName_GO_df <- 
  read.table(
    "OGG_DanioName_GO_df.csv",
    sep=",",
    header=TRUE,
    quote=""
  )



Allgenes.Hs <- 
  scan("Allgenes.Hs.txt",
       what="character")

#### Significance and R2 of OGGs -- ALL  ------

### Correct the pvalue with Bonferroni method

#First lets check non visual opsins 


unique_response <- opsins_omega_OGG_pgls_df %>% pull(Response) %>% unique()

opsins_omega_OGG_pgls_df_corr <- as.data.frame(NULL)
for(curr_response in unique_response){
  
  curr_df <- 
    opsins_omega_OGG_pgls_df %>%
    filter(Response == curr_response)
  
  
  initial_pvalues <- 
    as.numeric(
      curr_df %>% pull(pvalue)
    )
  
  #corrected_pvalues_FDR <- p.adjust(initial_pvalues, method = "fdr")
  corrected_pvalues_bonferroni <- p.adjust(initial_pvalues, method = "bonferroni")
  
  #curr_df <- 
  #  curr_df %>%
  #  mutate(FDR_pvalue = corrected_pvalues_FDR,
  #         Bonferroni_pvalue = corrected_pvalues_bonferroni)
  
  
  curr_df <- 
    curr_df %>%
    mutate(Bonferroni_pvalue = corrected_pvalues_bonferroni)
  
  opsins_omega_OGG_pgls_df_corr <- 
    rbind(opsins_omega_OGG_pgls_df_corr, curr_df)
  
}


#Extract significant OGG with non visual opsins  -- Positive 
positive_significant_OGG_nonvisual_df <- 
  opsins_omega_OGG_pgls_df_corr %>%
  filter(Response == "total_non_visual") %>%
  filter(slope > 0) %>%
  filter(Bonferroni_pvalue < 0.05) 
positive_significant_OGG_nonvisual_list <- 
  positive_significant_OGG_nonvisual_df$predictor

#Extract significant OGG with non visual opsins  -- Negative
negative_significant_OGG_nonvisual_df <- 
  opsins_omega_OGG_pgls_df_corr %>%
  filter(Response == "total_non_visual") %>%
  filter(slope < 0) %>%
  filter(Bonferroni_pvalue < 0.05) 
negative_significant_OGG_nonvisual_list <- 
  negative_significant_OGG_nonvisual_df$predictor


#Extract significant OGG with visual opsins  -- Positive 
positive_significant_OGG_visual_df <- 
  opsins_omega_OGG_pgls_df_corr %>%
  filter(Response == "total_visual") %>%
  filter(slope > 0) %>%
  filter(Bonferroni_pvalue < 0.05) 
positive_significant_OGG_visual_list <- 
  positive_significant_OGG_visual_df$predictor

#Extract significant OGG with visual opsins  -- Negative
negative_significant_OGG_visual_df <- 
  opsins_omega_OGG_pgls_df_corr %>%
  filter(Response == "total_visual") %>%
  filter(slope < 0) %>%
  filter(Bonferroni_pvalue < 0.05) 
negative_significant_OGG_visual_list <- 
  negative_significant_OGG_visual_df$predictor




#Extract significant OGG with ALL opsins  -- Positive 
positive_significant_OGG_ALL_df <- 
  opsins_omega_OGG_pgls_df_corr %>%
  filter(Response == "total_opsins") %>%
  filter(slope > 0) %>%
  filter(Bonferroni_pvalue < 0.05) 
positive_significant_OGG_ALL_list <- 
  positive_significant_OGG_ALL_df$predictor

#Extract significant OGG with ALL opsins  -- Negative
negative_significant_OGG_ALL_df <- 
  opsins_omega_OGG_pgls_df_corr %>%
  filter(Response == "total_opsins") %>%
  filter(slope < 0) %>%
  filter(Bonferroni_pvalue < 0.05) 
negative_significant_OGG_ALL_list <- 
  negative_significant_OGG_ALL_df$predictor





#### Plot the pGLS results - Non Visual -------

non_visual_pGLS_results <- 
  opsins_omega_OGG_pgls_df_corr %>%
  filter(Response == "total_non_visual")

non_visual_pGLS_results <- 
  non_visual_pGLS_results %>%
  mutate(slope_sign = if_else(
    slope > 0,
    "positive", 
    "negative"
  ))


non_visual_pGLS_results <- 
  non_visual_pGLS_results %>%
  mutate(R = if_else(
    slope > 0,
    sqrt(R2),
    -sqrt(R2)
  ))

non_visual_pGLS_results <- mutate_all(non_visual_pGLS_results, ~replace_na(.,0))

bonferoni_signif_treshhold <- 
  (0.05/nrow(non_visual_pGLS_results))

#nb_signif_pos <- nrow(non_visual_pGLS_results %>% filter(slope > 0) %>% filter(Bonferroni_pvalue < 0.05))
#top25_R2_tresh <- round(25*nb_signif_pos/100)
#top50_R2_tresh <- round(50*nb_signif_pos/100)


non_visual_pGLS_results %>%
  ggplot(., aes(x=-log(pvalue), y=R2, color=slope_sign)) +
  geom_point() +
  geom_vline(xintercept = -log(bonferoni_signif_treshhold), color="red", linetype="dashed") +
  #geom_hline(yintercept = top25_R2_tresh, color="#648FFF", linetype="dashed") +
  #geom_hline(yintercept = top50_R2_tresh, color="#648FFF", linetype="dashed") +
  scale_color_manual(values = slope_sign_color) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none") +
  ylim(0, 1)




### Extract name of top 10 positive genes


top_10_genes_pos <- 
  tail(non_visual_pGLS_results %>%
         filter(pvalue < bonferoni_signif_treshhold) %>%
         filter(slope > 0) %>%
         arrange(R2), 10) %>%
  pull(predictor)


corr_gene_names <- 
  OGG_DanioName_GO_df %>%
  filter(OGG %in% top_10_genes_pos) %>%
  dplyr::select(OGG, gene_name) %>%
  distinct() 
colnames(corr_gene_names) <- c("predictor", "gene_name")

left_join(corr_gene_names, non_visual_pGLS_results, by="predictor") %>%
  dplyr::select(predictor, gene_name, R2) %>%
  arrange(R2)


### Extract name of top 10 negative genes

top_10_genes_neg <- 
  tail(non_visual_pGLS_results %>%
         filter(pvalue < bonferoni_signif_treshhold) %>%
         filter(slope < 0) %>%
         arrange(R2), 10) %>%
  pull(predictor)


corr_gene_names <- 
  OGG_DanioName_GO_df %>%
  filter(OGG %in% top_10_genes_neg) %>%
  dplyr::select(OGG, gene_name) %>%
  distinct() 
colnames(corr_gene_names) <- c("predictor", "gene_name")

left_join(corr_gene_names, non_visual_pGLS_results, by="predictor") %>%
  dplyr::select(predictor, gene_name, R2, pvalue) %>%
  arrange(R2) %>%
  mutate(logpvalue = -log(pvalue))




#### Plot the pGLS results -  Visual -------

visual_pGLS_results <- 
  opsins_omega_OGG_pgls_df_corr %>%
  filter(Response == "total_visual")

visual_pGLS_results <- 
  visual_pGLS_results %>%
  mutate(slope_sign = if_else(
    slope > 0,
    "positive", 
    "negative"
  ))

bonferoni_signif_treshhold <- 
  (0.05/nrow(visual_pGLS_results))


visual_pGLS_results %>%
  ggplot(., aes(x=-log(pvalue), y=R2, color=slope_sign)) +
  geom_point() +
  geom_vline(xintercept = -log(bonferoni_signif_treshhold), color="red", linetype="dashed") +
  scale_color_manual(values = slope_sign_color) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none") +
  xlim(0, 25) +
  ylim(0, 1)




### Extract name of top 10 positive genes

top_10_genes_pos <- 
  tail(visual_pGLS_results %>%
         filter(pvalue < bonferoni_signif_treshhold) %>%
         filter(slope > 0) %>%
         arrange(R2), 30) %>%
  pull(predictor)


corr_gene_names <- 
  OGG_DanioName_GO_df %>%
  filter(OGG %in% top_10_genes_pos) %>%
  dplyr::select(OGG, gene_name) %>%
  distinct() 
colnames(corr_gene_names) <- c("predictor", "gene_name")

left_join(corr_gene_names, visual_pGLS_results, by="predictor") %>%
  dplyr::select(predictor, gene_name, R2, pvalue) %>%
  arrange(R2)


### Extract name of top 10 negative genes

top_10_genes_neg <- 
  tail(visual_pGLS_results %>%
         filter(pvalue < bonferoni_signif_treshhold) %>%
         filter(slope < 0) %>%
         arrange(R2), 30) %>%
  pull(predictor)


corr_gene_names <- 
  OGG_DanioName_GO_df %>%
  filter(OGG %in% top_10_genes_neg) %>%
  dplyr::select(OGG, gene_name) %>%
  distinct() 
colnames(corr_gene_names) <- c("predictor", "gene_name")

left_join(corr_gene_names, visual_pGLS_results, by="predictor") %>%
  dplyr::select(predictor, gene_name, R2, pvalue) %>%
  arrange(R2)


#### Plot the pGLS results -  ALL -------

all_pGLS_results <- 
  opsins_omega_OGG_pgls_df_corr %>%
  filter(Response == "total_opsins")

all_pGLS_results <- 
  all_pGLS_results %>%
  mutate(slope_sign = if_else(
    slope > 0,
    "positive", 
    "negative"
  ))

bonferoni_signif_treshhold <- 
  (0.05/nrow(all_pGLS_results))


all_pGLS_results %>%
  ggplot(., aes(x=-log(pvalue), y=R2, color=slope_sign)) +
  geom_point() +
  geom_vline(xintercept = -log(bonferoni_signif_treshhold), color="red", linetype="dashed") +
  scale_color_manual(values = slope_sign_color) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none") +
  xlim(0, 25) +
  ylim(0, 0.8)




### Extract name of top 10 positive genes


top_10_genes_pos <- 
  tail(all_pGLS_results %>%
         filter(pvalue < bonferoni_signif_treshhold) %>%
         filter(slope > 0) %>%
         arrange(R2), 10) %>%
  pull(predictor)


corr_gene_names <- 
  OGG_DanioName_GO_df %>%
  filter(OGG %in% top_10_genes_pos) %>%
  dplyr::select(OGG, gene_name) %>%
  distinct() 
colnames(corr_gene_names) <- c("predictor", "gene_name")

left_join(corr_gene_names, all_pGLS_results, by="predictor") %>%
  dplyr::select(predictor, gene_name, R2) %>%
  arrange(R2)


### Extract name of top 10 negative genes

top_10_genes_neg <- 
  tail(all_pGLS_results %>%
         filter(pvalue < bonferoni_signif_treshhold) %>%
         filter(slope < 0) %>%
         arrange(R2), 10) %>%
  pull(predictor)


corr_gene_names <- 
  OGG_DanioName_GO_df %>%
  filter(OGG %in% top_10_genes_neg) %>%
  dplyr::select(OGG, gene_name) %>%
  distinct() 
colnames(corr_gene_names) <- c("predictor", "gene_name")

left_join(corr_gene_names, all_pGLS_results, by="predictor") %>%
  dplyr::select(predictor, gene_name, R2) %>%
  arrange(R2)


#### GOST - Non visual -- Danio -------


#First extract gene name of significant OGG
nonvisual.positive.genename <- 
  OGG_DanioName_df %>%
  filter(OGG %in% positive_significant_OGG_nonvisual_list) %>%
  pull(gene_name)



nonvisual.negative.genename <- 
  OGG_DanioName_df %>%
  filter(OGG %in% negative_significant_OGG_nonvisual_list) %>%
  pull(gene_name)



#Now launch an enrichment analysis using gost


gost.nonvisual.positive <-  
  gost(query = nonvisual.positive.genename,
       organism = "drerio",
       custom_bg = OGG_DanioName_df$gene_name,
       significant=FALSE)


gost.nonvisual.negative <-  
  gost(query = nonvisual.negative.genename,
       organism = "drerio",
       custom_bg = OGG_DanioName_df$gene_name,
       significant=FALSE)


#Extract significant results

Significant_results_positive <- 
  gost.nonvisual.positive$result %>%
  arrange(p_value) %>%
  filter(p_value < 0.05) 


Significant_results_negative <- 
  gost.nonvisual.negative$result %>%
  arrange(p_value) %>%
  filter(p_value < 0.05) 



Significant_results_negative <- 
  Significant_results_negative %>%
  dplyr::select(term_id, source, term_name, p_value, intersection_size)

colnames(Significant_results_negative) <-
  c("term", "source", "name", "pvalue", "gene number")

Significant_results_negative$term <-
  gsub("REAC:", "", Significant_results_negative$term)
Significant_results_negative$pvalue <- 
  round(Significant_results_negative$pvalue, digits = 4)




Significant_results_negative <- 
  read.table("Significant_results_negative.csv",
             sep=",",
             header=TRUE)

Significant_results_negative <- 
  Significant_results_negative %>%
  filter(! name == "KEGG root term")
colnames(Significant_results_negative) <- c("term", "source", "name", "pvalue", "gene number")
gt_results_neg <- gt(Significant_results_negative)






#### GOST - Non visual -- Human -------

library(orthogene)

# Transform Danio gene name to Human gene name
nonvisual.positive.genename <- 
  OGG_DanioName_df %>%
  filter(OGG %in% positive_significant_OGG_nonvisual_list) %>%
  pull(gene_name)

nonvisual.negative.genename <- 
  OGG_DanioName_df %>%
  filter(OGG %in% negative_significant_OGG_nonvisual_list) %>%
  pull(gene_name)


mapped_genes_positive <- orthogene::map_genes(genes = nonvisual.positive.genename,
                                              species = "zebrafish", 
                                              drop_na = FALSE)
mapped_genes_negative <- orthogene::map_genes(genes = nonvisual.negative.genename,
                                              species = "zebrafish", 
                                              drop_na = FALSE)


mapped_genes_positive_df <- 
  as.data.frame(
    mapped_genes_positive %>%
      filter(! is.na(name)) %>%
      pull(name) %>%
      unique())

mapped_genes_negative_df <- 
  as.data.frame(
    mapped_genes_negative %>%
      filter(! is.na(name)) %>%
      pull(name) %>%
      unique())


colnames(mapped_genes_positive_df) <- c("Dr_gene")
rownames(mapped_genes_positive_df) <- mapped_genes_positive_df$Dr_gene

colnames(mapped_genes_negative_df) <- c("Dr_gene")
rownames(mapped_genes_negative_df) <- mapped_genes_negative_df$Dr_gene




OrthoMap_Pos <- orthogene::convert_orthologs(gene_df = mapped_genes_positive_df,
                                             gene_input = "rownames", 
                                             gene_output = "columns", 
                                             input_species = "zebrafish",
                                             output_species = "human",
                                             non121_strategy = "keep_both_species") 
OrthoMap_Neg <- orthogene::convert_orthologs(gene_df = mapped_genes_negative_df,
                                             gene_input = "rownames", 
                                             gene_output = "columns", 
                                             input_species = "zebrafish",
                                             output_species = "human",
                                             non121_strategy = "keep_both_species") 


nonvisual.positive.Hs <- OrthoMap_Pos$ortholog_gene
nonvisual.positive.Hs <- unname(nonvisual.positive.Hs)


nonvisual.negative.Hs <- OrthoMap_Neg$ortholog_gene
nonvisual.negative.Hs <- unname(nonvisual.negative.Hs)







gost.nonvisual.positive.HS <-  
  gost(query = nonvisual.positive.Hs,
       organism = "hsapiens",
       custom_bg = Allgenes.Hs,
       significant=FALSE)

gost.nonvisual.negative.HS <-  
  gost(query = nonvisual.negative.Hs,
       organism = "hsapiens",
       custom_bg = Allgenes.Hs,
       significant=FALSE)




Significant_results_positive_HS <- 
  gost.nonvisual.positive.HS$result %>%
  arrange(p_value) %>%
  filter(p_value < 0.05) ### NOTHING


Significant_results_negative_HS <- 
  gost.nonvisual.negative.HS$result %>%
  arrange(p_value) %>%
  filter(p_value < 0.05) 



Significant_results_positive <- 
  Significant_results_positive %>%
  dplyr::select(term_id, source, term_name, p_value, intersection_size)

colnames(Significant_results_positive) <-
  c("term", "source", "name", "pvalue", "gene number")

Significant_results_positive$term <-
  gsub("REAC:", "", Significant_results_positive$term)
Significant_results_positive$pvalue <- 
  round(Significant_results_positive$pvalue, digits = 4)


gt_results_pos <- gt(Significant_results_positive)



#### GOST - Visual -- Danio -------


#First extract gene name of significant OGG
visual.positive.genename <- 
  OGG_DanioName_df %>%
  filter(OGG %in% positive_significant_OGG_visual_list) %>%
  pull(gene_name)



visual.negative.genename <- 
  OGG_DanioName_df %>%
  filter(OGG %in% negative_significant_OGG_visual_list) %>%
  pull(gene_name)



#Now launch an enrichment analysis using gost


gost.visual.positive <-  
  gost(query = visual.positive.genename,
       organism = "drerio",
       custom_bg = OGG_DanioName_df$gene_name,
       significant=FALSE)


gost.visual.negative <-  
  gost(query = visual.negative.genename,
       organism = "drerio",
       custom_bg = OGG_DanioName_df$gene_name,
       significant=FALSE)


#Extract significant results

Significant_results_positive <- 
  gost.visual.positive$result %>%
  arrange(p_value) %>%
  filter(p_value < 0.05) ##nothing


Significant_results_negative <- 
  gost.visual.negative$result %>%
  arrange(p_value) %>%
  filter(p_value < 0.05)  ##nothing



Significant_results_positive <- 
  Significant_results_positive %>%
  dplyr::select(term_id, source, term_name, p_value, intersection_size)

colnames(Significant_results_positive) <-
  c("term", "source", "name", "pvalue", "gene number")

Significant_results_positive$term <-
  gsub("REAC:", "", Significant_results_positive$term)
Significant_results_positive$pvalue <- 
  round(Significant_results_positive$pvalue, digits = 4)


gt_results_pos <- gt(Significant_results_positive)




#### GOST - Visual -- Human -------


# Transform Danio gene name to Human gene name
visual.positive.genename <- 
  OGG_DanioName_df %>%
  filter(OGG %in% positive_significant_OGG_visual_list) %>%
  pull(gene_name)

visual.negative.genename <- 
  OGG_DanioName_df %>%
  filter(OGG %in% negative_significant_OGG_visual_list) %>%
  pull(gene_name)


mapped_genes_positive <- orthogene::map_genes(genes = visual.positive.genename,
                                              species = "zebrafish", 
                                              drop_na = FALSE)
mapped_genes_negative <- orthogene::map_genes(genes = visual.negative.genename,
                                              species = "zebrafish", 
                                              drop_na = FALSE)


mapped_genes_positive_df <- 
  as.data.frame(
    mapped_genes_positive %>%
      filter(! is.na(name)) %>%
      pull(name) %>%
      unique())

mapped_genes_negative_df <- 
  as.data.frame(
    mapped_genes_negative %>%
      filter(! is.na(name)) %>%
      pull(name) %>%
      unique())


colnames(mapped_genes_positive_df) <- c("Dr_gene")
rownames(mapped_genes_positive_df) <- mapped_genes_positive_df$Dr_gene

colnames(mapped_genes_negative_df) <- c("Dr_gene")
rownames(mapped_genes_negative_df) <- mapped_genes_negative_df$Dr_gene




OrthoMap_Pos <- orthogene::convert_orthologs(gene_df = mapped_genes_positive_df,
                                             gene_input = "rownames", 
                                             gene_output = "columns", 
                                             input_species = "zebrafish",
                                             output_species = "human",
                                             non121_strategy = "keep_both_species") 
OrthoMap_Neg <- orthogene::convert_orthologs(gene_df = mapped_genes_negative_df,
                                             gene_input = "rownames", 
                                             gene_output = "columns", 
                                             input_species = "zebrafish",
                                             output_species = "human",
                                             non121_strategy = "keep_both_species") 


visual.positive.Hs <- OrthoMap_Pos$ortholog_gene
visual.positive.Hs <- unname(visual.positive.Hs)


visual.negative.Hs <- OrthoMap_Neg$ortholog_gene
visual.negative.Hs <- unname(visual.negative.Hs)



gost.visual.positive.HS <-  
  gost(query = visual.positive.Hs,
       organism = "hsapiens",
       custom_bg = Allgenes.Hs,
       significant=FALSE)

gost.visual.negative.HS <-  
  gost(query = visual.negative.Hs,
       organism = "hsapiens",
       custom_bg = Allgenes.Hs,
       significant=FALSE)




Significant_results_positive_HS <- 
  gost.visual.positive.HS$result %>%
  arrange(p_value) %>%
  filter(p_value < 0.05) ### NOTHING


Significant_results_negative_HS <- 
  gost.visual.negative.HS$result %>%
  arrange(p_value) %>%
  filter(p_value < 0.05) ### NOTHING



Significant_results_positive <- 
  Significant_results_positive %>%
  dplyr::select(term_id, source, term_name, p_value, intersection_size)

colnames(Significant_results_positive) <-
  c("term", "source", "name", "pvalue", "gene number")

Significant_results_positive$term <-
  gsub("REAC:", "", Significant_results_positive$term)
Significant_results_positive$pvalue <- 
  round(Significant_results_positive$pvalue, digits = 4)


gt_results_pos <- gt(Significant_results_positive)



#### GOST - ALL -- Danio -------


#First extract gene name of significant OGG
ALL.positive.genename <- 
  OGG_DanioName_df %>%
  filter(OGG %in% positive_significant_OGG_ALL_list) %>%
  pull(gene_name)


ALL.negative.genename <- 
  OGG_DanioName_df %>%
  filter(OGG %in% negative_significant_OGG_ALL_list) %>%
  pull(gene_name)



#Now launch an enrichment analysis using gost


gost.ALL.positive <-  
  gost(query = ALL.positive.genename,
       organism = "drerio",
       custom_bg = OGG_DanioName_df$gene_name,
       significant=FALSE)


gost.ALL.negative <-  
  gost(query = ALL.negative.genename,
       organism = "drerio",
       custom_bg = OGG_DanioName_df$gene_name,
       significant=FALSE)


#Extract significant results

Significant_results_positive <- 
  gost.ALL.positive$result %>%
  arrange(p_value) %>%
  filter(p_value < 0.05)


Significant_results_negative <- 
  gost.ALL.negative$result %>%
  arrange(p_value) %>%
  filter(p_value < 0.05) ### NOTHING



Significant_results_positive <- 
  Significant_results_positive %>%
  dplyr::select(term_id, source, term_name, p_value, intersection_size)

colnames(Significant_results_positive) <-
  c("term", "source", "name", "pvalue", "gene number")

Significant_results_positive$term <-
  gsub("REAC:", "", Significant_results_positive$term)
Significant_results_positive$pvalue <- 
  round(Significant_results_positive$pvalue, digits = 4)


gt_results_pos <- gt(Significant_results_positive)




#### GOST - ALL -- Human -------


# Transform Danio gene name to Human gene name
ALL.positive.genename <- 
  OGG_DanioName_df %>%
  filter(OGG %in% positive_significant_OGG_ALL_list) %>%
  pull(gene_name)

ALL.negative.genename <- 
  OGG_DanioName_df %>%
  filter(OGG %in% negative_significant_OGG_ALL_list) %>%
  pull(gene_name)


mapped_genes_positive <- orthogene::map_genes(genes = ALL.positive.genename,
                                              species = "zebrafish", 
                                              drop_na = FALSE)
mapped_genes_negative <- orthogene::map_genes(genes = ALL.negative.genename,
                                              species = "zebrafish", 
                                              drop_na = FALSE)


mapped_genes_positive_df <- 
  as.data.frame(
    mapped_genes_positive %>%
      filter(! is.na(name)) %>%
      pull(name) %>%
      unique())

mapped_genes_negative_df <- 
  as.data.frame(
    mapped_genes_negative %>%
      filter(! is.na(name)) %>%
      pull(name) %>%
      unique())


colnames(mapped_genes_positive_df) <- c("Dr_gene")
rownames(mapped_genes_positive_df) <- mapped_genes_positive_df$Dr_gene

colnames(mapped_genes_negative_df) <- c("Dr_gene")
rownames(mapped_genes_negative_df) <- mapped_genes_negative_df$Dr_gene




OrthoMap_Pos <- orthogene::convert_orthologs(gene_df = mapped_genes_positive_df,
                                             gene_input = "rownames", 
                                             gene_output = "columns", 
                                             input_species = "zebrafish",
                                             output_species = "human",
                                             non121_strategy = "keep_both_species") 
OrthoMap_Neg <- orthogene::convert_orthologs(gene_df = mapped_genes_negative_df,
                                             gene_input = "rownames", 
                                             gene_output = "columns", 
                                             input_species = "zebrafish",
                                             output_species = "human",
                                             non121_strategy = "keep_both_species") 


ALL.positive.Hs <- OrthoMap_Pos$ortholog_gene
ALL.positive.Hs <- unname(ALL.positive.Hs)


ALL.negative.Hs <- OrthoMap_Neg$ortholog_gene
ALL.negative.Hs <- unname(ALL.negative.Hs)



gost.ALL.positive.HS <-  
  gost(query = ALL.positive.Hs,
       organism = "hsapiens",
       custom_bg = Allgenes.Hs,
       significant=FALSE)

gost.ALL.negative.HS <-  
  gost(query = ALL.negative.Hs,
       organism = "hsapiens",
       custom_bg = Allgenes.Hs,
       significant=FALSE)




Significant_results_positive_HS <- 
  gost.ALL.positive.HS$result %>%
  arrange(p_value) %>%
  filter(p_value < 0.05) ### NOTHING


Significant_results_negative_HS <- 
  gost.ALL.negative$result %>%
  arrange(p_value) %>%
  filter(p_value < 0.05) ### NOTHING



Significant_results_positive <- 
  Significant_results_positive %>%
  dplyr::select(term_id, source, term_name, p_value, intersection_size)

colnames(Significant_results_positive) <-
  c("term", "source", "name", "pvalue", "gene number")

Significant_results_positive$term <-
  gsub("REAC:", "", Significant_results_positive$term)
Significant_results_positive$pvalue <- 
  round(Significant_results_positive$pvalue, digits = 4)


gt_results_pos <- gt(Significant_results_positive)

