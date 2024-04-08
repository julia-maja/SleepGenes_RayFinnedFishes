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
library(treeio)
library(tidytree)
library(TDbook)
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

##### All Opsins Phylogeny  ---------------------------------



opsin_tree <- read.tree("Complete_opsins.aln.trimal.treefile.rooted") #rooted on iTOL


#Define monophyletic clades
exorh <- scan("Clades/exorh", what="character")
opn3 <- scan("Clades/opn3", what="character")
opn4m2 <- scan("Clades/opn4m2", what="character")
opn6 <- scan("Clades/opn6", what="character")
opn7b <- scan("Clades/opn7b", what="character")
opn8b <- scan("Clades/opn8b", what="character")
opn9 <- scan("Clades/opn9", what="character")
parapinopsin <- scan("Clades/parapinopsin", what="character")
pinopsin <- scan("Clades/pinopsin", what="character")
rh1 <- scan("Clades/rh1", what="character")
rrh <- scan("Clades/rrh", what="character")
sws2 <- scan("Clades/sws2", what="character")
tmt2 <- scan("Clades/tmt2", what="character")
va <- scan("Clades/va", what="character")
lws <- scan("Clades/lws", what="character")
opn4m1_3 <- scan("Clades/opn4m1_3", what="character")
opn4x <- scan("Clades/opn4x", what="character")
opn5 <- scan("Clades/opn5", what="character")
opn7a <- scan("Clades/opn7a", what="character")
opn8a <- scan("Clades/opn8a", what="character")
opn8c <- scan("Clades/opn8c", what="character")
parietopsin <- scan("Clades/parietopsin", what="character")
rgr <- scan("Clades/rgr", what="character")
rh2 <- scan("Clades/rh2", what="character")
sws1 <- scan("Clades/sws1", what="character")
tmt1 <- scan("Clades/tmt1", what="character")
tmt3 <- scan("Clades/tmt3", what="character")



#Define the MRCA of monophyletic clades
MRCA_exorh <- findMRCA(opsin_tree, tips=exorh, type="node")
MRCA_opn3 <- findMRCA(opsin_tree, tips=opn3, type="node")
MRCA_opn4m2 <- findMRCA(opsin_tree, tips=opn4m2, type="node")
MRCA_opn6 <- findMRCA(opsin_tree, tips=opn6, type="node")
MRCA_opn7b <- findMRCA(opsin_tree, tips=opn7b, type="node")
MRCA_opn8b <- findMRCA(opsin_tree, tips=opn8b, type="node")
MRCA_opn9 <- findMRCA(opsin_tree, tips=opn9, type="node")
MRCA_pinopsin <- findMRCA(opsin_tree, tips=pinopsin, type="node")
MRCA_rh1 <- findMRCA(opsin_tree, tips=rh1, type="node")
MRCA_rrh <- findMRCA(opsin_tree, tips=rrh, type="node")
MRCA_sws2 <- findMRCA(opsin_tree, tips=sws2, type="node")
MRCA_tmt2 <- findMRCA(opsin_tree, tips=tmt2, type="node")
MRCA_va <- findMRCA(opsin_tree, tips=va, type="node")
MRCA_lws <- findMRCA(opsin_tree, tips=lws, type="node")
MRCA_opn4m1_3 <- findMRCA(opsin_tree, tips=opn4m1_3, type="node")
MRCA_opn4x <- findMRCA(opsin_tree, tips=opn4x, type="node")
MRCA_opn5 <- findMRCA(opsin_tree, tips=opn5, type="node")
MRCA_opn7a <- findMRCA(opsin_tree, tips=opn7a, type="node")
MRCA_opn8a <- findMRCA(opsin_tree, tips=opn8a, type="node")
MRCA_opn8c <- findMRCA(opsin_tree, tips=opn8c, type="node")
MRCA_parapinopsin <- findMRCA(opsin_tree, tips=parapinopsin, type="node")
MRCA_parietopsin <- findMRCA(opsin_tree, tips=parietopsin, type="node")
MRCA_rgr <- findMRCA(opsin_tree, tips=rgr, type="node")
MRCA_rh2 <- findMRCA(opsin_tree, tips=rh2, type="node")
MRCA_sws1 <- findMRCA(opsin_tree, tips=sws1, type="node")
MRCA_tmt1 <- findMRCA(opsin_tree, tips=tmt1, type="node")
MRCA_tmt3 <- findMRCA(opsin_tree, tips=tmt3, type="node")



#Create a node number and node label correspondance table

opsin_tree_nodelabel <- makeNodeLabel(opsin_tree, method="number", prefix="Node")

species_test <- "Xiphias_gladius---rh1_exorh-CM029152.1-19693244-19696318---5_exons"
misc_table <- as.data.frame(species_test) %>% mutate(value = 1)
colnames(misc_table) <- c("node","value")
misc_table$node <- as.numeric(misc_table$node)
node_label_corresp <- 
  left_join(opsin_tree_nodelabel, misc_table, by = 'node') %>%
  dplyr::select(node, label)

node_nb <- 
  node_label_corresp %>%
  filter(., grepl("Node", label)) %>%
  pull(node)
species_names <- 
  node_label_corresp %>%
  filter(., !grepl("Node", label)) %>%
  pull(label)
species_and_nodes <- c(species_names, node_nb)
node_label_corresp <- node_label_corresp %>% mutate(tree_names = species_and_nodes)
node_label_corresp <- as.data.frame(node_label_corresp)



list_clades <- 
  c("exorh","lws","opn3","opn4m1_3","opn4m2","opn4x","opn5","opn6","opn7a",
    "opn7b","opn8a","opn8b","opn8c","opn9","parapinopsin","parietopsin",
    "pinopsin","rgr","rh1","rh2","rrh","sws1","sws2","tmt1","tmt2","tmt3","va")


exorh_df <- as.data.frame(exorh) %>% mutate(Clade = "exorh")
lws_df <- as.data.frame(lws) %>% mutate(Clade = "lws")
opn3_df <- as.data.frame(opn3) %>% mutate(Clade = "opn3")
opn4m1_3_df <- as.data.frame(opn4m1_3) %>% mutate(Clade = "opn4m1_3")
opn4m2_df <- as.data.frame(opn4m2) %>% mutate(Clade = "opn4m2")
opn4x_df <- as.data.frame(opn4x) %>% mutate(Clade = "opn4x")
opn5_df <- as.data.frame(opn5) %>% mutate(Clade = "opn5")
opn6_df <- as.data.frame(opn6) %>% mutate(Clade = "opn6")
opn7a_df <- as.data.frame(opn7a) %>% mutate(Clade = "opn7a")
opn7b_df <- as.data.frame(opn7b) %>% mutate(Clade = "opn7b")
opn8a_df <- as.data.frame(opn8a) %>% mutate(Clade = "opn8a")
opn8b_df <- as.data.frame(opn8b) %>% mutate(Clade = "opn8b")
opn8c_df <- as.data.frame(opn8c) %>% mutate(Clade = "opn8c")
opn9_df <- as.data.frame(opn9) %>% mutate(Clade = "opn9")
parapinopsin_df <- as.data.frame(parapinopsin) %>% mutate(Clade = "parapinopsin")
parietopsin_df <- as.data.frame(parietopsin) %>% mutate(Clade = "parietopsin")
pinopsin_df <- as.data.frame(pinopsin) %>% mutate(Clade = "pinopsin")
rgr_df <- as.data.frame(rgr) %>% mutate(Clade = "rgr")
rh1_df <- as.data.frame(rh1) %>% mutate(Clade = "rh1")
rh2_df <- as.data.frame(rh2) %>% mutate(Clade = "rh2")
rrh_df <- as.data.frame(rrh) %>% mutate(Clade = "rrh")
sws1_df <- as.data.frame(sws1) %>% mutate(Clade = "sws1")
sws2_df <- as.data.frame(sws2) %>% mutate(Clade = "sws2")
tmt1_df <- as.data.frame(tmt1) %>% mutate(Clade = "tmt1")
tmt2_df <- as.data.frame(tmt2) %>% mutate(Clade = "tmt2")
tmt3_df <- as.data.frame(tmt3) %>% mutate(Clade = "tmt3")
va_df <- as.data.frame(va) %>% mutate(Clade = "va")



colnames(exorh_df) <- c("tips_labels", "clade")
colnames(lws_df) <- c("tips_labels", "clade")
colnames(opn3_df) <- c("tips_labels", "clade")
colnames(opn4m1_3_df) <- c("tips_labels", "clade")
colnames(opn4m2_df) <- c("tips_labels", "clade")
colnames(opn4x_df) <- c("tips_labels", "clade")
colnames(opn5_df) <- c("tips_labels", "clade")
colnames(opn6_df) <- c("tips_labels", "clade")
colnames(opn7a_df) <- c("tips_labels", "clade")
colnames(opn7b_df) <- c("tips_labels", "clade")
colnames(opn8a_df) <- c("tips_labels", "clade")
colnames(opn8b_df) <- c("tips_labels", "clade")
colnames(opn8c_df) <- c("tips_labels", "clade")
colnames(opn9_df) <- c("tips_labels", "clade")
colnames(parapinopsin_df) <- c("tips_labels", "clade")
colnames(parietopsin_df) <- c("tips_labels", "clade")
colnames(pinopsin_df) <- c("tips_labels", "clade")
colnames(rgr_df) <- c("tips_labels", "clade")
colnames(rh1_df) <- c("tips_labels", "clade")
colnames(rh2_df) <- c("tips_labels", "clade")
colnames(rrh_df) <- c("tips_labels", "clade")
colnames(sws1_df) <- c("tips_labels", "clade")
colnames(sws2_df) <- c("tips_labels", "clade")
colnames(tmt1_df) <- c("tips_labels", "clade")
colnames(tmt2_df) <- c("tips_labels", "clade")
colnames(tmt3_df) <- c("tips_labels", "clade")
colnames(va_df) <- c("tips_labels", "clade")

All_tips_clades <- 
  rbind(exorh_df,
        lws_df,
        opn3_df,
        opn4m1_3_df,
        opn4m2_df,
        opn4x_df,
        opn5_df,
        opn6_df,
        opn7a_df,
        opn7b_df,
        opn8a_df,
        opn8b_df,
        opn8c_df,
        opn9_df,
        parapinopsin_df,
        parietopsin_df,
        pinopsin_df,
        rgr_df,
        rh1_df,
        rh2_df,
        rrh_df,
        sws1_df,
        sws2_df,
        tmt1_df,
        tmt2_df,
        tmt3_df,
        va_df)




df_nodes_clades <- as.data.frame(NULL)
for(curr_clade in list_clades){
  
  curr_nodes <- 
    getMRCA(opsin_tree_nodelabel, 
            All_tips_clades %>% 
              filter(clade == curr_clade) %>%
              pull(tips_labels) %>%
              unique())
  
  curr_df <-
    as.data.frame(
      getDescendants(opsin_tree_nodelabel, node=curr_nodes)) %>%
    mutate(Clade_curr = curr_clade)
  
  colnames(curr_df) <- c("node", "Clade")
  
  df_nodes_clades <- rbind(df_nodes_clades, curr_df)
}

Clades_for_tree <- 
  left_join(node_label_corresp, df_nodes_clades, by="node") %>%
  dplyr::select(label, Clade)




#Define a color palette for opsin families



ggtree(opsin_tree_nodelabel, layout="circular", size=0.6) %<+% Clades_for_tree +
  aes(color=Clade) +
  scale_color_manual(values = opsins_clades_colors) +
  theme(legend.position = "none")



## test collapse clades
exorh_random <- sample(exorh,1)
lws_random <- sample(lws,1)
opn3_random <- sample(opn3,1)
opn4m1_3_random <- sample(opn4m1_3,1)
opn4m2_random <- sample(opn4m2,1)
opn4x_random <- sample(opn4x,1)
opn5_random <- sample(opn5,1)
opn6_random <- sample(opn6,1)
opn7a_random <- sample(opn7a,1)
opn7b_random <- sample(opn7b,1)
opn8a_random <- sample(opn8a,1)
opn8b_random <- sample(opn8b,1)
opn8c_random <- sample(opn8c,1)
opn9_random <- sample(opn9,1)
parapinopsin_random <- sample(parapinopsin,1)
parietopsin_random <- sample(parietopsin,1)
pinopsin_random <- sample(pinopsin,1)
rgr_random <- sample(rgr,1)
rh1_random <- sample(rh1,1)
rh2_random <- sample(rh2,1)
rrh_random <- sample(rrh,1)
sws1_random <- sample(sws1,1)
sws2_random <- sample(sws2,1)
tmt1_random <- sample(tmt1,1)
tmt2_random <- sample(tmt2,1)
tmt3_random <- sample(tmt3,1)
va_random <- sample(va,1)


random_tips <- c(exorh_random,
                 lws_random,
                 opn3_random,
                 opn4m1_3_random,
                 opn4m2_random,
                 opn4x_random,
                 opn5_random,
                 opn6_random,
                 opn7a_random,
                 opn7b_random,
                 opn8a_random,
                 opn8b_random,
                 opn8c_random,
                 opn9_random,
                 parapinopsin_random,
                 parietopsin_random,
                 pinopsin_random,
                 rgr_random,
                 rh1_random,
                 rh2_random,
                 rrh_random,
                 sws1_random,
                 sws2_random,
                 tmt1_random,
                 tmt2_random,
                 tmt3_random,
                 va_random)


random_tips_space <- 
  gsub("_", " ", random_tips, ignore.case = FALSE, perl = FALSE,
       fixed = FALSE, useBytes = FALSE)

opsin_tree <- read.tree("converted.newick")
opsin_tree$tip.label <- gsub("'", "", opsin_tree$tip.label)


clade_tree <- keep.tip(opsin_tree, random_tips_space)



clade_tree$tip.label[clade_tree$tip.label==gsub("_", " ", exorh_random)] <-  "exorh"
clade_tree$tip.label[clade_tree$tip.label==gsub("_", " ", lws_random)] <-  "lws"
clade_tree$tip.label[clade_tree$tip.label==gsub("_", " ", opn3_random)] <-  "opn3"
clade_tree$tip.label[clade_tree$tip.label==gsub("_", " ", opn4m1_3_random)] <-  "opn4m1/3"
clade_tree$tip.label[clade_tree$tip.label==gsub("_", " ", opn4m2_random)] <-  "opn4m2"
clade_tree$tip.label[clade_tree$tip.label==gsub("_", " ", opn4x_random)] <-  "opn4x"
clade_tree$tip.label[clade_tree$tip.label==gsub("_", " ", opn5_random)] <-  "opn5"
clade_tree$tip.label[clade_tree$tip.label==gsub("_", " ", opn6_random)] <-  "opn6"
clade_tree$tip.label[clade_tree$tip.label==gsub("_", " ", opn7a_random)] <-  "opn7a"
clade_tree$tip.label[clade_tree$tip.label==gsub("_", " ", opn7b_random)] <-  "opn7b"
clade_tree$tip.label[clade_tree$tip.label==gsub("_", " ", opn8a_random)] <-  "opn8a"
clade_tree$tip.label[clade_tree$tip.label==gsub("_", " ", opn8b_random)] <-  "opn8b"
clade_tree$tip.label[clade_tree$tip.label==gsub("_", " ", opn8c_random)] <-  "opn8c"
clade_tree$tip.label[clade_tree$tip.label==gsub("_", " ", opn9_random)] <-  "opn9"
clade_tree$tip.label[clade_tree$tip.label==gsub("_", " ", parapinopsin_random)] <-  "parapinopsin"
clade_tree$tip.label[clade_tree$tip.label==gsub("_", " ", parietopsin_random)] <-  "parietopsin"
clade_tree$tip.label[clade_tree$tip.label==gsub("_", " ", pinopsin_random)] <-  "pinopsin"
clade_tree$tip.label[clade_tree$tip.label==gsub("_", " ", rgr_random)] <-  "rgr"
clade_tree$tip.label[clade_tree$tip.label==gsub("_", " ", rh1_random)] <-  "rh1"
clade_tree$tip.label[clade_tree$tip.label==gsub("_", " ", rh2_random)] <-  "rh2"
clade_tree$tip.label[clade_tree$tip.label==gsub("_", " ", rrh_random)] <-  "rrh"
clade_tree$tip.label[clade_tree$tip.label==gsub("_", " ", sws1_random)] <-  "sws1"
clade_tree$tip.label[clade_tree$tip.label==gsub("_", " ", sws2_random)] <-  "sws2"
clade_tree$tip.label[clade_tree$tip.label==gsub("_", " ", tmt1_random)] <-  "tmt1"
clade_tree$tip.label[clade_tree$tip.label==gsub("_", " ", tmt2_random)] <-  "tmt2"
clade_tree$tip.label[clade_tree$tip.label==gsub("_", " ", tmt3_random)] <-  "tmt3"
clade_tree$tip.label[clade_tree$tip.label==gsub("_", " ", va_random)] <- "va"


ggtree(clade_tree, branch.length="none", size=2) +
  xlim(0, 30) +
  geom_label2(aes(label=label, subset = !is.na(as.numeric(label)) & 
                    as.numeric(label) > 0)) +
  geom_tiplab() +
  coord_flip()



ggtree(clade_tree, size=2) +
  xlim(0, 5) +
  geom_label2(aes(label=label, subset = !is.na(as.numeric(label)) & 
                    as.numeric(label) > 0)) +
  geom_tiplab()





ggtree(clade_tree, branch.length="none", size=2) +
  xlim(0, 30) +
  geom_label2(aes(label=label, subset = !is.na(as.numeric(label)) & 
                    as.numeric(label) > 0),  label.size = 0.1) +
  geom_tiplab() +
  coord_flip()




ggtree(clade_tree, branch.length="none", size=2) +
  xlim(0, 30) +
  coord_flip()

