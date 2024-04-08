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
library(gggenomes)

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

##### Pinopsin - PAML analysis  --------------------------------


#FreeRatio
crst <- treeio::read.paml_rst("FreeRatio/rst")
mlcfile <- read.codeml_mlc("FreeRatio/FreeRatio.mcl")

mRCA_polypteriformes <- 
  findMRCA(mlcfile@phylo, 
           c("Erpetoichthys_calabaricus---pinopsin-LR536439.2-78697483-78707040---5_exons",
             "Polypterus_senegalus---pinopsin-CM029060.1-117343073-117355106---5_exons")
  )


mlcfile@phylo <- 
  root(mlcfile@phylo, node=mRCA_polypteriformes)
mlcfile@phylo$tip.label <- gsub(pattern="---.*", "", mlcfile@phylo$tip.label) 



ggtree(mlcfile, branch.length="none", aes(color=dN_vs_dS), size=2) +
  geom_tiplab(size = 2, fontface=3) + 
  xlim(0, 13) +
  scale_color_continuous(name='dN/dS', limits=c(0, 1.5),
                         oob=scales::squish, low="darkgreen", high="red") +
  theme_tree2(legend.position=c(0.07, 0.8))



phylo_crst <- crst@phylo
mlcfile_df <- as.data.frame(mlcfile@data)
mlcfile_df <- mlcfile_df 
misc_table <- as.data.frame("Pluvialis_apricaria") %>% mutate(value = 1)
colnames(misc_table) <- c("node","value")
misc_table$node <- as.numeric(misc_table$node)
node_label_corresp <- as.data.frame(left_join(phylo_crst, misc_table, by = 'node') %>% dplyr::select(node, label))

mlcfile_df <- 
  left_join(mlcfile_df, node_label_corresp, by="node")


#Two_Ratio

crst <- treeio::read.paml_rst("TwoRatio/rst")
mlcfile <- read.codeml_mlc("TwoRatio/TwoRatio.mcl")


mRCA_polypteriformes <- 
  findMRCA(mlcfile@phylo, 
           c("Erpetoichthys_calabaricus---pinopsin-LR536439.2-78697483-78707040---5_exons",
             "Polypterus_senegalus---pinopsin-CM029060.1-117343073-117355106---5_exons")
  )
mlcfile@phylo <- 
  root(mlcfile@phylo, node=mRCA_polypteriformes)
mlcfile@phylo$tip.label <- gsub(pattern="---.*", "", mlcfile@phylo$tip.label) 



ggtree(mlcfile, branch.length="none", aes(color=dN_vs_dS), size=2) +
  geom_tiplab(size = 2, fontface=3) + 
  xlim(0, 13) +
  scale_color_continuous(name='dN/dS', limits=c(0, 1.5),
                         oob=scales::squish, low="darkgreen", high="red") +
  theme_tree2(legend.position="none")


phylo_crst <- crst@phylo
mlcfile_df <- as.data.frame(mlcfile@data)
mlcfile_df <- mlcfile_df %>% filter(dN_vs_dS > 0) %>% filter(dN_vs_dS < 3)
misc_table <- as.data.frame("Pluvialis_apricaria") %>% mutate(value = 1)
colnames(misc_table) <- c("node","value")
misc_table$node <- as.numeric(misc_table$node)
node_label_corresp <- as.data.frame(left_join(phylo_crst, misc_table, by = 'node') %>% dplyr::select(node, label))

mlcfile_df <- 
  left_join(mlcfile_df, node_label_corresp, by="node")


#One Ratio

crst <- treeio::read.paml_rst("OneRatio/rst")
mlcfile <- read.codeml_mlc("OneRatio/OneRatio.mcl")

mRCA_polypteriformes <- 
  findMRCA(mlcfile@phylo, 
           c("Erpetoichthys_calabaricus---pinopsin-LR536439.2-78697483-78707040---5_exons",
             "Polypterus_senegalus---pinopsin-CM029060.1-117343073-117355106---5_exons")
  )
mlcfile@phylo <- 
  root(mlcfile@phylo, node=mRCA_polypteriformes)

mlcfile@phylo$tip.label <- gsub(pattern="---.*", "", mlcfile@phylo$tip.label) 


ggtree(mlcfile, branch.length="none", aes(color=dN_vs_dS), size=2) +
  geom_tiplab(size = 2, fontface=3) + 
  xlim(0, 13) +
  scale_color_continuous(name='dN/dS', limits=c(0, 1.5),
                         oob=scales::squish, low="darkgreen", high="red") +
  theme_tree2(legend.position="none")


### Compare likelihood of the different models

oneratio_np <- as.numeric(
  system('grep "lnL" OneRatio/OneRatio.mcl | sed "s/.*np: //g" | sed "s/).*//g"', intern = TRUE))
freeratio_np <- as.numeric(
  system('grep "lnL" FreeRatio/FreeRatio.mcl | sed "s/.*np: //g" | sed "s/).*//g"', intern = TRUE))
tworatio_np <- as.numeric(
  system('grep "lnL" TwoRatio/TwoRatio.mcl | sed "s/.*np: //g" | sed "s/).*//g"', intern = TRUE))

oneratio_lnl <- as.numeric(
  system('grep "lnL" OneRatio/OneRatio.mcl | sed "s/.*): *//g" | sed "s/ .*//g"', intern = TRUE))
freeratio_lnl <- as.numeric(
  system('grep "lnL" FreeRatio/FreeRatio.mcl | sed "s/.*): *//g" | sed "s/ .*//g"', intern = TRUE))
tworatio_lnl <- as.numeric(
  system('grep "lnL" TwoRatio/TwoRatio.mcl | sed "s/.*): *//g" | sed "s/ .*//g"', intern = TRUE))

#FreeRatio vs OneRatio

df_free_vs_one <- freeratio_np - oneratio_np
LRT_free_vs_one <- 2 * ((freeratio_lnl) - (oneratio_lnl))
pchisq(LRT_free_vs_one, df_free_vs_one, lower.tail=FALSE)

#TwoRatio vs OneRatio

df_two_vs_one <- tworatio_np - oneratio_np
LRT_two_vs_one <- 2 * ((tworatio_lnl) - (oneratio_lnl))
pchisq(LRT_two_vs_one, df_two_vs_one, lower.tail=FALSE)

#FreeRatio vs TwoRatio

df_free_vs_two <- freeratio_np - tworatio_np
LRT_free_vs_two <- 2 * ((freeratio_lnl) - (tworatio_lnl))
pchisq(LRT_free_vs_two, df_free_vs_two, lower.tail=FALSE)



##### Pinopsin - Synteny analysis  --------------------------------


my_seqs <- read.table("cluster_ID_ordered.csv", sep=",")
colnames(my_seqs) <- c("bin_id", "seq_id", "length")
my_seqs <- tibble(my_seqs)

my_genes <- read.table("seq_cluster_info.csv", sep=",")
colnames(my_genes) <-  c("seq_id", "start", "end", "gene_name", "strand")
my_genes <- tibble(my_genes)

my_links <- read.table("link_table.csv", sep=",")
colnames(my_links) <-  c("seq_id", "start", "end", "seq_id2", "start2","end2", "OGG")
my_links <- tibble(my_links)

p <-  
  gggenomes(seqs=my_seqs, genes=my_genes, links=my_links) %>%
  flip_seqs("NC_053159.1-116343073-118355106") 



p + geom_seq() + geom_gene() + geom_link()


