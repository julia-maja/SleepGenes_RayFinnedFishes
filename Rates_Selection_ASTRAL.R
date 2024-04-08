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
library("traitdataform")
library("traitdata")
library(janitor)
library(rfishbase)
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

VO_NVO_shape <-
  c("VO"=15,
    "NVO"=16)

##### NOTUNG results -- collapse 95 --------------------------------

#Duplications 
ASTRAL_duplications <- 
  read.table("batch.txt.AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel.duplication.txt",
             header=TRUE,
             sep="\t")
ASTRAL_duplications$GT.name <- sub(".prot.trimmed.aln.nooutgroup.collapsed95.treefile_recs.nwk.TREE", "", ASTRAL_duplications$GT.name)


ASTRAL_duplications_rev <-
  as.data.frame(t(ASTRAL_duplications))
ASTRAL_duplications_rev <- 
  ASTRAL_duplications_rev %>%
  row_to_names(row_number = 1)
ASTRAL_duplications_rev <- tibble::rownames_to_column(ASTRAL_duplications_rev, "label")
ASTRAL_duplications_rev <- 
  ASTRAL_duplications_rev %>%
  mutate(across(-c(label), as.numeric))



#Losses 
ASTRAL_losses <- 
  read.table("batch.txt.AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel.loss.txt",
             header=TRUE,
             sep="\t")
ASTRAL_losses$GT.name <- sub(".prot.trimmed.aln.nooutgroup.collapsed95.treefile_recs.nwk.TREE", "", ASTRAL_losses$GT.name)


ASTRAL_losses_rev <-
  as.data.frame(t(ASTRAL_losses))
ASTRAL_losses_rev <- 
  ASTRAL_losses_rev %>%
  row_to_names(row_number = 1)
ASTRAL_losses_rev <- tibble::rownames_to_column(ASTRAL_losses_rev, "label")
ASTRAL_losses_rev <- 
  ASTRAL_losses_rev %>%
  mutate(across(-c(label), as.numeric))



#gene Count 
ASTRAL_geneCount <- 
  read.table("batch.txt.AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel.geneCount.txt",
             header=TRUE,
             sep="\t")
ASTRAL_geneCount$GT.name <- sub(".prot.trimmed.aln.nooutgroup.collapsed95.treefile_recs.nwk.TREE", "", ASTRAL_geneCount$GT.name)


ASTRAL_geneCount_rev <-
  as.data.frame(t(ASTRAL_geneCount))
ASTRAL_geneCount_rev <- 
  ASTRAL_geneCount_rev %>%
  row_to_names(row_number = 1)
ASTRAL_geneCount_rev <- tibble::rownames_to_column(ASTRAL_geneCount_rev, "label")
ASTRAL_geneCount_rev <- 
  ASTRAL_geneCount_rev %>%
  mutate(across(-c(label), as.numeric))



#branch length
ASTRAL_branchlength <- 
  read.table("Species_branchlength.tsv",
             sep=" ",
             header=FALSE)
colnames(ASTRAL_branchlength) <- c("label", "length")

ASTRAL_branchlength <-
  ASTRAL_branchlength %>%
  mutate(real_length = length * 1000) %>% 
  dplyr::select(label, real_length)

colnames(ASTRAL_branchlength) <- c("label", "branch_length")


#Parse loss table
ASTRAL_losses_rev <- 
  ASTRAL_losses_rev %>%
  mutate(all_ospin_loss = 
           exorh+lws+opn3+opn4m1_3+opn4m2+opn4x+opn5+opn6+opn7a+opn7b+opn8a+
           opn8b+opn8c+opn9+parapinopsin+parietopsin+pinopsin+rgr+rh1+
           rh2+rrh+sws1+sws2+tmt1+tmt2+tmt3+va,
         VO_ospin_loss = 
           lws+rh1+rh2+sws1+sws2,
         NVO_ospin_loss = 
           exorh+opn3+opn4m1_3+opn4m2+opn4x+opn5+opn6+opn7a+opn7b+opn8a+
           opn8b+opn8c+opn9+parapinopsin+parietopsin+pinopsin+rgr+
           rrh+tmt1+tmt2+tmt3+va)

tail(ASTRAL_losses_rev %>% dplyr::select(label, all_ospin_loss) %>% arrange(all_ospin_loss), 10)


ASTRAL_geneCount_rev <- 
  ASTRAL_geneCount_rev %>%
  mutate(all_ospin_count = 
           exorh+lws+opn3+opn4m1_3+opn4m2+opn4x+opn5+opn6+opn7a+opn7b+opn8a+
           opn8b+opn8c+opn9+parapinopsin+parietopsin+pinopsin+rgr+rh1+
           rh2+rrh+sws1+sws2+tmt1+tmt2+tmt3+va,
         VO_ospin_count = 
           lws+rh1+rh2+sws1+sws2,
         NVO_ospin_count = 
           exorh+opn3+opn4m1_3+opn4m2+opn4x+opn5+opn6+opn7a+opn7b+opn8a+
           opn8b+opn8c+opn9+parapinopsin+parietopsin+pinopsin+rgr+
           rrh+tmt1+tmt2+tmt3+va)



ASTRAL_duplications_rev <- 
  ASTRAL_duplications_rev %>%
  mutate(all_ospin_dup = 
           exorh+lws+opn3+opn4m1_3+opn4m2+opn4x+opn5+opn6+opn7a+opn7b+opn8a+
           opn8b+opn8c+opn9+parapinopsin+parietopsin+pinopsin+rgr+rh1+
           rh2+rrh+sws1+sws2+tmt1+tmt2+tmt3+va,
         VO_ospin_dup = 
           lws+rh1+rh2+sws1+sws2,
         NVO_ospin_dup = 
           exorh+opn3+opn4m1_3+opn4m2+opn4x+opn5+opn6+opn7a+opn7b+opn8a+
           opn8b+opn8c+opn9+parapinopsin+parietopsin+pinopsin+rgr+
           rrh+tmt1+tmt2+tmt3+va)



ASTRAL_duplications_losses_rev <-
  left_join(ASTRAL_duplications_rev, ASTRAL_losses_rev, by="label", suffix = c(".dup", ".loss"))

ASTRAL_duplications_losses_count_rev <- 
  left_join(ASTRAL_duplications_losses_rev, ASTRAL_geneCount_rev, by="label")




#parse gene count df
ASTRAL_geneCount <-
  ASTRAL_geneCount %>%
  mutate(subfamily = GT.name) %>%
  dplyr::select(-GT.name)

ASTRAL_geneCount <- 
  as.data.frame(
    ASTRAL_geneCount %>%
      pivot_longer(!subfamily, 
                   names_to = "label", 
                   values_to = "Notung_count")
  )


ASTRAL_geneCount_total <- 
  ASTRAL_geneCount %>%
  group_by(label) %>%
  summarise(total_count = sum(Notung_count)) %>%
  mutate(subfamily = "Total") %>%
  dplyr::select(subfamily, label, total_count)
colnames(ASTRAL_geneCount_total) <- c("subfamily", "label", "Notung_count")
ASTRAL_geneCount_total <- as.data.frame(ASTRAL_geneCount_total)

ASTRAL_geneCount <- rbind(ASTRAL_geneCount, ASTRAL_geneCount_total)



#parse gene dup df
ASTRAL_duplications <-
  ASTRAL_duplications %>%
  mutate(subfamily = GT.name) %>%
  dplyr::select(-GT.name)

ASTRAL_duplications <- 
  as.data.frame(
    ASTRAL_duplications %>%
      pivot_longer(!subfamily, 
                   names_to = "label", 
                   values_to = "Notung_count")
  )


ASTRAL_duplications_total <- 
  ASTRAL_duplications %>%
  group_by(label) %>%
  summarise(total_count = sum(Notung_count)) %>%
  mutate(subfamily = "Total") %>%
  dplyr::select(subfamily, label, total_count)

colnames(ASTRAL_duplications_total) <- c("subfamily", "label", "duplication_nb")
ASTRAL_duplications_total <- as.data.frame(ASTRAL_duplications_total)


colnames(ASTRAL_duplications_total) <- c("subfamily", "label", "duplication_nb")
colnames(ASTRAL_duplications) <- c("subfamily", "label", "duplication_nb")


ASTRAL_duplications <- rbind(ASTRAL_duplications, ASTRAL_duplications_total)



#parse gene loss df
ASTRAL_losses <-
  ASTRAL_losses %>%
  mutate(subfamily = GT.name) %>%
  dplyr::select(-GT.name)

ASTRAL_losses <- 
  as.data.frame(
    ASTRAL_losses %>%
      pivot_longer(!subfamily, 
                   names_to = "label", 
                   values_to = "Notung_count")
  )


ASTRAL_losses_total <- 
  ASTRAL_losses %>%
  group_by(label) %>%
  summarise(total_count = sum(Notung_count)) %>%
  mutate(subfamily = "Total") %>%
  dplyr::select(subfamily, label, total_count)
colnames(ASTRAL_losses_total) <- c("subfamily", "label", "Notung_count")
ASTRAL_losses_total <- as.data.frame(ASTRAL_losses_total)


colnames(ASTRAL_losses_total) <- c("subfamily", "label", "loss_nb")
colnames(ASTRAL_losses) <- c("subfamily", "label", "loss_nb")

ASTRAL_losses <- rbind(ASTRAL_losses, ASTRAL_losses_total)



#Combine the three dataframes

Astral_duploss <- left_join(ASTRAL_losses, ASTRAL_duplications, by=c("subfamily", "label"))
Notung_Astral_df <- left_join(Astral_duploss, ASTRAL_geneCount, by=c("subfamily", "label"))

Notung_Astral_df_total <- 
  Notung_Astral_df %>%
  filter(subfamily == "Total")

#Now lets compute birth and death rates

Notung_Astral_df <- merge(ASTRAL_branchlength, Notung_Astral_df, by="label")

Notung_Astral_df <- Notung_Astral_df %>%
  mutate(Ancestral_gene_copy = Notung_count + loss_nb - duplication_nb) %>%
  mutate(Birth_rate = case_when(
    duplication_nb == 0 ~ 0,
    duplication_nb == loss_nb ~ (duplication_nb/(branch_length*Ancestral_gene_copy)),
    duplication_nb != loss_nb & duplication_nb != 0  ~ (duplication_nb/((duplication_nb-loss_nb)*branch_length)) * log(1+((duplication_nb-loss_nb)/(Ancestral_gene_copy)))
  )) %>%
  mutate(Death_rate = case_when(
    loss_nb == 0 ~ 0,
    duplication_nb != loss_nb & loss_nb == 1 & Ancestral_gene_copy == 1 ~ (loss_nb/(branch_length * Ancestral_gene_copy)),
    duplication_nb == loss_nb ~ (loss_nb/(branch_length*Ancestral_gene_copy)),
    duplication_nb != loss_nb & loss_nb != 0 & duplication_nb != 0 ~ (loss_nb/((duplication_nb-loss_nb)*branch_length)) * log(1+((duplication_nb-loss_nb)/(Ancestral_gene_copy))),
    duplication_nb != loss_nb & loss_nb != 0 & duplication_nb == 0 ~ (loss_nb/(branch_length * Ancestral_gene_copy))
    
  )) 

Notung_Astral_df <- 
  Notung_Astral_df %>%
  filter(Ancestral_gene_copy != 0)



Notung_Astral_df <- 
  Notung_Astral_df %>%
  mutate(Bi_bi = 
           Birth_rate * branch_length) %>%
  mutate(Li_bi =
           Death_rate * branch_length)





total_sum_Bi_bi <- sum(Notung_Astral_df %>%
                         filter(subfamily == "Total") %>%
                         pull(Bi_bi))

total_sum_Li_bi <- sum(Notung_Astral_df %>%
                         filter(subfamily == "Total") %>%
                         pull(Li_bi))

sum_branch_length <- sum(Notung_Astral_df %>% 
                           filter(subfamily == "Total") %>%
                           pull(branch_length))

Mean_Birth_rate <- total_sum_Bi_bi / sum_branch_length
Mean_Death_rate <- total_sum_Li_bi / sum_branch_length


Notung_Astral_df %>%
  filter(subfamily == "Total") %>%
  ggplot(., aes(x=Birth_rate)) +
  geom_density()


Notung_Astral_df %>%
  filter(subfamily != "Total") %>%
  ggplot(., aes(x=log(Birth_rate), color=subfamily)) +
  geom_density() 


Notung_Astral_mean_rates <- 
  as.data.frame(
    Notung_Astral_df %>%
      group_by(subfamily) %>%
      summarise(sum_Bi_bi = sum(Bi_bi),
                sum_Li_bi = sum(Li_bi),
                sum_branch_length = sum(branch_length)) %>%
      mutate(Mean_Birth_rate = sum_Bi_bi/sum_branch_length,
             Mean_Death_rate = sum_Li_bi/sum_branch_length))


Notung_Astral_mean_rates %>%
  filter(subfamily == "Total") %>%
  pull(Mean_Birth_rate) 


Notung_Astral_mean_rates %>%
  filter(subfamily == "Total") %>%
  pull(Mean_Death_rate) 



tail(
  Notung_Astral_df %>%
    filter(subfamily == "Total") %>%
    dplyr::select(label, duplication_nb, loss_nb) %>%
    arrange(loss_nb), 10)

tail(Notung_Astral_df %>%
       filter(subfamily == "Total") %>%
       dplyr::select(label, duplication_nb, loss_nb) %>%
       arrange(duplication_nb), 10)



#Notung_Astral_df %>%
#  filter(label == "Node298")


Notung_Astral_mean_rates %>%
  ggplot(., aes(x=subfamily, y=Mean_Birth_rate)) +
  geom_bar(stat="identity", fill="darkgreen")


Notung_Astral_mean_rates %>%
  ggplot(., aes(x=subfamily, y=Mean_Death_rate)) +
  geom_bar(stat="identity", fill="darkgreen")


Notung_Astral_mean_rates %>%
  ggplot(., aes(x=Mean_Birth_rate, y=Mean_Death_rate, color=subfamily)) +
  geom_point() +
  scale_color_manual(values = opsins_clades_colors) +
  theme_minimal()



Notung_Astral_mean_rates_collapse95 <- Notung_Astral_mean_rates



##### NOTUNG Summary --------------------------------



Notung_Astral_mean_rates_collapse95 <- 
  Notung_Astral_mean_rates_collapse95 %>%
  mutate(tree_parsing = "Collapsing 95")
Notung_Astral_mean_rates_all <-  Notung_Astral_mean_rates_collapse95


Notung_Astral_mean_rates_all %>%
  filter(subfamily == "Total")


Notung_Astral_mean_rates_all %>%
  filter(tree_parsing %in% c("Collapsing 95")) %>%
  arrange(Mean_Birth_rate) %>%
  dplyr::select(subfamily, Mean_Birth_rate, Mean_Death_rate)

Notung_Astral_mean_rates_all <- 
  Notung_Astral_mean_rates_all %>%
  mutate(VO_NVO = if_else(
    subfamily %in% c("rh2", "rh1", "lws", "sws1", "sws2"),
    "VO",
    "NVO"
  ))


Notung_Astral_mean_rates_all %>%
  filter(tree_parsing %in% c("Collapsing 95")) %>%
  filter(subfamily != "Total") %>%
  ggplot(., aes(x=Mean_Birth_rate, y=Mean_Death_rate, color=subfamily, shape=VO_NVO)) +
  geom_point(size = 3) +
  scale_color_manual(values = opsins_clades_colors) +
  scale_shape_manual(values = VO_NVO_shape) +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") 

Notung_Astral_mean_rates_all %>%
  filter(tree_parsing %in% c("Collapsing 95")) %>%
  filter(subfamily != "Total") %>%
  ggplot(., aes(x=Mean_Birth_rate, y=Mean_Death_rate, color=subfamily, shape=VO_NVO)) +
  geom_point(size = 3) +
  scale_color_manual(values = opsins_clades_colors) +
  scale_shape_manual(values = VO_NVO_shape) +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) 



Notung_Astral_mean_rates_all %>%
  filter(tree_parsing %in% c("Collapsing 95")) %>%
  filter(subfamily != "Total") %>%
  dplyr::select(subfamily, Mean_Birth_rate, Mean_Death_rate) %>%
  arrange(Mean_Birth_rate)



## Check if dynamic is linked to gene properties

Subfam_stats <- 
  read.table("Subfamilies_stats.csv",
             header=FALSE,
             sep=",")
colnames(Subfam_stats) <- c("subfamily", "mean_length", "mean_exon_nb")


Notung_Astral_mean_rates_all_famstat <- 
  left_join(Notung_Astral_mean_rates_all,
            Subfam_stats,
            by="subfamily")




cor_test_result <- 
  cor.test(
    Notung_Astral_mean_rates_all_famstat %>% 
      filter(tree_parsing %in% c("Collapsing 95")) %>%
      filter(subfamily != "Total") %>%
      pull(Mean_Birth_rate),
    Notung_Astral_mean_rates_all_famstat %>% 
      filter(tree_parsing %in% c("Collapsing 95")) %>%
      filter(subfamily != "Total") %>%
      pull(mean_length),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)


Notung_Astral_mean_rates_all_famstat %>%
  filter(tree_parsing %in% c("Collapsing 95")) %>%
  filter(subfamily != "Total") %>%
  ggplot(., aes(x=mean_length, y=Mean_Birth_rate, color=subfamily)) +
  geom_smooth(method = "lm", color = "black") +
  geom_point(size = 3 , aes(shape=VO_NVO)) +
  scale_color_manual(values = opsins_clades_colors) +
  scale_shape_manual(values = VO_NVO_shape) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  #annotate(x=400, y=0.003, 
  #         label=paste("R = ", estimate_p, "; P = ", pvalue_p), 
  #         geom="text", size=5) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  xlab("Mean gene size (AA)") +
  ylab("Mean birth rate") +
  theme(legend.position="none")
  








cor_test_result <- 
  cor.test(
    Notung_Astral_mean_rates_all_famstat %>% 
      filter(tree_parsing %in% c("Collapsing 95")) %>%
      filter(subfamily != "Total") %>%
      pull(Mean_Birth_rate),
    Notung_Astral_mean_rates_all_famstat %>% 
      filter(tree_parsing %in% c("Collapsing 95")) %>%
      filter(subfamily != "Total") %>%
      pull(mean_exon_nb),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)


Notung_Astral_mean_rates_all_famstat %>%
  filter(tree_parsing %in% c("Collapsing 95")) %>%
  filter(subfamily != "Total") %>%
  ggplot(., aes(x=mean_exon_nb, y=Mean_Birth_rate, color=subfamily)) +
  geom_smooth(method = "lm", color = "black") +
  geom_point(size = 3 , aes(shape=VO_NVO)) +
  scale_color_manual(values = opsins_clades_colors) +
  scale_shape_manual(values = VO_NVO_shape) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  #annotate(x=3, y=0.003, 
  #         label=paste("R = ", estimate_p, "; P = ", pvalue_p), 
  #         geom="text", size=5) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  xlab("Mean exon number") +
  ylab("Mean birth rate") +
  theme(legend.position="none")






cor_test_result <- 
  cor.test(
    Notung_Astral_mean_rates_all_famstat %>% 
      filter(tree_parsing %in% c("Collapsing 95")) %>%
      filter(subfamily != "Total") %>%
      pull(Mean_Death_rate),
    Notung_Astral_mean_rates_all_famstat %>% 
      filter(tree_parsing %in% c("Collapsing 95")) %>%
      filter(subfamily != "Total") %>%
      pull(mean_length),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)


Notung_Astral_mean_rates_all_famstat %>%
  filter(tree_parsing %in% c("Collapsing 95")) %>%
  filter(subfamily != "Total") %>%
  ggplot(., aes(x=mean_length, y=Mean_Death_rate, color=subfamily)) +
  geom_smooth(method = "lm", color = "black") +
  geom_point(size = 3 , aes(shape=VO_NVO)) +
  scale_color_manual(values = opsins_clades_colors) +
  scale_shape_manual(values = VO_NVO_shape) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  #annotate(x=350, y=0.01, 
  #         label=paste("R = ", estimate_p, "; P = ", pvalue_p), 
  #         geom="text", size=5) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  xlab("Mean gene size (AA)") +
  ylab("Mean death rate") +
  theme(legend.position="none")








cor_test_result <- 
  cor.test(
    Notung_Astral_mean_rates_all_famstat %>% 
      filter(tree_parsing %in% c("Collapsing 95")) %>%
      filter(subfamily != "Total") %>%
      pull(Mean_Death_rate),
    Notung_Astral_mean_rates_all_famstat %>% 
      filter(tree_parsing %in% c("Collapsing 95")) %>%
      filter(subfamily != "Total") %>%
      pull(mean_exon_nb),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)


Notung_Astral_mean_rates_all_famstat %>%
  filter(tree_parsing %in% c("Collapsing 95")) %>%
  filter(subfamily != "Total") %>%
  ggplot(., aes(x=mean_exon_nb, y=Mean_Death_rate, color=subfamily)) +
  geom_smooth(method = "lm", color = "black") +
  geom_point(size = 3 , aes(shape=VO_NVO)) +
  scale_color_manual(values = opsins_clades_colors) +
  scale_shape_manual(values = VO_NVO_shape) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  #annotate(x=3, y=0.01, 
  #         label=paste("R = ", estimate_p, "; P = ", pvalue_p), 
  #         geom="text", size=5) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  xlab("Mean exon number") +
  ylab("Mean death rate") +
  theme(legend.position="none")




##### NOTUNG results -- collapse 95 -- add VO / NVO --------------------------------

#Duplications 
ASTRAL_duplications <- 
  read.table("batch.txt.AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel.duplication.txt",
             header=TRUE,
             sep="\t")
ASTRAL_duplications$GT.name <- sub(".prot.trimmed.aln.nooutgroup.collapsed95.treefile_recs.nwk.TREE", "", ASTRAL_duplications$GT.name)


ASTRAL_duplications_rev <-
  as.data.frame(t(ASTRAL_duplications))
ASTRAL_duplications_rev <- 
  ASTRAL_duplications_rev %>%
  row_to_names(row_number = 1)
ASTRAL_duplications_rev <- tibble::rownames_to_column(ASTRAL_duplications_rev, "label")
ASTRAL_duplications_rev <- 
  ASTRAL_duplications_rev %>%
  mutate(across(-c(label), as.numeric))



#Losses 
ASTRAL_losses <- 
  read.table("batch.txt.AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel.loss.txt",
             header=TRUE,
             sep="\t")
ASTRAL_losses$GT.name <- sub(".prot.trimmed.aln.nooutgroup.collapsed95.treefile_recs.nwk.TREE", "", ASTRAL_losses$GT.name)


ASTRAL_losses_rev <-
  as.data.frame(t(ASTRAL_losses))
ASTRAL_losses_rev <- 
  ASTRAL_losses_rev %>%
  row_to_names(row_number = 1)
ASTRAL_losses_rev <- tibble::rownames_to_column(ASTRAL_losses_rev, "label")
ASTRAL_losses_rev <- 
  ASTRAL_losses_rev %>%
  mutate(across(-c(label), as.numeric))



#gene Count 
ASTRAL_geneCount <- 
  read.table("batch.txt.AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel.geneCount.txt",
             header=TRUE,
             sep="\t")
ASTRAL_geneCount$GT.name <- sub(".prot.trimmed.aln.nooutgroup.collapsed95.treefile_recs.nwk.TREE", "", ASTRAL_geneCount$GT.name)


ASTRAL_geneCount_rev <-
  as.data.frame(t(ASTRAL_geneCount))
ASTRAL_geneCount_rev <- 
  ASTRAL_geneCount_rev %>%
  row_to_names(row_number = 1)
ASTRAL_geneCount_rev <- tibble::rownames_to_column(ASTRAL_geneCount_rev, "label")
ASTRAL_geneCount_rev <- 
  ASTRAL_geneCount_rev %>%
  mutate(across(-c(label), as.numeric))



#branch length
ASTRAL_branchlength <- 
  read.table("Species_branchlength.tsv",
             sep=" ",
             header=FALSE)
colnames(ASTRAL_branchlength) <- c("label", "length")

ASTRAL_branchlength <-
  ASTRAL_branchlength %>%
  mutate(real_length = length * 1000) %>% 
  dplyr::select(label, real_length)

colnames(ASTRAL_branchlength) <- c("label", "branch_length")


#Parse loss table
ASTRAL_losses_rev <- 
  ASTRAL_losses_rev %>%
  mutate(all_ospin_loss = 
           exorh+lws+opn3+opn4m1_3+opn4m2+opn4x+opn5+opn6+opn7a+opn7b+opn8a+
           opn8b+opn8c+opn9+parapinopsin+parietopsin+pinopsin+rgr+rh1+
           rh2+rrh+sws1+sws2+tmt1+tmt2+tmt3+va,
         VO_ospin_loss = 
           lws+rh1+rh2+sws1+sws2,
         NVO_ospin_loss = 
           exorh+opn3+opn4m1_3+opn4m2+opn4x+opn5+opn6+opn7a+opn7b+opn8a+
           opn8b+opn8c+opn9+parapinopsin+parietopsin+pinopsin+rgr+
           rrh+tmt1+tmt2+tmt3+va)

tail(ASTRAL_losses_rev %>% dplyr::select(label, all_ospin_loss) %>% arrange(all_ospin_loss), 10)


ASTRAL_geneCount_rev <- 
  ASTRAL_geneCount_rev %>%
  mutate(all_ospin_count = 
           exorh+lws+opn3+opn4m1_3+opn4m2+opn4x+opn5+opn6+opn7a+opn7b+opn8a+
           opn8b+opn8c+opn9+parapinopsin+parietopsin+pinopsin+rgr+rh1+
           rh2+rrh+sws1+sws2+tmt1+tmt2+tmt3+va,
         VO_ospin_count = 
           lws+rh1+rh2+sws1+sws2,
         NVO_ospin_count = 
           exorh+opn3+opn4m1_3+opn4m2+opn4x+opn5+opn6+opn7a+opn7b+opn8a+
           opn8b+opn8c+opn9+parapinopsin+parietopsin+pinopsin+rgr+
           rrh+tmt1+tmt2+tmt3+va)



ASTRAL_duplications_rev <- 
  ASTRAL_duplications_rev %>%
  mutate(all_ospin_dup = 
           exorh+lws+opn3+opn4m1_3+opn4m2+opn4x+opn5+opn6+opn7a+opn7b+opn8a+
           opn8b+opn8c+opn9+parapinopsin+parietopsin+pinopsin+rgr+rh1+
           rh2+rrh+sws1+sws2+tmt1+tmt2+tmt3+va,
         VO_ospin_dup = 
           lws+rh1+rh2+sws1+sws2,
         NVO_ospin_dup = 
           exorh+opn3+opn4m1_3+opn4m2+opn4x+opn5+opn6+opn7a+opn7b+opn8a+
           opn8b+opn8c+opn9+parapinopsin+parietopsin+pinopsin+rgr+
           rrh+tmt1+tmt2+tmt3+va)



ASTRAL_duplications_losses_rev <-
  left_join(ASTRAL_duplications_rev, ASTRAL_losses_rev, by="label", suffix = c(".dup", ".loss"))

ASTRAL_duplications_losses_count_rev <- 
  left_join(ASTRAL_duplications_losses_rev, ASTRAL_geneCount_rev, by="label")




#parse gene count df
ASTRAL_geneCount <-
  ASTRAL_geneCount %>%
  mutate(subfamily = GT.name) %>%
  dplyr::select(-GT.name)

ASTRAL_geneCount <- 
  as.data.frame(
    ASTRAL_geneCount %>%
      pivot_longer(!subfamily, 
                   names_to = "label", 
                   values_to = "Notung_count")
  )


ASTRAL_geneCount_total <- 
  ASTRAL_geneCount %>%
  group_by(label) %>%
  summarise(total_count = sum(Notung_count)) %>%
  mutate(subfamily = "Total") %>%
  dplyr::select(subfamily, label, total_count)
colnames(ASTRAL_geneCount_total) <- c("subfamily", "label", "Notung_count")
ASTRAL_geneCount_total <- as.data.frame(ASTRAL_geneCount_total)



ASTRAL_geneCount_Visual <- 
  ASTRAL_geneCount %>%
  filter(subfamily %in% c("lws","rh1","rh2","sws1","sws2")) %>%
  group_by(label) %>%
  summarise(total_count = sum(Notung_count)) %>%
  mutate(subfamily = "Visual") %>%
  dplyr::select(subfamily, label, total_count)
colnames(ASTRAL_geneCount_Visual) <- c("subfamily", "label", "Notung_count")
ASTRAL_geneCount_Visual <- as.data.frame(ASTRAL_geneCount_Visual)



ASTRAL_geneCount_Non_Visual <- 
  ASTRAL_geneCount %>%
  filter(! subfamily %in% c("lws","rh1","rh2","sws1","sws2")) %>%
  group_by(label) %>%
  summarise(total_count = sum(Notung_count)) %>%
  mutate(subfamily = "Non_visual") %>%
  dplyr::select(subfamily, label, total_count)
colnames(ASTRAL_geneCount_Non_Visual) <- c("subfamily", "label", "Notung_count")
ASTRAL_geneCount_Non_Visual <- as.data.frame(ASTRAL_geneCount_Non_Visual)




ASTRAL_geneCount <- rbind(ASTRAL_geneCount, ASTRAL_geneCount_total)
ASTRAL_geneCount <- rbind(ASTRAL_geneCount, ASTRAL_geneCount_Visual)
ASTRAL_geneCount <- rbind(ASTRAL_geneCount, ASTRAL_geneCount_Non_Visual)



#parse gene dup df
ASTRAL_duplications <-
  ASTRAL_duplications %>%
  mutate(subfamily = GT.name) %>%
  dplyr::select(-GT.name)

ASTRAL_duplications <- 
  as.data.frame(
    ASTRAL_duplications %>%
      pivot_longer(!subfamily, 
                   names_to = "label", 
                   values_to = "Notung_count")
  )


ASTRAL_duplications_total <- 
  ASTRAL_duplications %>%
  group_by(label) %>%
  summarise(total_count = sum(Notung_count)) %>%
  mutate(subfamily = "Total") %>%
  dplyr::select(subfamily, label, total_count)
colnames(ASTRAL_duplications_total) <- c("subfamily", "label", "duplication_nb")
ASTRAL_duplications_total <- as.data.frame(ASTRAL_duplications_total)


ASTRAL_duplications_Visual <- 
  ASTRAL_duplications %>%
  filter(subfamily %in% c("lws","rh1","rh2","sws1","sws2")) %>%
  group_by(label) %>%
  summarise(total_count = sum(Notung_count)) %>%
  mutate(subfamily = "Visual") %>%
  dplyr::select(subfamily, label, total_count)
colnames(ASTRAL_duplications_Visual) <- c("subfamily", "label", "duplication_nb")
ASTRAL_duplications_Visual <- as.data.frame(ASTRAL_duplications_Visual)


ASTRAL_duplications_Non_Visual <- 
  ASTRAL_duplications %>%
  filter(! subfamily %in% c("lws","rh1","rh2","sws1","sws2")) %>%
  group_by(label) %>%
  summarise(total_count = sum(Notung_count)) %>%
  mutate(subfamily = "Non_visual") %>%
  dplyr::select(subfamily, label, total_count)
colnames(ASTRAL_duplications_Non_Visual) <- c("subfamily", "label", "duplication_nb")
ASTRAL_duplications_Non_Visual <- as.data.frame(ASTRAL_duplications_Non_Visual)



colnames(ASTRAL_duplications_total) <- c("subfamily", "label", "duplication_nb")
colnames(ASTRAL_duplications) <- c("subfamily", "label", "duplication_nb")
colnames(ASTRAL_duplications_Visual) <- c("subfamily", "label", "duplication_nb")
colnames(ASTRAL_duplications_Non_Visual) <- c("subfamily", "label", "duplication_nb")


ASTRAL_duplications <- rbind(ASTRAL_duplications, ASTRAL_duplications_total)
ASTRAL_duplications <- rbind(ASTRAL_duplications, ASTRAL_duplications_Visual)
ASTRAL_duplications <- rbind(ASTRAL_duplications, ASTRAL_duplications_Non_Visual)


#parse gene loss df
ASTRAL_losses <-
  ASTRAL_losses %>%
  mutate(subfamily = GT.name) %>%
  dplyr::select(-GT.name)

ASTRAL_losses <- 
  as.data.frame(
    ASTRAL_losses %>%
      pivot_longer(!subfamily, 
                   names_to = "label", 
                   values_to = "Notung_count")
  )


ASTRAL_losses_total <- 
  ASTRAL_losses %>%
  group_by(label) %>%
  summarise(total_count = sum(Notung_count)) %>%
  mutate(subfamily = "Total") %>%
  dplyr::select(subfamily, label, total_count)
colnames(ASTRAL_losses_total) <- c("subfamily", "label", "Notung_count")
ASTRAL_losses_total <- as.data.frame(ASTRAL_losses_total)



ASTRAL_losses_Visual <- 
  ASTRAL_losses %>%
  filter(subfamily %in% c("lws","rh1","rh2","sws1","sws2")) %>%
  group_by(label) %>%
  summarise(total_count = sum(Notung_count)) %>%
  mutate(subfamily = "Visual") %>%
  dplyr::select(subfamily, label, total_count)
colnames(ASTRAL_losses_Visual) <- c("subfamily", "label", "duplication_nb")
ASTRAL_losses_Visual <- as.data.frame(ASTRAL_losses_Visual)


ASTRAL_losses_Non_Visual <- 
  ASTRAL_losses %>%
  filter(! subfamily %in% c("lws","rh1","rh2","sws1","sws2")) %>%
  group_by(label) %>%
  summarise(total_count = sum(Notung_count)) %>%
  mutate(subfamily = "Non_visual") %>%
  dplyr::select(subfamily, label, total_count)
colnames(ASTRAL_losses_Non_Visual) <- c("subfamily", "label", "duplication_nb")
ASTRAL_losses_Non_Visual <- as.data.frame(ASTRAL_losses_Non_Visual)




colnames(ASTRAL_losses_total) <- c("subfamily", "label", "loss_nb")
colnames(ASTRAL_losses) <- c("subfamily", "label", "loss_nb")
colnames(ASTRAL_losses_Visual) <- c("subfamily", "label", "loss_nb")
colnames(ASTRAL_losses_Non_Visual) <- c("subfamily", "label", "loss_nb")

ASTRAL_losses <- rbind(ASTRAL_losses, ASTRAL_losses_total)
ASTRAL_losses <- rbind(ASTRAL_losses, ASTRAL_losses_Visual)
ASTRAL_losses <- rbind(ASTRAL_losses, ASTRAL_losses_Non_Visual)



#Combine the three dataframes

Astral_duploss <- left_join(ASTRAL_losses, ASTRAL_duplications, by=c("subfamily", "label"))
Notung_Astral_df <- left_join(Astral_duploss, ASTRAL_geneCount, by=c("subfamily", "label"))

Notung_Astral_df_total <- 
  Notung_Astral_df %>%
  filter(subfamily == "Total")

#Now lets compute birth and death rates

Notung_Astral_df <- merge(ASTRAL_branchlength, Notung_Astral_df, by="label")

Notung_Astral_df <- Notung_Astral_df %>%
  mutate(Ancestral_gene_copy = Notung_count + loss_nb - duplication_nb) %>%
  mutate(Birth_rate = case_when(
    duplication_nb == 0 ~ 0,
    duplication_nb == loss_nb ~ (duplication_nb/(branch_length*Ancestral_gene_copy)),
    duplication_nb != loss_nb & duplication_nb != 0  ~ (duplication_nb/((duplication_nb-loss_nb)*branch_length)) * log(1+((duplication_nb-loss_nb)/(Ancestral_gene_copy)))
  )) %>%
  mutate(Death_rate = case_when(
    loss_nb == 0 ~ 0,
    duplication_nb != loss_nb & loss_nb == 1 & Ancestral_gene_copy == 1 ~ (loss_nb/(branch_length * Ancestral_gene_copy)),
    duplication_nb == loss_nb ~ (loss_nb/(branch_length*Ancestral_gene_copy)),
    duplication_nb != loss_nb & loss_nb != 0 & duplication_nb != 0 ~ (loss_nb/((duplication_nb-loss_nb)*branch_length)) * log(1+((duplication_nb-loss_nb)/(Ancestral_gene_copy))),
    duplication_nb != loss_nb & loss_nb != 0 & duplication_nb == 0 ~ (loss_nb/(branch_length * Ancestral_gene_copy))
    
  )) 

Notung_Astral_df <- 
  Notung_Astral_df %>%
  filter(Ancestral_gene_copy != 0)



Notung_Astral_df <- 
  Notung_Astral_df %>%
  mutate(Bi_bi = 
           Birth_rate * branch_length) %>%
  mutate(Li_bi =
           Death_rate * branch_length)





total_sum_Bi_bi <- sum(Notung_Astral_df %>%
                         filter(subfamily == "Total") %>%
                         pull(Bi_bi))

total_sum_Li_bi <- sum(Notung_Astral_df %>%
                         filter(subfamily == "Total") %>%
                         pull(Li_bi))

sum_branch_length <- sum(Notung_Astral_df %>% 
                           filter(subfamily == "Total") %>%
                           pull(branch_length))

Mean_Birth_rate <- total_sum_Bi_bi / sum_branch_length
Mean_Death_rate <- total_sum_Li_bi / sum_branch_length


Notung_Astral_df %>%
  filter(subfamily == "Total") %>%
  ggplot(., aes(x=Birth_rate)) +
  geom_density()


Notung_Astral_df %>%
  filter(subfamily != "Total") %>%
  ggplot(., aes(x=log(Birth_rate), color=subfamily)) +
  geom_density() 


Notung_Astral_mean_rates <- 
  as.data.frame(
    Notung_Astral_df %>%
      group_by(subfamily) %>%
      summarise(sum_Bi_bi = sum(Bi_bi),
                sum_Li_bi = sum(Li_bi),
                sum_branch_length = sum(branch_length)) %>%
      mutate(Mean_Birth_rate = sum_Bi_bi/sum_branch_length,
             Mean_Death_rate = sum_Li_bi/sum_branch_length))


Notung_Astral_mean_rates %>%
  filter(subfamily == "Total") %>%
  pull(Mean_Birth_rate) 


Notung_Astral_mean_rates %>%
  filter(subfamily == "Total") %>%
  pull(Mean_Death_rate) 



tail(
  Notung_Astral_df %>%
    filter(subfamily == "Total") %>%
    dplyr::select(label, duplication_nb, loss_nb) %>%
    arrange(loss_nb), 10)

tail(Notung_Astral_df %>%
       filter(subfamily == "Total") %>%
       dplyr::select(label, duplication_nb, loss_nb) %>%
       arrange(duplication_nb), 10)



#Notung_Astral_df %>%
#  filter(label == "Node298")


Notung_Astral_mean_rates %>%
  ggplot(., aes(x=subfamily, y=Mean_Birth_rate)) +
  geom_bar(stat="identity", fill="darkgreen")


Notung_Astral_mean_rates %>%
  ggplot(., aes(x=subfamily, y=Mean_Death_rate)) +
  geom_bar(stat="identity", fill="darkgreen")


Notung_Astral_mean_rates %>%
  ggplot(., aes(x=Mean_Birth_rate, y=Mean_Death_rate, color=subfamily)) +
  geom_point() +
  scale_color_manual(values = opsins_clades_colors) +
  theme_minimal()


#wilcox.test(
#  Notung_Astral_df %>% filter(subfamily == "lws") %>% pull(Birth_rate),
#  Notung_Astral_df %>% filter(subfamily == "rh2") %>% pull(Birth_rate)
#)
Notung_Astral_mean_rates_collapse95 <- Notung_Astral_mean_rates




##### Correlation dup / loss nb ?  --------------------------------

## First lets check duplications 

Notung_Astral_df_parsed_duplications <- 
  Notung_Astral_df %>%
  filter(branch_length > 2) %>%
  filter(subfamily %in% c("Visual", "Non_visual")) %>%
  dplyr::select(label, subfamily, Birth_rate)

Notung_Astral_df_parsed_losses <- 
  Notung_Astral_df %>%
  filter(branch_length > 2) %>%
  filter(subfamily %in% c("Visual", "Non_visual")) %>%
  dplyr::select(label, subfamily, Death_rate)


Notung_Astral_dup_wide <- 
  Notung_Astral_df_parsed_duplications %>%
  pivot_wider(names_from = subfamily, values_from = Birth_rate)
Notung_Astral_loss_wide <- 
  Notung_Astral_df_parsed_losses %>%
  pivot_wider(names_from = subfamily, values_from = Death_rate)




cor_test_result <- 
  cor.test(
    Notung_Astral_dup_wide %>% pull(Visual),
    Notung_Astral_dup_wide %>% pull(Non_visual),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)


Notung_Astral_dup_wide %>%
  ggplot(., aes(x=Non_visual, y=Visual)) +
  geom_smooth(method = "lm", color = "black") +
  geom_point()




Notung_Astral_dup_wide %>%
  ggplot(., aes(x=Non_visual, y=Visual)) +
  geom_smooth(method = "lm", color = "black") +
  geom_point() +
  scale_color_manual(values = opsins_clades_colors) +
  theme_minimal() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  annotate(x=3, y=0.003, 
           label=paste("R = ", estimate_p, "; P = ", pvalue_p), 
           geom="text", size=5) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  xlab("Proportion of positive sites") +
  ylab("Mean birth rate")



## Now lets check losses 



##### Dup vs Sec -- BUSTED --------------------------------

#Load pvalues from BUSTED results for each subfamily

BUSTED_pvalues <- 
  read.table(
    "BUSTED/Table_subfam_bustedPvalue.csv",
    sep=",",
    header=FALSE)
colnames(BUSTED_pvalues) <- c("subfamily", "pvalue")


#Extract dN/dS from BUSTED results

all_genes <- BUSTED_pvalues %>% pull(subfamily)

ALL_subfam_BUSTED_df <- as.data.frame(NULL)
for(curr_gene in all_genes){
  

  setwd(paste0("/Users/maximepolicarpo/Desktop/Manuscript_Opsins_RayFinnedFishes.nosync/FigShare_Files/R_scripts_related_files/Rates_Selection_ASTRAL/Reconciled_trees_95/",curr_gene, "/"))

  dup_omega <- 
    as.numeric(
      system("grep -A3 'Improving branch lengths' BUSTED/slurm* | grep 'test' | sed 's/.* //g'", intern = TRUE))
    
  spec_omega <- 
    as.numeric(
      system("grep -A3 'Improving branch lengths' BUSTED/slurm* | grep 'background' | sed 's/.* //g'", intern = TRUE))
  
  
  #Make a summary table
  
  curr_df <- 
    as.data.frame(cbind(curr_gene, spec_omega, dup_omega))
  colnames(curr_df) <- 
    c("subfamily","Spec_omega", "Dup_omega")
  

  curr_df$Spec_omega <- as.numeric(curr_df$Spec_omega)
  curr_df$Dup_omega <- as.numeric(curr_df$Dup_omega)
  
  ALL_subfam_BUSTED_df <- 
    rbind(ALL_subfam_BUSTED_df, curr_df)
  
} 


BUSTED_df <- 
  left_join(BUSTED_pvalues, ALL_subfam_BUSTED_df, by="subfamily")




BUSTED_df <- 
  BUSTED_df %>%
  mutate(significant = if_else(
    pvalue <= 0.05,
    "Significant",
    "Non-significant"
  ))



signif_shape <-
  c("Significant"=17,
    "Non-significant"=19)



BUSTED_df %>%
  ggplot(., aes(x=Spec_omega, y=Dup_omega, color=subfamily, shape=significant)) +
  geom_point(size = 3) + 
  geom_abline(slope=1, linetype="dashed", color = "gray") +
  scale_color_manual(values = opsins_clades_colors) +
  scale_shape_manual(values = signif_shape) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  xlim(0, 0.7) +
  ylim(0, 0.7) + 
  theme(legend.position="none")

BUSTED_df %>%
  ggplot(., aes(x=Spec_omega, y=Dup_omega, color=subfamily, shape=significant)) +
  geom_point(size = 3) + 
  geom_abline(slope=1, linetype="dashed", color = "gray") +
  scale_color_manual(values = opsins_clades_colors) +
  scale_shape_manual(values = signif_shape) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  xlim(0, 0.7) +
  ylim(0, 0.7) 




# Any link with birth and death rates ? 

BUSTED_df <- 
  BUSTED_df %>%
  mutate(delta_omega = Dup_omega - Spec_omega)

BUSTED_df_rates <- 
  left_join(BUSTED_df, Notung_Astral_mean_rates_all_famstat, by="subfamily")



## Birth rate 


cor_test_result <- 
  cor.test(
    BUSTED_df_rates %>% 
      filter(tree_parsing %in% c("Collapsing 95")) %>%
      filter(subfamily != "Total") %>%
      pull(Mean_Birth_rate),
    BUSTED_df_rates %>% 
      filter(tree_parsing %in% c("Collapsing 95")) %>%
      filter(subfamily != "Total") %>%
      pull(delta_omega),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)



BUSTED_df_rates %>%
  filter(tree_parsing == "Collapsing 95") %>%
  ggplot(., aes(x=Mean_Birth_rate, y=delta_omega, color=subfamily)) +
  geom_point(size = 3 , aes(shape = VO_NVO)) + 
  scale_color_manual(values = opsins_clades_colors) +
  scale_shape_manual(values = VO_NVO_shape) + 
  geom_smooth(method = "lm", color = "black") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none") 


## Death rate 


cor_test_result <- 
  cor.test(
    BUSTED_df_rates %>% 
      filter(tree_parsing %in% c("Collapsing 95")) %>%
      filter(subfamily != "Total") %>%
      pull(Mean_Death_rate),
    BUSTED_df_rates %>% 
      filter(tree_parsing %in% c("Collapsing 95")) %>%
      filter(subfamily != "Total") %>%
      pull(delta_omega),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)



BUSTED_df_rates %>%
  filter(tree_parsing == "Collapsing 95") %>%
  ggplot(., aes(x=Mean_Death_rate, y=delta_omega, color=subfamily)) +
  geom_point(size = 3 , aes(shape = VO_NVO)) + 
  scale_color_manual(values = opsins_clades_colors) +
  scale_shape_manual(values = VO_NVO_shape) + 
  geom_smooth(method = "lm", color = "black") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none")


##### Omega per branch -- Data load  --------------------------------

omega_sp_tree_df <- 
  read.table("OMEGA_per_branch_speciestree/Opsins.omega.csv",
             header=FALSE,
             sep=",")
colnames(omega_sp_tree_df) <- c("LB", "MLE", "UB", "subfamily", "label", "dN", "dS")
omega_gene_tree_df <- 
  read.table("OMEGA_per_branch_genetree/Opsins.omega.csv",
             header=FALSE,
             sep=",")
colnames(omega_gene_tree_df) <- c("LB", "MLE", "UB", "subfamily", "label", "dN", "dS")


species_tree_corr_df <- 
  read.table("OMEGA_per_branch_speciestree/Corresponding_labels.csv",
             header=FALSE,
             sep=",")
colnames(species_tree_corr_df) <- c("subfamily", "old_label", "label")

omega_sp_tree_df <- 
  omega_sp_tree_df %>%
  mutate(VO_NVO = if_else(
    subfamily %in% c("lws", "sws1", "sws2", "rh1", "rh2"),
    "visual",
    "non-visual"
  ))

omega_gene_tree_df <- 
  omega_gene_tree_df %>%
  mutate(VO_NVO = if_else(
    subfamily %in% c("lws", "sws1", "sws2", "rh1", "rh2"),
    "visual",
    "non-visual"
  ))



##### Omega per branch -- Comparison species/gene tree  --------------------------------

### Compute the median omega per subfamily

median_omega_df_sp <- 
  as.data.frame(
    omega_sp_tree_df %>% 
      filter(dS > 0.01) %>% filter(dS < 1) %>% filter(dN < 1) %>%
      #filter(MLE < 2.5) %>%
      #filter(LB != 0) %>%
      #filter(UB != 10000) %>%
      group_by(subfamily) %>%
      summarise(median_omega = median(MLE),
                mean_omega = mean(MLE))
  )


median_omega_df_gene <- 
  as.data.frame(
    omega_gene_tree_df %>% 
      filter(dS > 0.01) %>% filter(dS < 1) %>% filter(dN < 1) %>%
      #filter(MLE < 2.5) %>%
      #filter(LB != 0) %>%
      #filter(UB != 10000) %>%
      group_by(subfamily) %>%
      summarise(median_omega = median(MLE),
                mean_omega = mean(MLE))
  )


### Compare results between gene tree and species tree


omega_df_both <- 
  left_join(median_omega_df_sp, median_omega_df_gene, 
            by="subfamily",
            suffix = c(".sp", ".gene"))


cor_test_result <- 
  cor.test(
    omega_df_both %>% 
      pull(mean_omega.gene),
    omega_df_both %>% 
      pull(mean_omega.sp),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)



omega_df_both %>%
  ggplot(., aes(x=mean_omega.sp, y=mean_omega.gene, color=subfamily)) +
  geom_point() + 
  geom_abline(slope=1, linetype="dashed", color = "gray") +
  geom_smooth(method = "lm", color = "black") +
  scale_color_manual(values = opsins_clades_colors) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none")





### We can compare also the w on terminal branches

## First we have to rename the branch labels for the species tree

omega_sp_tree_df_newnames <-
  left_join(species_tree_corr_df, omega_sp_tree_df, by=c("label", "subfamily")) %>%
  dplyr::select(subfamily, old_label, LB, MLE, UB, dN, dS, VO_NVO)
colnames(omega_sp_tree_df_newnames) <-
  c("subfamily", "label", "LB", "MLE", "UB","dN", "dS", "VO_NVO")

omega_sp_tree_df_newnames$label <- gsub("\\.", "_", omega_sp_tree_df_newnames$label)
omega_sp_tree_df_newnames$label <- gsub("---", "_", omega_sp_tree_df_newnames$label)
omega_sp_tree_df_newnames$label <- gsub("-", "_", omega_sp_tree_df_newnames$label)


omega_sp_tree_vs_gene_tree <- 
  left_join(omega_sp_tree_df_newnames, omega_gene_tree_df, by="label", suffix = c(".sptree", ".genetree"))

omega_sp_tree_vs_gene_tree_OKbranches <-
  omega_sp_tree_vs_gene_tree %>%
  filter(! is.na(MLE.genetree)) %>%
  filter(! is.na(MLE.sptree)) %>%
  #filter(MLE.genetree < 2.5) %>%
  #filter(LB.genetree != 0) %>%
  #filter(UB.genetree != 10000) %>%
  #filter(MLE.sptree < 2.5) %>%
  #filter(LB.sptree != 0) %>%
  #filter(UB.sptree != 10000) 
  filter(dS.genetree > 0.01) %>% filter(dS.genetree < 1) %>% filter(dN.genetree < 1) %>%
  filter(dS.sptree > 0.01) %>% filter(dS.sptree < 1) %>% filter(dN.sptree < 1) 
  

cor_test_result <- 
  cor.test(
    omega_sp_tree_vs_gene_tree_OKbranches %>% 
      pull(MLE.genetree),
    omega_sp_tree_vs_gene_tree_OKbranches %>% 
      pull(MLE.sptree),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)




omega_sp_tree_vs_gene_tree_OKbranches %>%
  ggplot(., aes(x=MLE.sptree, y=MLE.genetree, color=subfamily.sptree)) +
  geom_point() + 
  geom_abline(slope=1, linetype="dashed", color = "gray") +
  geom_smooth(method = "lm", color = "black") +
  scale_color_manual(values = opsins_clades_colors) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none") +
  xlab("dN/dS - Reconciled tree") +
  ylab("dN/dS - ML gene tree") +
  ylim(0, 3) +
  xlim(0, 3)






##### Omega per branch -- Gene level  --------------------------------


## First make a boxplot with the w distribution per subfamily


omega_sp_tree_df %>%
  #filter(MLE < 2.5) %>%
  #filter(LB != 0) %>%
  #filter(UB != 10000) %>%
  filter(dS > 0.01) %>%
  filter(dS < 1) %>%
  filter(dN < 1) %>%
  ggplot(., aes(x=reorder(subfamily, MLE, mean), y=MLE, color=subfamily)) +
  geom_boxplot() +
  scale_color_manual(values = opsins_clades_colors) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none") +
  ylab("dN/dS") +
  ylim(0, 2.5)



## Now make a boxplot of the w distribution in VO vs NVO ... 



omega_sp_tree_df %>%
  #filter(MLE < 2.5) %>%
  #filter(LB != 0) %>%
  #filter(UB != 10000) %>%
  filter(dS > 0.01) %>%
  filter(dS < 1) %>%
  filter(dN < 1) %>%
  ggplot(., aes(x=VO_NVO, y=MLE)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none") +
  ylab("dN/dS") +
  ylim(0, 2.5)


  
  

all_wilcox_test_subfams <- 
  pairwise.wilcox.test(
    x=omega_sp_tree_df %>% filter(dS > 0.01) %>% filter(dS < 1) %>% filter(dN < 1) %>% pull(MLE),
    g=omega_sp_tree_df %>% filter(dS > 0.01) %>% filter(dS < 1) %>% filter(dN < 1) %>% pull(subfamily),
    p.adjust.method="bonferroni"
  )

my_pvals <- all_wilcox_test_subfams$p.value

my_pvals[my_pvals > 0.05 & !is.na(my_pvals)] <- NA


my_transformed_pvals=-log10(my_pvals)






#corrplot(as.matrix(my_transformed_pvals),is.corr=F, type = 'lower', tl.col="black")

library(ggcorrplot)

ggcorrplot(as.matrix(my_transformed_pvals)) +
  scale_fill_gradient2(low = "#6D9EC1", high =  "darkred", mid = "white", midpoint = 1)



#Wiloxon test between non visual and visual opsins



wilcox.test(
  omega_sp_tree_df %>% 
    filter(dS > 0.01) %>%  filter(dS < 1) %>% filter(dN < 1) %>%
    filter(VO_NVO == "visual") %>%
    pull(MLE),
  omega_sp_tree_df %>% 
    filter(dS > 0.01) %>%  filter(dS < 1) %>% filter(dN < 1) %>%
    filter(VO_NVO == "non-visual") %>%
    pull(MLE)
)






omega_df_both %>%
  arrange(mean_omega.sp) %>%
  dplyr::select(subfamily, median_omega.sp, mean_omega.sp)
  
  
  
  
  
##### Omega per branch -- Link with birth / death rates  --------------------------------

Notung_Astral_mean_rates_omega <- 
  left_join(Notung_Astral_mean_rates_all,
            omega_df_both,
            by="subfamily")

Notung_Astral_mean_rates_omega_95 <- 
  Notung_Astral_mean_rates_omega %>%
  filter(subfamily != "Total") %>%
  filter(tree_parsing == "Collapsing 95")



## First with Birth Rate


cor_test_result <- 
  cor.test(
    Notung_Astral_mean_rates_omega_95 %>% 
      pull(Mean_Birth_rate),
    Notung_Astral_mean_rates_omega_95 %>% 
      pull(mean_omega.sp),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)





Notung_Astral_mean_rates_omega_95 %>%
  ggplot(., aes(x=mean_omega.sp, y=Mean_Birth_rate, color=subfamily)) +
  geom_point(size = 3, aes(shape = VO_NVO)) + 
  geom_smooth(method = "lm", color = "black") +
  scale_color_manual(values = opsins_clades_colors) +
  scale_shape_manual(values = VO_NVO_shape) + 
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none") +
  xlab("Mean omega") +
  ylab("Mean Birth Rate")

Notung_Astral_mean_rates_omega_95 %>%
  ggplot(., aes(x=median_omega.sp, y=Mean_Birth_rate, color=subfamily)) +
  geom_point(size = 3, aes(shape = VO_NVO)) +
  geom_smooth(method = "lm", color = "black") +
  scale_color_manual(values = opsins_clades_colors) +
  scale_shape_manual(values = VO_NVO_shape) + 
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none") +
  xlab("Median omega") +
  ylab("Mean Birth Rate")



## Now with Death Rate

cor_test_result <- 
  cor.test(
    Notung_Astral_mean_rates_omega_95 %>% 
      pull(Mean_Death_rate),
    Notung_Astral_mean_rates_omega_95 %>% 
      pull(mean_omega.sp),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)





Notung_Astral_mean_rates_omega_95 %>%
  ggplot(., aes(x=mean_omega.sp, y=Mean_Birth_rate, color=subfamily)) +
  geom_point(size = 3, aes(shape = VO_NVO)) +
  geom_smooth(method = "lm", color = "black") +
  scale_color_manual(values = opsins_clades_colors) +
  scale_shape_manual(values = VO_NVO_shape) + 
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none") +
  xlab("Mean omega") +
  ylab("Mean Death Rate")


Notung_Astral_mean_rates_omega_95 %>%
  ggplot(., aes(x=median_omega.sp, y=Mean_Birth_rate, color=subfamily)) +
  geom_point(size = 3, aes(shape = VO_NVO)) +
  geom_smooth(method = "lm", color = "black") +
  scale_color_manual(values = opsins_clades_colors) +
  scale_shape_manual(values = VO_NVO_shape) + 
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none") +
  xlab("Median omega") +
  ylab("Mean Death Rate")







##### Omega per branch -- Species level  --------------------------------

omega_sp_tree_df_newnames <-
  left_join(species_tree_corr_df, omega_sp_tree_df, by=c("label", "subfamily")) %>%
  dplyr::select(subfamily, old_label, LB, MLE, UB, dN, dS, VO_NVO)
colnames(omega_sp_tree_df_newnames) <-
  c("subfamily", "label", "LB", "MLE", "UB", "dN", "dS","VO_NVO")

omega_sp_tree_df_newnames$species <- gsub("---.*", "", omega_sp_tree_df_newnames$label)
  



omega_sp_tree_df_newnames %>%
  #filter(MLE < 2.5) %>%
  #filter(LB != 0) %>%
  #filter(UB != 10000) %>%
  filter(dS > 0.01) %>%
  filter(dS < 1) %>%
  filter(dN < 1) %>%
  ggplot(., aes(x=reorder(species, MLE, mean), y=MLE)) +
  geom_boxplot() +
  #scale_color_manual(values = opsins_clades_colors) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(legend.position="none") +
  ylab("dN/dS") +
  xlab("species") +
  ylim(0, 3)



#Wiloxon test between all species



all_wilcox_test_species  <- 
  pairwise.wilcox.test(
    x=omega_sp_tree_df_newnames %>% filter(dS > 0.01) %>% filter(dS < 1) %>% filter(dN < 1) %>% pull(MLE),
    g=omega_sp_tree_df_newnames %>% filter(dS > 0.01) %>% filter(dS < 1) %>% filter(dN < 1) %>% pull(species),
    p.adjust.method="bonferroni"
  )


## check some stats

omega_sp_tree_df_newnames_OKvalues <-
  omega_sp_tree_df_newnames %>% 
  filter(dS > 0.01) %>%
  filter(dS < 1) %>%
  filter(dN < 1) 


omega_per_sp_summary <- 
  as.data.frame(
    omega_sp_tree_df_newnames_OKvalues %>%
  group_by(species) %>%
  summarise(median_omega = median(MLE),
            mean_omega = mean(MLE),
            nb_genes = n())
  ) %>%
  filter(nb_genes >= 3)


tail(omega_per_sp_summary %>% arrange(mean_omega), 20)
head(omega_per_sp_summary %>% arrange(mean_omega), 20)



omega_per_sp_summary %>%
  ggplot(., aes(x=median_omega)) +
  geom_histogram(bins=80, color="black", fill="gray") +
  theme_classic() +
  xlab("median dN/dS") +
  ylab("Number of species") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none") 


omega_per_sp_summary %>%
  ggplot(., aes(x=mean_omega)) +
  geom_histogram(bins=80, color="black", fill="gray") +
  theme_classic() +
  xlab("mean dN/dS") +
  ylab("Number of species") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none") 



cor_test_result <- 
  cor.test(
    omega_per_sp_summary %>% 
      pull(nb_genes),
    omega_per_sp_summary %>% 
      pull(median_omega),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)


omega_per_sp_summary %>%
  ggplot(., aes(x=median_omega, y=nb_genes)) +
  geom_point() +
  geom_smooth(method = "lm", color = "black") +
  theme_classic() +
  xlab("median dN/dS") +
  ylab("Number of genes") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none") 



cor_test_result <- 
  cor.test(
    omega_per_sp_summary %>% 
      pull(nb_genes),
    omega_per_sp_summary %>% 
      pull(mean_omega),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)


omega_per_sp_summary %>%
  ggplot(., aes(x=mean_omega, y=nb_genes)) +
  geom_point() +
  geom_smooth(method = "lm", color = "black") +
  theme_classic() +
  xlab("mean dN/dS") +
  ylab("Number of genes") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none")





omega_per_sp_summary_info <-
  left_join(omega_per_sp_summary, species_table, by='species')

cor_test_result <- 
  cor.test(
    omega_per_sp_summary_info %>% 
      filter(infraclass == "Teleostei") %>%
      filter(WGD == "No") %>%
      pull(nb_genes),
    omega_per_sp_summary_info %>% 
      filter(infraclass == "Teleostei") %>%
      filter(WGD == "No") %>%
      pull(median_omega),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)


omega_per_sp_summary_info %>%
  ggplot(., aes(x=median_omega, y=nb_genes)) +
  geom_point() +
  geom_smooth(method = "lm", color = "black") +
  theme_classic() +
  xlab("median dN/dS") +
  ylab("Number of genes") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none") 



cor_test_result <- 
  cor.test(
    omega_per_sp_summary_info %>% 
      pull(nb_genes),
    omega_per_sp_summary_info %>% 
      pull(mean_omega),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)


omega_per_sp_summary_info %>%
  ggplot(., aes(x=mean_omega, y=nb_genes)) +
  geom_point() +
  geom_smooth(method = "lm", color = "black") +
  theme_classic() +
  xlab("mean dN/dS") +
  ylab("Number of genes") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none")



## Compute the median or mean omega per subfamily

omega_per_sp_per_gene <- 
  as.data.frame(
    omega_sp_tree_df_newnames_OKvalues %>%
      group_by(species, subfamily) %>%
      summarise(median_omega = median(MLE),
                mean_omega = mean(MLE))
  ) 

omega_per_sp_per_gene_median <- 
  as.data.frame(
    omega_per_sp_per_gene %>%
  pivot_wider(names_from = subfamily, values_from = c(median_omega, mean_omega)))

## Compute the median or mean omega per VO/VNO


omega_per_sp_summary <- 
  as.data.frame(
    omega_sp_tree_df_newnames_OKvalues %>%
      group_by(species) %>%
      summarise(median_omega = median(MLE),
                mean_omega = mean(MLE),
                nb_genes = n())
  ) %>%
  filter(nb_genes >= 3)



omega_per_sp_summary_NVO <- 
  as.data.frame(
    omega_sp_tree_df_newnames_OKvalues %>%
      filter(VO_NVO == "non-visual") %>%
      group_by(species) %>%
      summarise(median_omega_non_visual = median(MLE),
                mean_omega_non_visual = mean(MLE),
                nb_genes = n())
  ) %>%
  filter(nb_genes >= 3)

omega_per_sp_summary_VO <- 
  as.data.frame(
    omega_sp_tree_df_newnames_OKvalues %>%
      filter(VO_NVO == "visual") %>%
      group_by(species) %>%
      summarise(median_omega_visual = median(MLE),
                mean_omega_visual = mean(MLE),
                nb_genes = n())
  ) %>%
  filter(nb_genes >= 3)


omega_per_sp_summary %>%
  ggplot(., aes(x=mean_omega, y=nb_genes)) +
  geom_point() +
  geom_smooth(method = "lm", color = "black") +
  theme_classic() +
  xlab("mean dN/dS") +
  ylab("Number of genes") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none")



##### pGLS between NVO and VO omega values   -------------------

teleost_noWGD <- 
  all_ospins_complete_wide_ecology %>%
  filter(infraclass == "Teleostei") %>%
  filter(WGD == "No") %>%
  pull(species)

teleost_WGD <- 
  all_ospins_complete_wide_ecology %>%
  filter(infraclass == "Teleostei") %>%
  filter(WGD == "Yes") %>%
  pull(species)



omega_per_sp_summary_VO_NVO <- 
  left_join(omega_per_sp_summary_NVO, omega_per_sp_summary_VO, by="species")

curr_astral_caper <- 
  comparative.data(phy = ASTRAL_species_tree_nodelabel, data = omega_per_sp_summary_VO_NVO,
                   names.col = species, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)

curr_astral_caper_onlyteleost <- 
  comparative.data(phy = ASTRAL_species_tree_nodelabel, 
                   data = omega_per_sp_summary_VO_NVO %>% filter(species %in% teleost_noWGD),
                   names.col = species, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)



pgls_omega <-
  pgls(mean_omega_visual ~ mean_omega_non_visual, 
       data = curr_astral_caper, 
       lambda = "ML",
       bounds=list(lambda=c(0,1)))

pgls_omega_teleostNoWGD <-
  pgls(mean_omega_visual ~ mean_omega_non_visual, 
       data = curr_astral_caper_onlyteleost, 
       lambda = "ML",
       bounds=list(lambda=c(0,1)))

summary(pgls_omega)
summary(pgls_omega_teleostNoWGD)

curr_lm <- 
  lm(data = omega_per_sp_summary_VO_NVO %>% 
       filter(! is.na(mean_omega_non_visual)) %>%
       filter(! is.na(mean_omega_visual)),
     formula = mean_omega_visual ~ mean_omega_non_visual)
curr_function <- GLS_function(curr_lm)


curr_lm_teleost <- 
  lm(data = omega_per_sp_summary_VO_NVO %>% 
       filter(! is.na(mean_omega_non_visual)) %>%
       filter(! is.na(mean_omega_visual)) %>% 
       filter(species %in% teleost_noWGD),
     formula = mean_omega_visual ~ mean_omega_non_visual)
curr_function_teleost <- GLS_function(curr_lm_teleost)




omega_per_sp_summary_VO_NVO_info <- 
  omega_per_sp_summary_VO_NVO %>%
  mutate(teleost_WGD = if_else(
    species %in% teleost_WGD,
    "Teleost_WGD",
    "Other"
  ))





omega_per_sp_summary_VO_NVO_info %>%
  ggplot(., aes(x=mean_omega_non_visual, y=mean_omega_visual, color=teleost_WGD)) +
  geom_point() +
  theme_classic() +
  stat_function(fun = curr_function, color="gray") +
  stat_function(fun = curr_function_teleost, color="black") +
  #labs(subtitle = bquote(R^2 ~ "=" ~ .(curr_R2) ~ ";" ~ FDR.pvalue ~ "=" ~.(curr_pvalue) ~ ";" ~ lambda ~ "=" ~ .(curr_lambda))) +
  xlab("Mean visual dN/dS") +
  ylab("Mean non-visual dN/dS") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") +
  scale_color_manual(values = c("Teleost_WGD" = "gray", "Other" = "black"))



nrow(omega_per_sp_summary_VO_NVO %>% 
       filter(! is.na(mean_omega_non_visual)) %>% 
       filter(! is.na(mean_omega_visual))
)


nrow(omega_per_sp_summary_VO_NVO %>% 
       filter(! is.na(mean_omega_non_visual)) %>% 
       filter(! is.na(mean_omega_visual)) %>%
       filter(species %in% teleost_noWGD)
)

##### pGLS between ecology and omega   -------------------


all_ospins_complete_wide_ecology <- 
  read.table("all_ospins_complete_wide_ecology.tsv",
             header=TRUE,
             sep="\t")
  
  
  
all_ospins_complete_wide_omega <- 
  left_join(omega_per_sp_summary, all_ospins_complete_wide_ecology, by="species")

all_ospins_complete_wide_omega <- 
  left_join(omega_per_sp_per_gene_median, all_ospins_complete_wide_omega, by="species")

all_ospins_complete_wide_omega <- 
  left_join(all_ospins_complete_wide_omega, omega_per_sp_summary_VO_NVO, by="species")


responses_variables <- 
  c("mean_omega", colnames(omega_per_sp_per_gene_median)[29:55], colnames(omega_per_sp_summary_VO_NVO)[c(3,6)])


categorical_predictor_variables <- 
  c("Fresh_salt", "depth_category", "Food_cat", "Diel_pattern_simplified",
    "Electrogenic_simplified", "EnvTemp", "AnaCat")

numeric_pred_variables <- 
  c("ED_SL","Absolute_ED", "Absolute_SL")


ASTRAL_omega_pgls_cat <- as.data.frame(NULL)
ASTRAL_omega_pgls_num <- as.data.frame(NULL)
for(curr_response in responses_variables){
  for(curr_predictor in categorical_predictor_variables){
    
    
    current_test_name <- paste(curr_response, curr_predictor, sep="_")
    print(current_test_name)
    
    curr_data_opsins <- 
      all_ospins_complete_wide_omega %>%
      dplyr::select(species, curr_response, curr_predictor) %>%
      filter(! is.na(curr_predictor))
    colnames(curr_data_opsins) <- c("species", "reponse", "predictor")
    
    
    
    #remove categories with 2 or less species ...
    
    good_preds <- 
      curr_data_opsins %>%
      group_by(predictor) %>%
      summarise(count = n()) %>%
      filter(count > 2) %>%
      pull(predictor)
    
    curr_data_opsins <-
      curr_data_opsins %>%
      filter(predictor %in% good_preds)
    
    
    nb_factors <- 
      nrow(curr_data_opsins %>% 
             filter(! is.na(reponse)) %>% 
             group_by(predictor) %>% 
             summarise(n()) %>%
             filter(! is.na(predictor)))
    
    curr_astral_caper <- 
      comparative.data(phy = ASTRAL_species_tree_nodelabel, data = curr_data_opsins,
                       names.col = species, vcv = TRUE,
                       na.omit = FALSE, warn.dropped = TRUE)
    
    
    
    if(nb_factors >  1){
      curr_pgls_result <- "none"
    }
    
    if(nb_factors >=  2){
      if(current_test_name %in% c("mean_omega_tmt2_Fresh_salt","mean_omega_rgr_Electrogenic_simplified","mean_omega_rgr_DemersPelag","mean_omega_opn4m2_EnvTemp","median_omega_rgr_Electrogenic_simplified","median_omega_rgr_DemersPelag","median_omega_rh2_EnvTemp","median_omega_opn7b_EnvTemp","median_omega_opn4m2_EnvTemp")){ #reduce the lambda search boundary for failed computations
        curr_pgls_result <-
          pgls(reponse ~ predictor, 
               data = curr_astral_caper, 
               lambda = "ML",
               bounds=list(lambda=c(0.5,1)))
      } else {
        curr_pgls_result <-
          pgls(reponse ~ predictor, 
               data = curr_astral_caper, 
               lambda = "ML",
               bounds=list(lambda=c(0,1)))
      }
      
      
      
      
      
      curr_pgls_summary <- summary(curr_pgls_result)
      
      pGLS_R2 <- as.numeric(curr_pgls_summary$adj.r.squared)
      pGLS_lambda <- as.numeric(curr_pgls_summary$param[2])
      
      pGLS_pvalue <- pf(curr_pgls_summary$fstatistic[1], 
                        curr_pgls_summary$fstatistic[2], 
                        curr_pgls_summary$fstatistic[3], 
                        lower.tail = FALSE)
      
      
      curr_df <- 
        as.data.frame(
          cbind(
            curr_response,
            pGLS_R2,
            pGLS_pvalue,
            pGLS_lambda,
            curr_predictor
          )
        )
      
      colnames(curr_df) <- c("Response", "R2", "pvalue", "lambda","predictor")
      
      
      ASTRAL_omega_pgls_cat <- 
        rbind(ASTRAL_omega_pgls_cat,
              curr_df)
      
    }
  }
  
  
  
  for(curr_predictor_num in numeric_pred_variables){
    
    current_test_name <- paste(curr_response, curr_predictor_num, sep="_")
    print(current_test_name)
    
    curr_data_opsins <- 
      all_ospins_complete_wide_omega %>%
      dplyr::select(species, curr_response, curr_predictor_num)
    
    
    colnames(curr_data_opsins) <- c("species", "reponse", "predictor")
    
    curr_astral_caper <- 
      comparative.data(phy = ASTRAL_species_tree_nodelabel, data = curr_data_opsins,
                       names.col = species, vcv = TRUE,
                       na.omit = FALSE, warn.dropped = TRUE)
    
    
    if(current_test_name %in% c("mean_omega_tmt2_Fresh_salt","mean_omega_rgr_Electrogenic_simplified","mean_omega_rgr_DemersPelag","mean_omega_opn4m2_EnvTemp","median_omega_rgr_Electrogenic_simplified","median_omega_rgr_DemersPelag","median_omega_rh2_EnvTemp","median_omega_opn7b_EnvTemp","median_omega_opn4m2_EnvTemp")){ #reduce the lambda search boundary for failed computations
      curr_pgls_result <-
        pgls(reponse ~ predictor, 
             data = curr_astral_caper, 
             lambda = "ML",
             bounds=list(lambda=c(0.7,1)))
    } else {
      curr_pgls_result <-
        pgls(reponse ~ predictor, 
             data = curr_astral_caper, 
             lambda = "ML",
             bounds=list(lambda=c(0,1)))
    }
    
    sum_fit_phy <- summary(curr_pgls_result)
    PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
    PGLS_r2 =      formatC(sum_fit_phy$adj.r.squared, digits = 2)
    PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
    if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
    lambda <- as.numeric(sum_fit_phy$param[2])
    slope <- formatC(sum_fit_phy$coefficients[2], digits = 3)
    
    print(PGLS_pvalue)
    
    curr_df <- 
      as.data.frame(
        cbind(
          curr_response,
          PGLS_r2,
          PGLS_pvalue,
          lambda,
          slope,
          curr_predictor_num
        )
      )
    
    colnames(curr_df) <- c("Response", "R2", "pvalue", "lambda", "slope","predictor")
    
    ASTRAL_omega_pgls_num <- 
      rbind(ASTRAL_omega_pgls_num, curr_df)
  }
  
}




curr_astral_caper <- 
  comparative.data(phy = ASTRAL_species_tree_nodelabel, data = all_ospins_complete_wide_omega,
                   names.col = species, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)
omega_vs_WGD <- 
  pgls(mean_omega ~ WGD, 
     data = curr_astral_caper, 
     lambda = "ML",
     bounds=list(lambda=c(0,1)))

summary(omega_vs_WGD)




all_ospins_complete_wide_omega %>%
  filter(! is.na(WGD)) %>%
  ggplot(., aes(y=mean_omega, x=WGD)) +
  geom_boxplot() +
  theme_classic() +
  xlab("") +
  ylab("mean dN/dS") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") 


all_ospins_complete_wide_omega %>%
  arrange(desc(mean_omega)) %>%
  dplyr::select(species, mean_omega , WGD)





##### pGLS between ecology and omega -- teleost no WGD  -------------------

all_ospins_complete_wide_omega <- 
  left_join(omega_per_sp_summary, all_ospins_complete_wide_ecology, by="species")

all_ospins_complete_wide_omega <- 
  left_join(omega_per_sp_per_gene_median, all_ospins_complete_wide_omega, by="species")

all_ospins_complete_wide_omega <- 
  left_join(all_ospins_complete_wide_omega, omega_per_sp_summary_VO_NVO, by="species")



responses_variables <- 
  c("mean_omega", colnames(omega_per_sp_per_gene_median)[29:55], colnames(omega_per_sp_summary_VO_NVO)[c(3,6)])
responses_variables <- 
  responses_variables[! responses_variables %in% c("mean_omega_pinopsin")]


categorical_predictor_variables <- 
  c("Fresh_salt", "depth_category", "Food_cat", "Diel_pattern_simplified",
    "Electrogenic_simplified", "EnvTemp", "AnaCat")

numeric_pred_variables <- 
  c("ED_SL","Absolute_ED", "Absolute_SL")


ASTRAL_omega_pgls_teleost_cat <- as.data.frame(NULL)
ASTRAL_omega_pgls_teleost_num <- as.data.frame(NULL)
for(curr_response in responses_variables){
  for(curr_predictor in categorical_predictor_variables){
    
    
    current_test_name <- paste(curr_response, curr_predictor, sep="_")
    print(current_test_name)
    
    curr_data_opsins <- 
      all_ospins_complete_wide_omega %>%
      filter(infraclass == "Teleostei") %>%
      filter(WGD == "No") %>%
      dplyr::select(species, curr_response, curr_predictor) %>%
      filter(! is.na(curr_predictor))
    colnames(curr_data_opsins) <- c("species", "reponse", "predictor")
    
    
    
    #remove categories with 2 or less species ...
    
    good_preds <- 
      curr_data_opsins %>%
      group_by(predictor) %>%
      summarise(count = n()) %>%
      filter(count > 2) %>%
      pull(predictor)
    
    curr_data_opsins <-
      curr_data_opsins %>%
      filter(predictor %in% good_preds)
    
    
    nb_factors <- 
      nrow(curr_data_opsins %>% 
             filter(! is.na(reponse)) %>% 
             group_by(predictor) %>% 
             summarise(n()) %>%
             filter(! is.na(predictor)))
    
    curr_astral_caper <- 
      comparative.data(phy = ASTRAL_species_tree_nodelabel, data = curr_data_opsins,
                       names.col = species, vcv = TRUE,
                       na.omit = FALSE, warn.dropped = TRUE)
    
    
    
    if(nb_factors >  1){
      curr_pgls_result <- "none"
    }
    
    if(nb_factors >=  2){
      if(current_test_name %in% c("mean_omega_rgr_Electrogenic_simplified","mean_omega_rgr_DemersPelag","mean_omega_opn4m2_EnvTemp","median_omega_rgr_Electrogenic_simplified","median_omega_rgr_DemersPelag","median_omega_rh2_EnvTemp","median_omega_opn7b_EnvTemp","median_omega_opn4m2_EnvTemp")){ #reduce the lambda search boundary for failed computations
        curr_pgls_result <-
          pgls(reponse ~ predictor, 
               data = curr_astral_caper, 
               lambda = "ML",
               bounds=list(lambda=c(0.5,1)))
      } else {
        curr_pgls_result <-
          pgls(reponse ~ predictor, 
               data = curr_astral_caper, 
               lambda = "ML",
               bounds=list(lambda=c(0,1)))
      }
      
      
      
      
      
      curr_pgls_summary <- summary(curr_pgls_result)
      
      pGLS_R2 <- as.numeric(curr_pgls_summary$adj.r.squared)
      pGLS_lambda <- as.numeric(curr_pgls_summary$param[2])
      
      pGLS_pvalue <- pf(curr_pgls_summary$fstatistic[1], 
                        curr_pgls_summary$fstatistic[2], 
                        curr_pgls_summary$fstatistic[3], 
                        lower.tail = FALSE)
      
      
      curr_df <- 
        as.data.frame(
          cbind(
            curr_response,
            pGLS_R2,
            pGLS_pvalue,
            pGLS_lambda,
            curr_predictor
          )
        )
      
      colnames(curr_df) <- c("Response", "R2", "pvalue", "lambda","predictor")
      
      
      ASTRAL_omega_pgls_teleost_cat <- 
        rbind(ASTRAL_omega_pgls_teleost_cat,
              curr_df)
      
    }
  }
  
  
  
  for(curr_predictor_num in numeric_pred_variables){
    
    current_test_name <- paste(curr_response, curr_predictor_num, sep="_")
    print(current_test_name)
    
    curr_data_opsins <- 
      all_ospins_complete_wide_omega %>%
      filter(infraclass == "Teleostei") %>%
      filter(WGD == "No") %>%
      dplyr::select(species, curr_response, curr_predictor_num)
    
    
    colnames(curr_data_opsins) <- c("species", "reponse", "predictor")
    
    curr_astral_caper <- 
      comparative.data(phy = ASTRAL_species_tree_nodelabel, data = curr_data_opsins,
                       names.col = species, vcv = TRUE,
                       na.omit = FALSE, warn.dropped = TRUE)
    
    
    if(current_test_name %in% c("mean_omega_rgr_Electrogenic_simplified","mean_omega_rgr_DemersPelag","mean_omega_opn4m2_EnvTemp","median_omega_rgr_Electrogenic_simplified","median_omega_rgr_DemersPelag","median_omega_rh2_EnvTemp","median_omega_opn7b_EnvTemp","median_omega_opn4m2_EnvTemp")){ #reduce the lambda search boundary for failed computations
      curr_pgls_result <-
        pgls(reponse ~ predictor, 
             data = curr_astral_caper, 
             lambda = "ML",
             bounds=list(lambda=c(0.7,1)))
    } else {
      curr_pgls_result <-
        pgls(reponse ~ predictor, 
             data = curr_astral_caper, 
             lambda = "ML",
             bounds=list(lambda=c(0,1)))
    }
    
    sum_fit_phy <- summary(curr_pgls_result)
    PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
    PGLS_r2 =      formatC(sum_fit_phy$adj.r.squared, digits = 2)
    PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
    if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
    lambda <- as.numeric(sum_fit_phy$param[2])
    slope <- formatC(sum_fit_phy$coefficients[2], digits = 3)
    
    print(PGLS_pvalue)
    
    curr_df <- 
      as.data.frame(
        cbind(
          curr_response,
          PGLS_r2,
          PGLS_pvalue,
          lambda,
          slope,
          curr_predictor_num
        )
      )
    
    colnames(curr_df) <- c("Response", "R2", "pvalue", "lambda", "slope","predictor")
    
    ASTRAL_omega_pgls_teleost_num <- 
      rbind(ASTRAL_omega_pgls_teleost_num, curr_df)
  }
  
}






#####  OMEGA vs Ecology Summary and analysis   -------------------

#Import tables

ASTRAL_omega_pgls_num <- 
  read.table("ASTRAL_omega_pgls_num",
             sep="\t",
             header=TRUE)


ASTRAL_omega_pgls_cat <- 
  read.table("ASTRAL_omega_pgls_cat", 
             sep="\t",
             header=TRUE)


ASTRAL_omega_pgls_cat <- 
  ASTRAL_omega_pgls_cat %>%
  mutate(slope = NA)


ASTRAL_omega_pgls_cat <- 
  ASTRAL_omega_pgls_cat %>%
  dplyr::select(Response,R2, pvalue, lambda, slope, predictor)
colnames(ASTRAL_omega_pgls_cat) <- colnames(ASTRAL_omega_pgls_num)




#Merge tables


astral_omega_ecology <- 
  rbind(ASTRAL_omega_pgls_num, ASTRAL_omega_pgls_cat)


astral_omega_ecology <- 
  astral_omega_ecology %>%
  filter(! predictor %in% c("Amphibious", "Terrestrial_exploration", "Diadromous_cat"))


as.data.frame(
  astral_omega_ecology %>%
    group_by(Response) %>%
    summarise(n()))

as.data.frame(
  astral_omega_ecology %>%
    group_by(predictor) %>%
    summarise(n()))


#Add the FDR and Bonferroni p-value 

unique_response <- astral_omega_ecology$Response %>% unique()
astral_ecology_omega_corr <- as.data.frame(NULL)
for(curr_response in unique_response){
  
  curr_df <- 
    astral_omega_ecology %>%
    filter(Response == curr_response)
  
  
  initial_pvalues <- 
    as.numeric(
      curr_df %>% pull(pvalue)
    )
  
  corrected_pvalues_FDR <- p.adjust(initial_pvalues, method = "fdr")
  corrected_pvalues_bonferroni <- p.adjust(initial_pvalues, method = "bonferroni")
  
  curr_df <- 
    curr_df %>%
    mutate(FDR_pvalue = corrected_pvalues_FDR,
           Bonferroni_pvalue = corrected_pvalues_bonferroni)
  
  astral_ecology_omega_corr <- 
    rbind(astral_ecology_omega_corr, curr_df)
  
  
}




astral_ecology_omega_corr <- 
  astral_ecology_omega_corr %>%
  mutate(FDR_significant = if_else(
    as.numeric(FDR_pvalue) < 0.05,
    "Significant", 
    "N.S"
  ))

astral_ecology_omega_corr <- 
  astral_ecology_omega_corr %>%
  mutate(Bonferroni_significant = if_else(
    as.numeric(Bonferroni_pvalue) < 0.05,
    "Significant", 
    "N.S"
  ))

astral_ecology_omega_corr <- 
  astral_ecology_omega_corr %>%
  mutate(significant = if_else(
    as.numeric(pvalue) < 0.05,
    "Significant", 
    "N.S"
  ))


astral_ecology_omega_corr <- 
  astral_ecology_omega_corr %>%
  mutate(slope_sign = case_when(
    slope <= 0 ~ "minus",
    slope > 0 ~ "plus",
    is.na(slope) ~ "cat",
  ))



slope_sign_shapes <- 
  c("minus"="\u25BC",
    "plus"="\u25B2",
    "cat"="\u25FE")

#Lets make a summary graphic 


##Fill missing corr because of not enough data...

uniq_resp <- astral_ecology_omega_corr$Response %>% unique()
uniq_pred <- astral_ecology_omega_corr$predictor %>% unique()

for(curr_resp in uniq_resp){
  for(curr_pred in uniq_pred){
    
    nb_row <- 
      nrow(astral_ecology_omega_corr %>%
             filter(Response == curr_resp) %>%
             filter(predictor == curr_pred))
    
    if(nb_row < 1){
      
      new_line <- 
        as.data.frame(
          cbind(curr_resp, 0.1, 0.1, 0.1, 1, curr_pred, 1, 1, "N.S", "N.S", "N.S", "minus"))
      colnames(new_line) <- colnames(astral_ecology_omega_corr)
      astral_ecology_omega_corr <- rbind(astral_ecology_omega_corr, new_line)
      
    }
    
  }
}


astral_ecology_omega_corr$R2 <- as.numeric(astral_ecology_omega_corr$R2)


astral_ecology_omega_corr[(astral_ecology_omega_corr$Response == 'mean_omega_pinopsin'),"R2"] <- rep(0.1, 10)




astral_ecology_omega_corr$Response <-
  factor(astral_ecology_omega_corr$Response ,
         levels=c("mean_omega_rh1", "mean_omega_exorh", "mean_omega_rh2", "mean_omega_sws2", 
                  "mean_omega_sws1", "mean_omega_lws", "mean_omega_pinopsin", "mean_omega_va", 
                  "mean_omega_parapinopsin","mean_omega_parietopsin", "mean_omega_tmt1", "mean_omega_tmt2", 
                  "mean_omega_tmt3","mean_omega_opn3", "mean_omega_opn8b", "mean_omega_opn8a", 
                  "mean_omega_opn8c", "mean_omega_opn7b", "mean_omega_opn7a", "mean_omega_opn5", 
                  "mean_omega_opn9", "mean_omega_opn6", "mean_omega_rgr", "mean_omega_rrh", "mean_omega_opn4m1_3", 
                  "mean_omega_opn4m2", "mean_omega_opn4x", "mean_omega_non_visual", "mean_omega_visual", "mean_omega"))


astral_ecology_omega_corr$predictor <-
  factor(astral_ecology_omega_corr$predictor ,
         levels=c("Absolute_SL", "Absolute_ED", 
                  "ED_SL", "Electrogenic_simplified", "Food_cat","AnaCat", 
                  "Diel_pattern_simplified",
                  "EnvTemp", "Fresh_salt",
                  "depth_category"
         ))



Signif_colors <-
  c("N.S" = 0,
    "Significant" = 1)




ggplot(astral_ecology_omega_corr, 
       aes(y=predictor, x=Response, fill= R2, alpha=Bonferroni_significant, color="black")) + 
  geom_tile(color = "black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis_c(option = "plasma") +
  scale_alpha_manual(values = Signif_colors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("")

ggplot(astral_ecology_omega_corr, 
       aes(y=predictor, x=Response, fill= R2, alpha=Bonferroni_significant, color="black")) + 
  geom_tile(color = "black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_minimal() +
  scale_fill_viridis_c(option = "plasma") +
  scale_alpha_manual(values = Signif_colors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("")




#####  OMEGA vs Ecology Summary and analysis -- Teleost NO WGD  -------------------

#Import tables

ASTRAL_omega_pgls_teleost_num <- 
  read.table("ASTRAL_omega_pgls_teleost_num",
             sep="\t",
             header=TRUE)


ASTRAL_omega_pgls_teleost_cat <- 
  read.table("ASTRAL_omega_pgls_teleost_cat", 
             sep="\t",
             header=TRUE)


ASTRAL_omega_pgls_teleost_cat <- 
  ASTRAL_omega_pgls_teleost_cat %>%
  mutate(slope = NA)


ASTRAL_omega_pgls_teleost_cat <- 
  ASTRAL_omega_pgls_teleost_cat %>%
  dplyr::select(Response,R2, pvalue, lambda, slope, predictor)
colnames(ASTRAL_omega_pgls_teleost_cat) <- colnames(ASTRAL_omega_pgls_teleost_num)




#Merge tables


astral_omega_ecology_teleost <- 
  rbind(ASTRAL_omega_pgls_teleost_num, ASTRAL_omega_pgls_teleost_cat)


astral_omega_ecology_teleost <- 
  astral_omega_ecology_teleost %>%
  filter(! predictor %in% c("Amphibious", "Terrestrial_exploration", "Diadromous_cat"))


as.data.frame(
  astral_omega_ecology_teleost %>%
    group_by(Response) %>%
    summarise(n()))

as.data.frame(
  astral_omega_ecology_teleost %>%
    group_by(predictor) %>%
    summarise(n()))


#Add the FDR and Bonferroni p-value 

unique_response <- astral_omega_ecology_teleost$Response %>% unique()
astral_ecology_omega_corr_teleost <- as.data.frame(NULL)
for(curr_response in unique_response){
  
  curr_df <- 
    astral_omega_ecology_teleost %>%
    filter(Response == curr_response)
  
  
  initial_pvalues <- 
    as.numeric(
      curr_df %>% pull(pvalue)
    )
  
  corrected_pvalues_FDR <- p.adjust(initial_pvalues, method = "fdr")
  corrected_pvalues_bonferroni <- p.adjust(initial_pvalues, method = "bonferroni")
  
  curr_df <- 
    curr_df %>%
    mutate(FDR_pvalue = corrected_pvalues_FDR,
           Bonferroni_pvalue = corrected_pvalues_bonferroni)
  
  astral_ecology_omega_corr_teleost <- 
    rbind(astral_ecology_omega_corr_teleost, curr_df)
  
  
}




astral_ecology_omega_corr_teleost <- 
  astral_ecology_omega_corr_teleost %>%
  mutate(FDR_significant = if_else(
    as.numeric(FDR_pvalue) < 0.05,
    "Significant", 
    "N.S"
  ))

astral_ecology_omega_corr_teleost <- 
  astral_ecology_omega_corr_teleost %>%
  mutate(Bonferroni_significant = if_else(
    as.numeric(Bonferroni_pvalue) < 0.05,
    "Significant", 
    "N.S"
  ))

astral_ecology_omega_corr_teleost <- 
  astral_ecology_omega_corr_teleost %>%
  mutate(significant = if_else(
    as.numeric(pvalue) < 0.05,
    "Significant", 
    "N.S"
  ))


astral_ecology_omega_corr_teleost <- 
  astral_ecology_omega_corr_teleost %>%
  mutate(slope_sign = case_when(
    slope <= 0 ~ "minus",
    slope > 0 ~ "plus",
    is.na(slope) ~ "cat",
  ))



slope_sign_shapes <- 
  c("minus"="\u25BC",
    "plus"="\u25B2",
    "cat"="\u25FE")

#Lets make a summary graphic 

## add one dummy line for pinopsin
pinopsin_line <- 
  as.data.frame(
    cbind("mean_omega_pinopsin", 0.1, 0.1, 0.1, 1, "ED_SL", 1, 1, "N.S", "N.S", "N.S", "minus"))
colnames(pinopsin_line) <- colnames(astral_ecology_omega_corr_teleost)
astral_ecology_omega_corr_teleost <- rbind(astral_ecology_omega_corr_teleost, pinopsin_line)
astral_ecology_omega_corr_teleost$R2 <- as.numeric(astral_ecology_omega_corr_teleost$R2)

##Fill missing corr because of not enough data...

uniq_resp <- astral_ecology_omega_corr_teleost$Response %>% unique()
uniq_pred <- astral_ecology_omega_corr_teleost$predictor %>% unique()

for(curr_resp in uniq_resp){
  for(curr_pred in uniq_pred){
    
    nb_row <- 
      nrow(astral_ecology_omega_corr_teleost %>%
        filter(Response == curr_resp) %>%
        filter(predictor == curr_pred))
      
    if(nb_row < 1){
      
      new_line <- 
        as.data.frame(
          cbind(curr_resp, 0.1, 0.1, 0.1, 1, curr_pred, 1, 1, "N.S", "N.S", "N.S", "minus"))
      colnames(new_line) <- colnames(astral_ecology_omega_corr_teleost)
      astral_ecology_omega_corr_teleost <- rbind(astral_ecology_omega_corr_teleost, new_line)
      
    }
    
  }
}


astral_ecology_omega_corr_teleost$R2 <- as.numeric(astral_ecology_omega_corr_teleost$R2)






astral_ecology_omega_corr_teleost$Response <-
  factor(astral_ecology_omega_corr_teleost$Response ,
         levels=c("mean_omega_rh1", "mean_omega_exorh", "mean_omega_rh2", "mean_omega_sws2", 
                  "mean_omega_sws1", "mean_omega_lws", "mean_omega_pinopsin", "mean_omega_va", 
                  "mean_omega_parapinopsin","mean_omega_parietopsin", "mean_omega_tmt1", "mean_omega_tmt2", 
                  "mean_omega_tmt3","mean_omega_opn3", "mean_omega_opn8b", "mean_omega_opn8a", 
                  "mean_omega_opn8c", "mean_omega_opn7b", "mean_omega_opn7a", "mean_omega_opn5", 
                  "mean_omega_opn9", "mean_omega_opn6", "mean_omega_rgr", "mean_omega_rrh", "mean_omega_opn4m1_3", 
                  "mean_omega_opn4m2", "mean_omega_opn4x", "mean_omega_non_visual", "mean_omega_visual", "mean_omega"))


astral_ecology_omega_corr_teleost$predictor <-
  factor(astral_ecology_omega_corr_teleost$predictor ,
         levels=c("Absolute_SL", "Absolute_ED", 
                  "ED_SL", "Electrogenic_simplified", "Food_cat","AnaCat", 
                  "Diel_pattern_simplified",
                  "EnvTemp", "Fresh_salt",
                  "depth_category"
         ))



Signif_colors <-
  c("N.S" = 0,
    "Significant" = 1)




ggplot(astral_ecology_omega_corr_teleost, 
       aes(y=predictor, x=Response, fill= R2, alpha=Bonferroni_significant, color="black")) + 
  geom_tile(color = "black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis_c(option = "plasma") +
  scale_alpha_manual(values = Signif_colors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("")

ggplot(astral_ecology_omega_corr_teleost, 
       aes(y=predictor, x=Response, fill= R2, alpha=Bonferroni_significant, color="black")) + 
  geom_tile(color = "black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_minimal() +
  scale_fill_viridis_c(option = "plasma") +
  scale_alpha_manual(values = Signif_colors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("")





##### Omega vs Ecology  -- all graphics  -------------------



significant_correlations <- 
  astral_ecology_omega_corr %>%
  filter(Bonferroni_significant == "Significant")



plot_list <- list()
all_plot_names <- c()
for (row in 1:nrow(significant_correlations)){
  
  curr_lambda <- significant_correlations[row, "lambda"]
  curr_R2  <- significant_correlations[row, "R2"]
  curr_response  <- significant_correlations[row, "Response"]
  curr_predictor  <- significant_correlations[row, "predictor"]
  
  if(curr_predictor == "deep_depth"){
    curr_predictor <- "depth_deep" 
  }
  if(curr_predictor == "mean_depth"){
    curr_predictor <- "depth_mean" 
  }
  if(curr_predictor == "shallow_depth"){
    curr_predictor <- "depth_shallow" 
  }
  
  
  curr_pvalue <- significant_correlations[row, "FDR_pvalue"]
  slope <- significant_correlations[row, "slope"]
  
  
  curr_all_ospins_complete_wide_omega <- 
    all_ospins_complete_wide_omega %>%
    dplyr::select(species, curr_predictor, curr_response)
  colnames(curr_all_ospins_complete_wide_omega) <- c("species", "predictor", "response")
  
  
  curr_list_name <- paste(curr_response, curr_predictor, sep="_")
  
  if(is.na(slope)){
    
    good_preds <- 
      curr_all_ospins_complete_wide_omega %>%
      group_by(predictor) %>%
      summarise(count = n()) %>%
      filter(count > 2) %>%
      pull(predictor)
    
    curr_all_ospins_complete_wide_omega <-
      curr_all_ospins_complete_wide_omega %>%
      filter(predictor %in% good_preds)
    
    
    p = curr_all_ospins_complete_wide_omega %>%
      filter(! is.na(predictor)) %>%
      ggplot(., aes(x=predictor, y=response)) +
      geom_violin() + 
      stat_summary(fun=mean, geom="point", shape=18, size=4) +
      theme_classic() +
      labs(subtitle = bquote(FDR.pvalue ~ "=" ~.(curr_pvalue) ~ ";" ~ lambda ~ "=" ~ .(curr_lambda))) +
      xlab(curr_predictor) +
      ylab(curr_response) +
      theme(axis.text=element_text(size=16),
            axis.title=element_text(size=18),
            plot.subtitle=element_text(size=16),
            legend.position="none") 
    
  } else {
    
    
    curr_lm <- 
      lm(data = curr_all_ospins_complete_wide_omega %>% filter(! is.na(predictor)),
         formula = response ~ as.numeric(predictor))
    curr_function <- GLS_function(curr_lm)
    
    p = curr_all_ospins_complete_wide_omega %>%
      filter(! is.na(predictor)) %>%
      ggplot(., aes(x=predictor, y=response)) +
      geom_point() +
      theme_classic() +
      stat_function(fun = curr_function, color="black") +
      labs(subtitle = bquote(R^2 ~ "=" ~ .(curr_R2) ~ ";" ~ FDR.pvalue ~ "=" ~.(curr_pvalue) ~ ";" ~ lambda ~ "=" ~ .(curr_lambda))) +
      xlab(curr_predictor) +
      ylab(curr_response) +
      theme(axis.text=element_text(size=16),
            axis.title=element_text(size=18),
            plot.subtitle=element_text(size=16),
            legend.position="none") 
    
    
  }
  
  
  
  plot_list[[curr_list_name]] = p
  all_plot_names <- c(all_plot_names, curr_list_name)
}


for (curr_plot in all_plot_names) {
  print(plot_list[[curr_plot]])
}




##### Omega vs Ecology  -- all graphics -- Teleost NO WGD -------------------



significant_correlations <- 
  astral_ecology_omega_corr_teleost %>%
  filter(Bonferroni_significant == "Significant")



plot_list <- list()
all_plot_names <- c()
for (row in 1:nrow(significant_correlations)){
  
  curr_lambda <- significant_correlations[row, "lambda"]
  curr_R2  <- significant_correlations[row, "R2"]
  curr_response  <- significant_correlations[row, "Response"]
  curr_predictor  <- significant_correlations[row, "predictor"]
  
  if(curr_predictor == "deep_depth"){
    curr_predictor <- "depth_deep" 
  }
  if(curr_predictor == "mean_depth"){
    curr_predictor <- "depth_mean" 
  }
  if(curr_predictor == "shallow_depth"){
    curr_predictor <- "depth_shallow" 
  }
  
  
  curr_pvalue <- significant_correlations[row, "FDR_pvalue"]
  slope <- significant_correlations[row, "slope"]
  
  
  curr_all_ospins_complete_wide_omega <- 
    all_ospins_complete_wide_omega %>%
    filter(infraclass == "Teleostei") %>%
    filter(WGD == "No") %>%
    dplyr::select(species, curr_predictor, curr_response)
  colnames(curr_all_ospins_complete_wide_omega) <- c("species", "predictor", "response")
  
  
  curr_list_name <- paste(curr_response, curr_predictor, sep="_")
  
  if(is.na(slope)){
    
    good_preds <- 
      curr_all_ospins_complete_wide_omega %>%
      group_by(predictor) %>%
      summarise(count = n()) %>%
      filter(count > 2) %>%
      pull(predictor)
    
    curr_all_ospins_complete_wide_omega <-
      curr_all_ospins_complete_wide_omega %>%
      filter(predictor %in% good_preds)
    
    
    p = curr_all_ospins_complete_wide_omega %>%
      filter(! is.na(predictor)) %>%
      ggplot(., aes(x=predictor, y=response)) +
      geom_violin() + 
      stat_summary(fun=mean, geom="point", shape=18, size=4) +
      theme_classic() +
      #labs(subtitle = bquote(FDR.pvalue ~ "=" ~.(curr_pvalue) ~ ";" ~ lambda ~ "=" ~ .(curr_lambda))) +
      xlab(curr_predictor) +
      ylab(curr_response) +
      theme(axis.text=element_text(size=16),
            axis.title=element_text(size=18),
            plot.subtitle=element_text(size=16),
            legend.position="none") 
    
  } else {
    
    
    curr_lm <- 
      lm(data = curr_all_ospins_complete_wide_omega %>% filter(! is.na(predictor)),
         formula = response ~ as.numeric(predictor))
    curr_function <- GLS_function(curr_lm)
    
    p = curr_all_ospins_complete_wide_omega %>%
      filter(! is.na(predictor)) %>%
      ggplot(., aes(x=predictor, y=response)) +
      geom_point() +
      theme_classic() +
      stat_function(fun = curr_function, color="black") +
      #labs(subtitle = bquote(R^2 ~ "=" ~ .(curr_R2) ~ ";" ~ FDR.pvalue ~ "=" ~.(curr_pvalue) ~ ";" ~ lambda ~ "=" ~ .(curr_lambda))) +
      xlab(curr_predictor) +
      ylab(curr_response) +
      theme(axis.text=element_text(size=16),
            axis.title=element_text(size=18),
            plot.subtitle=element_text(size=16),
            legend.position="none") 
    
    
  }
  
  
  
  plot_list[[curr_list_name]] = p
  all_plot_names <- c(all_plot_names, curr_list_name)
}


for (curr_plot in all_plot_names) {
  print(plot_list[[curr_plot]])
}




 
all_ospins_complete_wide_omega %>%
  filter(infraclass == "Teleostei") %>%
  filter(WGD == "No") %>%
  dplyr::select(species, mean_omega_rh2, Diel_pattern_simplified) %>%
  filter(! is.na(Diel_pattern_simplified)) %>%
  group_by(Diel_pattern_simplified) %>%
  summarise(count=n(),
            mymean = mean(mean_omega_rh2, na.rm = TRUE))




all_ospins_complete_wide_omega %>%
  filter(infraclass == "Teleostei") %>%
  filter(WGD == "No") %>%
  dplyr::select(species, EnvTemp, mean_omega_opn9) %>%
  filter(! is.na(EnvTemp)) %>%
  group_by(EnvTemp) %>%
  summarise(count=n(),
            mymean = mean(mean_omega_opn9, na.rm = TRUE)) %>%
  arrange(mymean)


all_ospins_complete_wide_omega %>%
  filter(infraclass == "Teleostei") %>%
  filter(WGD == "No") %>%
  dplyr::select(species, depth_category, mean_omega_lws) %>%
  filter(! is.na(depth_category)) %>%
  group_by(depth_category) %>%
  summarise(count=n(),
            mymean = mean(mean_omega_lws, na.rm = TRUE)) %>%
  arrange(mymean)

all_ospins_complete_wide_omega %>%
  filter(infraclass == "Teleostei") %>%
  filter(WGD == "No") %>%
  dplyr::select(species, depth_category, mean_omega_exorh) %>%
  filter(! is.na(depth_category)) %>%
  group_by(depth_category) %>%
  summarise(count=n(),
            mymean = mean(mean_omega_exorh, na.rm = TRUE)) %>%
  arrange(mymean)


all_ospins_complete_wide_omega %>%
  filter(infraclass == "Teleostei") %>%
  filter(WGD == "No") %>%
  dplyr::select(species, depth_category, mean_omega_opn8b) %>%
  filter(! is.na(depth_category)) %>%
  group_by(depth_category) %>%
  summarise(count=n(),
            mymean = mean(mean_omega_opn8b, na.rm = TRUE)) %>%
  arrange(mymean)


all_ospins_complete_wide_omega %>%
  filter(infraclass == "Teleostei") %>%
  filter(WGD == "No") %>%
  dplyr::select(species, depth_category, mean_omega_rh2) %>%
  filter(! is.na(depth_category)) %>%
  group_by(depth_category) %>%
  summarise(count=n(),
            mymean = mean(mean_omega_rh2, na.rm = TRUE)) %>%
  arrange(mymean)


all_ospins_complete_wide_omega %>%
  filter(infraclass == "Teleostei") %>%
  filter(WGD == "No") %>%
  dplyr::select(species, Diel_pattern_simplified, mean_omega_rgr) %>%
  filter(! is.na(Diel_pattern_simplified)) %>%
  group_by(Diel_pattern_simplified) %>%
  summarise(count=n(),
            mymean = mean(mean_omega_rgr, na.rm = TRUE)) %>%
  arrange(mymean)


all_ospins_complete_wide_omega %>%
  filter(infraclass == "Teleostei") %>%
  filter(WGD == "No") %>%
  dplyr::select(species, AnaCat, mean_omega_rh2) %>%
  filter(! is.na(AnaCat)) %>%
  group_by(AnaCat) %>%
  summarise(count=n(),
            mymean = mean(mean_omega_rh2, na.rm = TRUE)) %>%
  arrange(mymean)

all_ospins_complete_wide_omega %>%
  filter(infraclass == "Teleostei") %>%
  filter(WGD == "No") %>%
  dplyr::select(species, Electrogenic_simplified, mean_omega_opn5) %>%
  filter(! is.na(Electrogenic_simplified)) %>%
  group_by(Electrogenic_simplified) %>%
  summarise(count=n(),
            mymean = mean(mean_omega_opn5, na.rm = TRUE)) %>%
  arrange(mymean)


all_ospins_complete_wide_omega %>%
  filter(infraclass == "Teleostei") %>%
  filter(WGD == "No") %>%
  dplyr::select(species, EnvTemp, mean_omega_visual) %>%
  filter(! is.na(EnvTemp)) %>%
  group_by(EnvTemp) %>%
  summarise(count=n(),
            mymean = mean(mean_omega_visual, na.rm = TRUE)) %>%
  arrange(mymean)

all_ospins_complete_wide_omega %>%
  filter(infraclass == "Teleostei") %>%
  filter(WGD == "No") %>%
  dplyr::select(species, EnvTemp, mean_omega) %>%
  filter(! is.na(EnvTemp)) %>%
  group_by(EnvTemp) %>%
  summarise(count=n(),
            mymean = mean(mean_omega, na.rm = TRUE)) %>%
  arrange(mymean)

#rh2 = Boreal >> Polar >> Deep-water >>>>>> tropical and subtrop the lowest
#sws1 =  Polar >> temperate >> Deep-water >>>>>> subtrop and tropical   (boreal) the lowest
#rh1 =  Polar  >> Deep-water >>>>>> subtrop and tropical   (boreal) the lowest
#opn8 and opn9 =  Polar  >> Deep-water >>>>>> subtrop and tropical   (boreal) the lowest


##### Omega vs Ecology  -- all graphics -- Teleost NO WGD - Better  -------------------


significant_correlations <- 
  astral_ecology_omega_corr_teleost %>%
  filter(Bonferroni_significant == "Significant")


### First depth and subfamilies

significant_correlation_depth <- 
  significant_correlations %>%
  filter(predictor == "depth_category") 


plot_list <- list()
all_plot_names <- c()
for (row in 1:nrow(significant_correlation_depth)){
  
  curr_lambda <- significant_correlation_depth[row, "lambda"]
  curr_R2  <- significant_correlation_depth[row, "R2"]
  curr_response  <- significant_correlation_depth[row, "Response"]
  curr_predictor  <- significant_correlation_depth[row, "predictor"]
  
  if(curr_predictor == "deep_depth"){
    curr_predictor <- "depth_deep" 
  }
  if(curr_predictor == "mean_depth"){
    curr_predictor <- "depth_mean" 
  }
  if(curr_predictor == "shallow_depth"){
    curr_predictor <- "depth_shallow" 
  }
  
  
  curr_pvalue <- significant_correlation_depth[row, "FDR_pvalue"]
  slope <- significant_correlation_depth[row, "slope"]
  
  curr_all_ospins_complete_wide_ecology <- 
    all_ospins_complete_wide_omega %>%
    filter(infraclass == "Teleostei") %>%
    filter(WGD == "No") %>%
    dplyr::select(species, curr_predictor, curr_response)
  colnames(curr_all_ospins_complete_wide_ecology) <- c("species", "predictor", "response")
  
  
  curr_list_name <- paste(curr_response, curr_predictor, sep="_")
  
  good_preds <- 
    curr_all_ospins_complete_wide_ecology %>%
    group_by(predictor) %>%
    summarise(count = n()) %>%
    filter(count > 2) %>%
    pull(predictor)
  
  curr_all_ospins_complete_wide_ecology <-
    curr_all_ospins_complete_wide_ecology %>%
    filter(predictor %in% good_preds)
  
  curr_all_ospins_complete_wide_ecology$predictor <-
    factor(curr_all_ospins_complete_wide_ecology$predictor ,
           levels=c("shallow", "mesophotic", "rariphotic", 
                    "mesopelagic", "bathypelagic"
           ))
  

  p = curr_all_ospins_complete_wide_ecology %>%
    filter(! is.na(predictor)) %>%
    ggplot(., aes(x=predictor, y=response)) +
    geom_violin() + 
    geom_jitter(shape=16, color="gray", alpha=0.8, position=position_jitter(0.2)) +
    stat_summary(fun=mean, geom="point", shape=18, size=2, color="red") +
    #stat_summary(fun=mean, geom="line", color="red", size=4) +
    theme_classic() +
    #labs(subtitle = bquote(FDR.pvalue ~ "=" ~.(curr_pvalue) ~ ";" ~ lambda ~ "=" ~ .(curr_lambda))) +
    xlab("") +
    ylab(curr_response) +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          plot.subtitle=element_text(size=16),
          legend.position="none")
  
  
  plot_list[[curr_list_name]] = p
  all_plot_names <- c(all_plot_names, curr_list_name)
}


for (curr_plot in all_plot_names) {
  print(plot_list[[curr_plot]])
}



### electrogenic 

significant_correlation_electro <- 
  significant_correlations %>%
  filter(predictor == "Electrogenic_simplified") 


plot_list <- list()
all_plot_names <- c()
for (row in 1:nrow(significant_correlation_electro)){
  
  curr_lambda <- significant_correlation_electro[row, "lambda"]
  curr_R2  <- significant_correlation_electro[row, "R2"]
  curr_response  <- significant_correlation_electro[row, "Response"]
  curr_predictor  <- significant_correlation_electro[row, "predictor"]
  
  
  curr_pvalue <- significant_correlation_electro[row, "FDR_pvalue"]
  slope <- significant_correlation_electro[row, "slope"]
  
  curr_all_ospins_complete_wide_ecology <- 
    all_ospins_complete_wide_omega %>%
    filter(infraclass == "Teleostei") %>%
    filter(WGD == "No") %>%
    dplyr::select(species, curr_predictor, curr_response)
  colnames(curr_all_ospins_complete_wide_ecology) <- c("species", "predictor", "response")
  
  
  curr_list_name <- paste(curr_response, curr_predictor, sep="_")
  
  good_preds <- 
    curr_all_ospins_complete_wide_ecology %>%
    group_by(predictor) %>%
    summarise(count = n()) %>%
    filter(count > 2) %>%
    pull(predictor)
  
  curr_all_ospins_complete_wide_ecology <-
    curr_all_ospins_complete_wide_ecology %>%
    filter(predictor %in% good_preds)
  
  curr_all_ospins_complete_wide_ecology$predictor <-
    factor(curr_all_ospins_complete_wide_ecology$predictor ,
           levels=c("no special ability", "discharging"
           ))
  
  p = curr_all_ospins_complete_wide_ecology %>%
    filter(! is.na(predictor)) %>%
    ggplot(., aes(x=predictor, y=response)) +
    geom_violin() + 
    geom_jitter(shape=16, color="gray", alpha=0.8, position=position_jitter(0.2)) +
    stat_summary(fun=mean, geom="point", shape=18, size=2, color="red") +
    #stat_summary(fun=mean, geom="line", color="red", size=4) +
    theme_classic() +
    #labs(subtitle = bquote(FDR.pvalue ~ "=" ~.(curr_pvalue) ~ ";" ~ lambda ~ "=" ~ .(curr_lambda))) +
    xlab("") +
    ylab(curr_response) +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          plot.subtitle=element_text(size=16),
          legend.position="none")
  
  
  plot_list[[curr_list_name]] = p
  all_plot_names <- c(all_plot_names, curr_list_name)
}


for (curr_plot in all_plot_names) {
  print(plot_list[[curr_plot]])
}






### Climate zone 

significant_correlation_Env <- 
  significant_correlations %>%
  filter(predictor == "EnvTemp") 


plot_list <- list()
all_plot_names <- c()
for (row in 1:nrow(significant_correlation_Env)){
  
  curr_lambda <- significant_correlation_Env[row, "lambda"]
  curr_R2  <- significant_correlation_Env[row, "R2"]
  curr_response  <- significant_correlation_Env[row, "Response"]
  curr_predictor  <- significant_correlation_Env[row, "predictor"]
  
  
  curr_pvalue <- significant_correlation_Env[row, "FDR_pvalue"]
  slope <- significant_correlation_Env[row, "slope"]
  
  curr_all_ospins_complete_wide_ecology <- 
    all_ospins_complete_wide_omega %>%
    filter(infraclass == "Teleostei") %>%
    filter(WGD == "No") %>%
    dplyr::select(species, curr_predictor, curr_response)
  colnames(curr_all_ospins_complete_wide_ecology) <- c("species", "predictor", "response")
  
  
  curr_list_name <- paste(curr_response, curr_predictor, sep="_")
  
  good_preds <- 
    curr_all_ospins_complete_wide_ecology %>%
    group_by(predictor) %>%
    summarise(count = n()) %>%
    filter(count > 2) %>%
    pull(predictor)
  
  curr_all_ospins_complete_wide_ecology <-
    curr_all_ospins_complete_wide_ecology %>%
    filter(predictor %in% good_preds)
  
  curr_all_ospins_complete_wide_ecology$predictor <-
    factor(curr_all_ospins_complete_wide_ecology$predictor ,
           levels=c("tropical","subtropical", "temperate", "boreal", "polar", "deep-water"
           ))
  
  p = curr_all_ospins_complete_wide_ecology %>%
    filter(! is.na(predictor)) %>%
    ggplot(., aes(x=predictor, y=response)) +
    geom_violin() + 
    geom_jitter(shape=16, color="gray", alpha=0.8, position=position_jitter(0.2)) +
    stat_summary(fun=mean, geom="point", shape=18, size=2, color="red") +
    #stat_summary(fun=mean, geom="line", color="red", size=4) +
    theme_classic() +
    #labs(subtitle = bquote(FDR.pvalue ~ "=" ~.(curr_pvalue) ~ ";" ~ lambda ~ "=" ~ .(curr_lambda))) +
    xlab("") +
    ylab(curr_response) +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          plot.subtitle=element_text(size=16),
          legend.position="none")
  
  
  plot_list[[curr_list_name]] = p
  all_plot_names <- c(all_plot_names, curr_list_name)
}


for (curr_plot in all_plot_names) {
  print(plot_list[[curr_plot]])
}




### Diel

significant_correlation_Diel <- 
  significant_correlations %>%
  filter(predictor == "Diel_pattern_simplified") 


plot_list <- list()
all_plot_names <- c()
for (row in 1:nrow(significant_correlation_Diel)){
  
  curr_lambda <- significant_correlation_Diel[row, "lambda"]
  curr_R2  <- significant_correlation_Diel[row, "R2"]
  curr_response  <- significant_correlation_Diel[row, "Response"]
  curr_predictor  <- significant_correlation_Diel[row, "predictor"]
  
  
  curr_pvalue <- significant_correlation_Diel[row, "FDR_pvalue"]
  slope <- significant_correlation_Diel[row, "slope"]
  
  curr_all_ospins_complete_wide_ecology <- 
    all_ospins_complete_wide_omega %>%
    filter(infraclass == "Teleostei") %>%
    filter(WGD == "No") %>%
    dplyr::select(species, curr_predictor, curr_response)
  colnames(curr_all_ospins_complete_wide_ecology) <- c("species", "predictor", "response")
  
  
  curr_list_name <- paste(curr_response, curr_predictor, sep="_")
  
  good_preds <- 
    curr_all_ospins_complete_wide_ecology %>%
    group_by(predictor) %>%
    summarise(count = n()) %>%
    filter(count > 2) %>%
    pull(predictor)
  
  curr_all_ospins_complete_wide_ecology <-
    curr_all_ospins_complete_wide_ecology %>%
    filter(predictor %in% good_preds)
  
  curr_all_ospins_complete_wide_ecology$predictor <-
    factor(curr_all_ospins_complete_wide_ecology$predictor ,
           levels=c("Diurnal","Crepuscular", "Nocturnal", "Cathemeral"
           ))
  
  p = curr_all_ospins_complete_wide_ecology %>%
    filter(! is.na(predictor)) %>%
    ggplot(., aes(x=predictor, y=response)) +
    geom_violin() + 
    geom_jitter(shape=16, color="gray", alpha=0.8, position=position_jitter(0.2)) +
    stat_summary(fun=mean, geom="point", shape=18, size=2, color="red") +
    #stat_summary(fun=mean, geom="line", color="red", size=4) +
    theme_classic() +
    #labs(subtitle = bquote(FDR.pvalue ~ "=" ~.(curr_pvalue) ~ ";" ~ lambda ~ "=" ~ .(curr_lambda))) +
    xlab("") +
    ylab(curr_response) +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          plot.subtitle=element_text(size=16),
          legend.position="none")
  
  
  plot_list[[curr_list_name]] = p
  all_plot_names <- c(all_plot_names, curr_list_name)
}


for (curr_plot in all_plot_names) {
  print(plot_list[[curr_plot]])
}




### Migration

significant_correlation_AnaCat <- 
  significant_correlations %>%
  filter(predictor == "AnaCat") 


plot_list <- list()
all_plot_names <- c()
for (row in 1:nrow(significant_correlation_AnaCat)){
  
  curr_lambda <- significant_correlation_AnaCat[row, "lambda"]
  curr_R2  <- significant_correlation_AnaCat[row, "R2"]
  curr_response  <- significant_correlation_AnaCat[row, "Response"]
  curr_predictor  <- significant_correlation_AnaCat[row, "predictor"]
  
  
  curr_pvalue <- significant_correlation_AnaCat[row, "FDR_pvalue"]
  slope <- significant_correlation_AnaCat[row, "slope"]
  
  curr_all_ospins_complete_wide_ecology <- 
    all_ospins_complete_wide_omega %>%
    filter(infraclass == "Teleostei") %>%
    filter(WGD == "No") %>%
    dplyr::select(species, curr_predictor, curr_response)
  colnames(curr_all_ospins_complete_wide_ecology) <- c("species", "predictor", "response")
  
  
  curr_list_name <- paste(curr_response, curr_predictor, sep="_")
  
  good_preds <- 
    curr_all_ospins_complete_wide_ecology %>%
    group_by(predictor) %>%
    summarise(count = n()) %>%
    filter(count > 2) %>%
    pull(predictor)
  
  curr_all_ospins_complete_wide_ecology <-
    curr_all_ospins_complete_wide_ecology %>%
    filter(predictor %in% good_preds)
  
  curr_all_ospins_complete_wide_ecology$predictor <-
    factor(curr_all_ospins_complete_wide_ecology$predictor ,
           levels=c("non-migratory","potamodromous", "catadromous", "oceanodromous",
                    "anadromous", "amphidromous"
           ))
  
  p = curr_all_ospins_complete_wide_ecology %>%
    filter(! is.na(predictor)) %>%
    ggplot(., aes(x=predictor, y=response)) +
    geom_violin() + 
    geom_jitter(shape=16, color="gray", alpha=0.8, position=position_jitter(0.2)) +
    stat_summary(fun=mean, geom="point", shape=18, size=2, color="red") +
    #stat_summary(fun=mean, geom="line", color="red", size=4) +
    theme_classic() +
    #labs(subtitle = bquote(FDR.pvalue ~ "=" ~.(curr_pvalue) ~ ";" ~ lambda ~ "=" ~ .(curr_lambda))) +
    xlab("") +
    ylab(curr_response) +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          plot.subtitle=element_text(size=16),
          legend.position="none")
  
  
  plot_list[[curr_list_name]] = p
  all_plot_names <- c(all_plot_names, curr_list_name)
}


for (curr_plot in all_plot_names) {
  print(plot_list[[curr_plot]])
}


