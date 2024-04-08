#### Libraries  ---------------------------------

rm(list = ls())

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
library(ggtree)
library("castor")
library(plotly)
library(ggfortify)
library("erer")
library(ggtreeExtra)
library("ggstar")
library(MASS)
library(PCAtest)
library(treemapify)
library(umap)
library(gt)
#library("MoreTreeTools")

split_tibble <- function(tibble, column = 'col') {
  tibble %>% split(., .[,column]) %>% lapply(., function(x) x[,setdiff(names(x),column)])
}

GLS_function <- function(gls_rslt) {
  GLS_cc <- coef(gls_rslt)
  GLS_f <- function(x) GLS_cc[1] + GLS_cc[2]*x
  return(GLS_f)
}

confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}

##### Data load - Colors -----------------


library(RColorBrewer)
pal <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
         "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
         "#920000","#924900","#db6d00","#24ff24","#ffff6d")

pie(rep(1,15), col=pal)


NVO_opsins_clades_colors <- 
  c("exorh"="#E6F5D0",
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
    "rrh"="#ABD9E9",
    "tmt1"="#FC4E2A",
    "tmt2"="#045A8D",
    "tmt3"="#D1E5F0",
    "va"="#8C510A"
  )

tissues_colors <- 
  c("brain"="#000000",
    "bones"="#004949",
    "eye"="#ffff6d",
    "fin"="#009292",
    "gill"="#ff6db6",
    "gut"="#ffb6db",
    "heart"="#490092",
    "kidney"="#006ddb",
    "liver"="#b66dff",
    "muscle"="#6db6ff",
    "ovary"="#b6dbff",
    "skin"="#920000",
    "spleen"="#924900",
    "testis"="#db6d00"
  )

VO_NVO_colors <- 
  c("NVO"="#A48528",
    "VO"="#0472D2"
  )


Cone_vs_Rho_colors <- 
  c("cone"="gray",
    "rho"="black"
  )


Cone_colors <- 
  c("rh2"="#229C88",
    "lws"="#FBB5B5",
    "sws2"="#7EC3FF",
    "sws1"="#9F1BD8")


##### Data load - Species info table -----------------


species_table <- 
  read.table("Actino_species_table_info_v2.tsv",
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


##### Data load - Species tree -----------------

#load phylogenetic tree

ASTRAL_species_tree_nodelabel <- 
  read.tree("AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel")
ASTRAL_species_tree_nodelabel$edge.length <- ASTRAL_species_tree_nodelabel$edge.length * 1000


##### Data load - Opsins table   ---------------------------------

opsins_table <- 
  read.table("ALL_opsins_classifications.csv",
             sep=",",
             header=FALSE)

colnames(opsins_table) <- c("gene_state", "gene_name", "clade", "species")




#Count the total nb of genes per species and per clade

all_opsins_table_count <- 
  as.data.frame(
    opsins_table %>%
      group_by(species, gene_state, clade) %>%
      summarise(count = n())
  )




###



visual_opsins_table_count <- 
  as.data.frame(
    opsins_table %>%
      filter(clade %in% c("lws", "rh1", "rh2", "sws1", "sws2")) %>%
      group_by(species, gene_state, clade) %>%
      summarise(count = n())
  )

non_visual_opsins_table_count <- 
  as.data.frame(
    opsins_table %>%
      filter(! clade %in% c("lws", "rh1", "rh2", "sws1", "sws2")) %>%
      group_by(species, gene_state, clade) %>%
      summarise(count = n())
  )








#Transform visual opsin table from long to wide format

visual_opsins_table_count_total <- 
  as.data.frame(
    visual_opsins_table_count %>%
      group_by(species, clade) %>%
      summarise(total = sum(count)))

visual_opsins_table_count_wide <- 
  as.data.frame(
    visual_opsins_table_count_total %>%
      pivot_wider(names_from = clade, values_from = total)
  )
visual_opsins_table_count_wide[is.na(visual_opsins_table_count_wide)] <- 0



#Transform non visual opsin table from long to wide format

non_visual_opsins_table_count_total <- 
  as.data.frame(
    non_visual_opsins_table_count %>%
      group_by(species, clade) %>%
      summarise(total = sum(count)))

non_visual_opsins_table_count_wide <- 
  as.data.frame(
    non_visual_opsins_table_count_total %>%
      pivot_wider(names_from = clade, values_from = total)
  )
non_visual_opsins_table_count_wide[is.na(non_visual_opsins_table_count_wide)] <- 0




# Make table with only Complete genes 
visual_opsins_table_count_complete <- 
  visual_opsins_table_count %>%
  filter(gene_state == "Complete")

non_visual_opsins_table_count_complete <- 
  non_visual_opsins_table_count %>%
  filter(gene_state == "Complete")


all_opsins_table_count_complete <- 
  rbind(visual_opsins_table_count_complete,
        non_visual_opsins_table_count_complete)

#Transform visual opsin table from long to wide format (only complete)
visual_opsins_table_count_complete_wide <- 
  as.data.frame(
    visual_opsins_table_count_complete %>%
      dplyr::select(! gene_state) %>%
      pivot_wider(names_from = clade, values_from = count)
  )

visual_opsins_table_count_complete_wide <- 
  replace(visual_opsins_table_count_complete_wide, is.na(visual_opsins_table_count_complete_wide), 0)

Prietella_table <- as.data.frame(cbind("Prietella_phreatophila", 0, 0, 0, 0, 0)) 
colnames(Prietella_table) <- colnames(visual_opsins_table_count_complete_wide)
visual_opsins_table_count_complete_wide <- 
  rbind(visual_opsins_table_count_complete_wide,
        Prietella_table)

visual_opsins_table_count_complete_wide$lws <- as.numeric(visual_opsins_table_count_complete_wide$lws)
visual_opsins_table_count_complete_wide$rh1 <- as.numeric(visual_opsins_table_count_complete_wide$rh1)
visual_opsins_table_count_complete_wide$rh2 <- as.numeric(visual_opsins_table_count_complete_wide$rh2)
visual_opsins_table_count_complete_wide$sws1 <- as.numeric(visual_opsins_table_count_complete_wide$sws1)
visual_opsins_table_count_complete_wide$sws2 <- as.numeric(visual_opsins_table_count_complete_wide$sws2)



#Transform non visual opsin table from long to wide format (only complete)

non_visual_opsins_table_count_complete_wide <- 
  as.data.frame(
    non_visual_opsins_table_count_complete %>%
      dplyr::select(! gene_state) %>%
      pivot_wider(names_from = clade, values_from = count)
  )

non_visual_opsins_table_count_complete_wide <- 
  replace(non_visual_opsins_table_count_complete_wide, is.na(non_visual_opsins_table_count_complete_wide), 0)


visual_opsins_table_count_complete_wide <- 
  visual_opsins_table_count_complete_wide %>%
  mutate(total_visual = rowSums(across(where(is.numeric))))


non_visual_opsins_table_count_complete_wide <- 
  non_visual_opsins_table_count_complete_wide %>%
  mutate(total_non_visual = rowSums(across(where(is.numeric))))




#Combine COMPLETE visual and non visual opsins tables
all_ospins_complete_wide <- 
  left_join(non_visual_opsins_table_count_complete_wide, visual_opsins_table_count_complete_wide, by="species")

all_ospins_complete_wide <- 
  replace(all_ospins_complete_wide, is.na(all_ospins_complete_wide), 0)


#Add the total of cone opsins genes


all_ospins_complete_wide <- 
  all_ospins_complete_wide %>%
  mutate(total_cone_visual = 
           sws1 + sws2 + lws + rh2)


#Add the total of opsins

all_ospins_complete_wide <- 
  all_ospins_complete_wide %>%
  mutate(total_opsins = 
           total_visual + total_non_visual)


#long to wide to long again to have all the info in the long tables

all_opsins_table_count_complete_wide <- 
  as.data.frame(
    all_opsins_table_count_complete %>%
      dplyr::select(! gene_state) %>%
      pivot_wider(names_from = clade, values_from = count)
  )
all_opsins_table_count_complete_wide[is.na(all_opsins_table_count_complete_wide)] <- 0

all_opsins_table_count_complete <- 
  as.data.frame(
    all_opsins_table_count_complete_wide %>%
      pivot_longer(!species, names_to = "clade", values_to = "count")
  )





##### Data load - OGG - gene   ---------------------------------

orthogroups_df <- 
  read.table("Table_corresp_OGG_gene.csv",
             sep=",",
             header=FALSE)
colnames(orthogroups_df) <- c("species", "OGG", "gene_name")

orthogroups_df <- 
  orthogroups_df %>%
  filter(! is.na(gene_name)) %>%
  filter(gene_name != "")


orthogroups_df %>%
  group_by(species) %>%
  summarise(n())



#### Import RNA-seq data  ---------------------------------

#Import read count on opsins
opsins_read_count_df <- 
  read.table("ALL_species_tissues_RNA.csv", sep=",")
colnames(opsins_read_count_df) <- c("species", "tissue", "gene_name", "exon_nb", "read_count")

#sum of read number per opsin

opsins_read_count_df <- 
  as.data.frame(
    opsins_read_count_df %>%
      group_by(species, tissue, gene_name) %>%
      summarise(sum_read_count = sum(read_count))
  )



#Import read count on all genes and all exons
all_genes_read_count_df <- 
  read.table("Annotated_species_tissues_all_genes_RNA.csv", sep=",")
colnames(all_genes_read_count_df) <- c("species", "tissue", "gene_name", "read_count")




#sum of read number per gene
all_genes_read_count_df <- 
  as.data.frame(
    all_genes_read_count_df %>%
      group_by(species, tissue, gene_name) %>%
      summarise(sum_read_count = sum(read_count))
  )


as.data.frame(
  all_genes_read_count_df %>%
  group_by(species, tissue) %>%
  summarise(count = n())) %>%
  dplyr::select(species, count) %>%
  distinct()





#### Filter RNA-seq data based on metrics  ---------------------------------


STAR_metrics_df <- 
  read.table("STAR_all_statistics.csv", sep=",", header=FALSE)
colnames(STAR_metrics_df) <- 
  c("species", "tissue", "stat", "value")


#First lets look at the percentage of mapped reads
perc_mapped_reads_df <- 
  STAR_metrics_df %>%
  filter(stat %in% c("Uniquely_mapped_reads_perc", "perc__of_reads_mapped_to_multiple_loci")) %>%
  group_by(species, tissue) %>%
  summarise(
    percentage_reads_mapped = 
      sum(value))



perc_mapped_reads_df %>%
  ungroup() %>%
  ggplot(., aes(x=species, y=percentage_reads_mapped)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.text.x=element_blank()) +
  theme(legend.position="none") +
  ylab("% of mapped reads") +
  xlab("species") +
  geom_hline(yintercept = 70, color="red", linetype="dashed")
  



bad_sample_perc_df <- 
  perc_mapped_reads_df %>%
  filter(percentage_reads_mapped < 70)
bad_sample_perc_df$sample <- 
  paste(bad_sample_perc_df$species,
        bad_sample_perc_df$tissue,
        sep="_")
bad_sample_perc <- bad_sample_perc_df$sample

#Now lets look at the number of reads


number_mapped_reads_df <- 
  STAR_metrics_df %>%
  filter(stat %in% c("Uniquely_mapped_reads_number", "Number_of_reads_mapped_to_multiple_loci")) %>%
  group_by(species, tissue) %>%
  summarise(
    number_reads_mapped = 
      sum(value))

number_mapped_reads_df$sample <- 
  paste(number_mapped_reads_df$species,
        number_mapped_reads_df$tissue,
        sep="_")



number_mapped_reads_df %>%
  ungroup() %>%
  mutate(! sample %in% bad_sample_perc) %>%
  ggplot(., aes(x=species, y=number_reads_mapped)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.text.x=element_blank()) +
  theme(legend.position="none") +
  ylab("Number of mapped reads") +
  xlab("species") +
  geom_hline(yintercept = 5000000, color="red", linetype="dashed")


bad_sample_number <-
  number_mapped_reads_df %>%
  filter(number_reads_mapped < 5000000) %>%
  pull(sample)


all_bad_sample <- c(bad_sample_number, bad_sample_perc) %>% unique()


#### Import genes lengths  ---------------------------------

#Import gene lengths for opsins

opsins_length_bp <- 
  read.table("TOTAL_opsins.fa.fai", sep="\t")
colnames(opsins_length_bp) <- c("gene_name", "bp", "offset", "linebases", "linewidth")
opsins_length_bp$gene_name <- gsub("---FALSE.*", "" ,opsins_length_bp$gene_name)
opsins_length_bp$gene_name <- gsub("---TRUE.*", "" ,opsins_length_bp$gene_name)
opsins_length_bp$gene_name <- gsub("-FALSE.*", "" ,opsins_length_bp$gene_name)
opsins_length_bp$gene_name <- gsub("-TRUE.*", "" ,opsins_length_bp$gene_name)



opsins_length_kbp <- 
  opsins_length_bp %>% 
  mutate(kbp = bp/1000) %>%
  dplyr::select(gene_name, bp, kbp)

#Import gene lengths for all annotated genes

allgenes_length_bp <- 
  read.table("Table_all_gene_length.csv", sep=",")
colnames(allgenes_length_bp) <- c("species", "gene_name", "aa_length")
allgenes_length_kbp <- 
  allgenes_length_bp %>% 
  mutate(bp = aa_length*3) %>%
  mutate(kbp = bp/1000) %>%
  dplyr::select(gene_name, bp, kbp)



#### Combine tables  ---------------------------------

#Add the gene subfamily and state to the opsins read counts
opsins_table_p <- opsins_table %>% dplyr::select(gene_state, gene_name, clade)
opsins_table_p$gene_name <- gsub("---FALSE.*", "" ,opsins_table_p$gene_name)
opsins_table_p$gene_name <- gsub("---TRUE.*", "" ,opsins_table_p$gene_name)
opsins_table_p$gene_name <- gsub("-FALSE.*", "" ,opsins_table_p$gene_name)
opsins_table_p$gene_name <- gsub("-TRUE.*", "" ,opsins_table_p$gene_name)


opsins_read_count_df_names <- 
  left_join(opsins_read_count_df, opsins_table_p, by="gene_name")

#Add the gene length to the opsins read count 
opsins_count_final_df <- 
  left_join(opsins_read_count_df_names, opsins_length_kbp, by="gene_name")

#Add the OGG to whole genes read counts
all_genes_read_count_df_OGG <- 
  left_join(orthogroups_df, all_genes_read_count_df, by=c("species", "gene_name"))


as.data.frame(
  all_genes_read_count_df_OGG %>%
  filter(! is.na(tissue)) %>%
  group_by(species, tissue) %>%
  summarise(count = n()) %>%
  dplyr::select(species, count) %>%
  distinct()
)
  
#Add the gene length to whole genes read counts

all_genes_read_count_df_OGG_length <- 
  left_join(all_genes_read_count_df_OGG, allgenes_length_kbp, by="gene_name")


as.data.frame(
  all_genes_read_count_df_OGG_length %>%
    filter(! is.na(tissue)) %>%
    filter(! is.na(bp)) %>%
    group_by(species, tissue) %>%
    summarise(count = n()) %>%
    dplyr::select(species, count) %>%
    distinct()
)

all_genes_count_final_df <- all_genes_read_count_df_OGG_length

rm(all_genes_read_count_df_OGG_length)
rm(all_genes_read_count_df_OGG)
rm(all_genes_read_count_df)


#### Compute RPK and TPM values -- Opsins ---------------------------------

## First compute RPK
opsins_count_final_df <- 
  opsins_count_final_df %>%
  mutate(RPK = sum_read_count/kbp)

opsins_count_final_onlyrpk <- opsins_count_final_df

## Now compute TPM

scaling_factor_opsins_df <- 
  as.data.frame(
    opsins_count_final_df %>%
  group_by(species, tissue) %>%
  summarise(sum_RPK = sum(RPK))
  )
scaling_factor_opsins_df <- 
  scaling_factor_opsins_df %>%
  mutate(scaling_factor = sum_RPK/1000000)
opsins_count_final_df <- 
  left_join(opsins_count_final_df, scaling_factor_opsins_df, by=c("species", "tissue"))
opsins_count_final_df <- 
  opsins_count_final_df %>%
  mutate(TPM = RPK/scaling_factor)
opsins_count_final_df[is.na(opsins_count_final_df$TPM), 'TPM'] <- 0



opsins_count_final_df %>%
  group_by(species, tissue) %>%
  summarise(sum_TPM = sum(TPM))

#### Compute RPK and TPM values -- ALL genes ---------------------------------

## First compute RPK
all_genes_count_final_df <- 
  all_genes_count_final_df %>%
  mutate(RPK = sum_read_count/kbp)

## Replace automatic opsins by manual opsins

Annot_opsins <- 
  read.table("Species_AnnotatedOpsins.csv",
             sep=",",
             header=FALSE)
colnames(Annot_opsins) <- c("species", "opsins_ID")
id_annotated_opsins <- Annot_opsins %>% pull(opsins_ID)

all_genes_count_final_df <- 
  all_genes_count_final_df %>%
  filter(! gene_name %in% id_annotated_opsins) 

good_sp <- all_genes_count_final_df %>% pull(species) %>% unique()

opsins_count_final_onlyrpk <- 
  opsins_count_final_onlyrpk %>%
  filter(species %in% good_sp) %>%
  dplyr::select(species, clade, gene_name, tissue, sum_read_count, bp, kbp, RPK) 
colnames(opsins_count_final_onlyrpk) <- colnames(all_genes_count_final_df)

all_genes_count_final_df <- 
  rbind(all_genes_count_final_df, opsins_count_final_onlyrpk)


opsins_count_final_onlyrpk %>%
  filter(species %in% good_sp)

## Now compute TPM

scaling_factor_allgenes_df <- 
  as.data.frame(
    all_genes_count_final_df %>%
      group_by(species, tissue) %>%
      summarise(sum_RPK = sum(RPK))
  ) %>% filter(! is.na(tissue))

scaling_factor_allgenes_df <- 
  scaling_factor_allgenes_df %>%
  mutate(scaling_factor = sum_RPK/1000000)
all_genes_count_final_df <- 
  left_join(all_genes_count_final_df, scaling_factor_allgenes_df, by=c("species", "tissue"))
all_genes_count_final_df <- 
  all_genes_count_final_df %>%
  mutate(TPM = RPK/scaling_factor)
all_genes_count_final_df[is.na(all_genes_count_final_df$TPM), 'TPM'] <- 0

all_genes_count_final_df <- all_genes_count_final_df %>% filter(! is.na(tissue))


all_genes_count_final_df %>%
  group_by(species, tissue) %>%
  summarise(sum_TPM = sum(TPM))


#### Remove bad samples -- Opsins + ALL genes ---------------------------------

#First create a sample column

all_genes_count_final_df$sample <- 
  paste(
    all_genes_count_final_df$species,
    all_genes_count_final_df$tissue,
    sep="_"
  )


opsins_count_final_df$sample <-
  paste(
    opsins_count_final_df$species,
    opsins_count_final_df$tissue,
    sep="_"
  )


#Now remove bad samples ...

all_genes_count_final_df <- 
  all_genes_count_final_df %>%
  filter(! sample %in% all_bad_sample)
all_genes_count_final_df <- all_genes_count_final_df %>% dplyr::select(-sample)

opsins_count_final_df <- 
  opsins_count_final_df %>%
  filter(! sample %in% all_bad_sample)
opsins_count_final_df <- opsins_count_final_df %>% dplyr::select(-sample)




#### Draw phylo tree of samples analysed ---------------------------------

RNA_sp <- 
  opsins_count_final_df %>%
  pull(species) %>% unique()
RNA_sp_tree <- 
  keep.tip(ASTRAL_species_tree_nodelabel, RNA_sp)

Annotated_sp <- all_genes_count_final_df %>% pull(species) %>% unique()
annotated_table <- 
  as.data.frame(RNA_sp) %>%
  mutate(annotation = if_else(
    RNA_sp %in% Annotated_sp,
    "Annotated",
    "Non_annotated"
  )) %>%
  mutate(species_space = gsub("_", " ", RNA_sp))

annotated_table <- annotated_table %>% dplyr::select("species_space", "RNA_sp", "annotation")
colnames(annotated_table) <- c("species_space","species", "annotation")

RNA_sp_tree_space <- RNA_sp_tree
RNA_sp_tree_space$tip.label <- 
  gsub("_", " ", RNA_sp_tree_space$tip.label)


Avail_df <- 
  opsins_count_final_df %>%
  dplyr::select(species, tissue) %>%
  distinct()

Avail_df$species_space <- gsub("_", " ", Avail_df$species)


RNA_sp_tree_space$tip.label[RNA_sp_tree_space$tip.label == "Gasterosteus aculeatus aculeatus"] <- "Gasterosteus aculeatus"
annotated_table$species_space[annotated_table$species_space == "Gasterosteus aculeatus aculeatus"] <- "Gasterosteus aculeatus"
Avail_df$species_space[Avail_df$species_space == "Gasterosteus aculeatus aculeatus"] <- "Gasterosteus aculeatus"


Avail_df$tissue <- 
  factor(Avail_df$tissue ,
         levels=c("eye", "brain", 
                  "fin", "gill", "skin","bones", 
                  "muscle","heart", "spleen",
                  "gut", "liver", "kidney", "testis",
                  "ovary"
         ))



ggtree(RNA_sp_tree_space, layout="circular", size=0.8) %<+% annotated_table  +
  geom_tiplab(size=2, offset=0.1, fontface=3) +
  aes(color=annotation) +
  scale_color_manual(values = c("Annotated" = "red", "Non_annotated" = "black")) +
  ggnewscale::new_scale_colour() + 
  geom_fruit(
    data=Avail_df,
    geom=geom_star,
    mapping=aes(x=tissue, y=species_space, fill=tissue, starshape=tissue),
    size=1,
    pwidth=2,
    offset=0.8,
    inherit.aes = FALSE,
    grid.params=list(
      linetype=3,
      size=0.2
    )
    
  ) +
  scale_starshape_manual(
    values=c(rep(15, 14)),
    guide="none"
    ) + 
  scale_fill_manual(
    values=tissues_colors
    ) +
  theme(legend.position = "none") 






#### All genes to All OGG  ---------------------------------

OGG_count_final_df <- 
  as.data.frame(
    all_genes_count_final_df %>%
      group_by(species, OGG, tissue) %>%
      summarise(OGG_TPM = sum(TPM))
  )

OGG_count_final_df$sample <- 
  paste(OGG_count_final_df$species, OGG_count_final_df$tissue, sep = "_")


#### All opsins to All Subfamilies  ---------------------------------

opsin_subfam_count_final_df <- 
  as.data.frame(
    opsins_count_final_df %>%
      group_by(species, clade, tissue) %>%
      summarise(clade_TPM = sum(TPM))
  )

opsin_subfam_count_final_df$sample <- 
  paste(opsin_subfam_count_final_df$species, opsin_subfam_count_final_df$tissue, sep = "_")



#### PCA using shared OGG  ---------------------------------


# pivot the dataframe

OGG_count_final_df_wide <- 
  as.data.frame(
    pivot_wider(OGG_count_final_df, 
                names_from = c(OGG), values_from = OGG_TPM)
  )

OGG_count_final_df_wide <- OGG_count_final_df_wide %>% 
  mutate_if(is.numeric, ~replace_na(., 0))


#Remove info columns not used by the PCA command

OGG_count_final_df_wide_woinfo <- 
  OGG_count_final_df_wide %>%
  dplyr::select(- c(species, tissue))

rownames(OGG_count_final_df_wide_woinfo) <- OGG_count_final_df_wide_woinfo[,1]
OGG_count_final_df_wide_woinfo <- OGG_count_final_df_wide_woinfo[,-1]


#remove 0 vaiance OGG (usually all set to 0)
OGG_count_final_df_wide_woinfo <- 
  OGG_count_final_df_wide_woinfo[ , which(apply(OGG_count_final_df_wide_woinfo, 2, var) != 0)]


#Lets extract only shared OGG

shared_OGG <- 
  orthogroups_df %>%
  dplyr::select(species, OGG) %>%
  distinct() %>%
  group_by(OGG) %>%
  summarise(nb_species = n()) %>%
  filter(nb_species == 45) %>%
  pull(OGG)
shared_OGG <- intersect(shared_OGG, OGG_count_final_df %>% pull(OGG) %>% unique())
shared_OGG <- c(shared_OGG, opsins_count_final_df %>% pull(clade) %>% unique())


OGG_count_final_df_wide_shared <- 
  OGG_count_final_df_wide %>%
  dplyr::select(species, tissue, sample, shared_OGG)

OGG_count_final_df_wide_woinfo_shared <- 
  OGG_count_final_df_wide_woinfo %>%
  dplyr::select(shared_OGG)


#Lets perform the PCA !

#PCA_shared_genes <- prcomp(OGG_count_final_df_wide_woinfo_shared, scale = TRUE)
#autoplot(PCA_shared_genes, data = OGG_count_final_df_wide_shared, colour = 'tissue') +
#  theme_classic()

PCA_shared_genes_log <- prcomp(log(OGG_count_final_df_wide_woinfo_shared + 1), scale = TRUE)

autoplot(PCA_shared_genes_log, data = OGG_count_final_df_wide_shared, colour = 'tissue') +
  theme_classic()


PC1 <- PCA_shared_genes_log$x[,1]
PC2 <- PCA_shared_genes_log$x[,2]


OGG_count_final_df_wide_shared$tissue <- 
  factor(OGG_count_final_df_wide_shared$tissue ,
         levels=c("eye", "brain", 
                  "fin", "gill", "skin","bones", 
                  "muscle","heart", "spleen",
                  "gut", "liver", "kidney", "testis",
                  "ovary"
         ))


ggplot(OGG_count_final_df_wide_shared,
       aes(x = PC1,
           y = PC2,
           color=tissue)) +
  geom_point() +
  scale_color_manual(values = tissues_colors) +
  theme_classic() +
  xlab("PC1 (17.53%)") +
  ylab("PC2 (7.78%)")



#### Compute the relative expression of opsins -- Annotated species  ---------------------------------


#Define opsins clades
vo_clades <- 
  opsins_count_final_df %>% 
  filter(clade %in% c("sws1", "sws2", "rh1", "rh2", "lws")) %>%
  pull(clade) %>% unique()
nvo_clades <- 
  opsins_count_final_df %>% 
  filter(! clade %in% c("sws1", "sws2", "rh1", "rh2", "lws")) %>%
  pull(clade) %>% unique()


VO_NVO_OGG_count_final_df <- 
  OGG_count_final_df %>%
  mutate(VO_VNO = case_when(
    OGG %in% vo_clades ~ "VO",
    OGG %in% nvo_clades ~ "NVO",
    (! OGG %in% nvo_clades) & (! OGG %in% vo_clades) ~ "non_opsin",
  ))


#compute TPM on opsins vs all other genes
VO_NVO_OGG_count_final_df_summary <- 
  as.data.frame(
    VO_NVO_OGG_count_final_df %>%
  group_by(species, tissue, VO_VNO) %>%
  summarise(VO_VNO_TPM = sum(OGG_TPM))
  )


#transform TPM in proportion (/1000000)

VO_NVO_OGG_count_final_df_summary <- 
  VO_NVO_OGG_count_final_df_summary %>%
  mutate(prop_expr = VO_VNO_TPM/1000000)


#lets quickly see how to looks

VO_NVO_OGG_count_final_df_summary %>%
  filter(VO_VNO != "non_opsin") %>%
  ggplot(., aes(x=tissue, y=prop_expr, fill=VO_VNO)) +
  geom_boxplot() +
  theme_classic()

VO_NVO_OGG_count_final_df_summary %>%
  filter(VO_VNO != "non_opsin") %>%
  filter(tissue != "eye") %>%
  filter(prop_expr < 0.005) %>%
  ggplot(., aes(x=tissue, y=prop_expr, fill=VO_VNO)) +
  geom_boxplot() +
  theme_classic()


#### Is there a correlation on relative expression methods -- Annotated species  ---------------------------------

#lets compute the relative read number mapped to opsins vs non opsins ...

opsins_count_final_df_annot <- 
  opsins_count_final_df %>%
  filter(species %in% good_sp)

opsins_count_final_df_annot <- 
  opsins_count_final_df_annot %>%
  mutate(VO_VNO = case_when(
    clade %in% vo_clades ~ "VO",
    clade %in% nvo_clades ~ "NVO",
    (! clade %in% nvo_clades) & (! clade %in% vo_clades) ~ "non_opsin",
  ))

opsins_count_final_df_annot_summary <- 
  as.data.frame(opsins_count_final_df_annot %>%
  group_by(species, tissue, VO_VNO) %>%
  summarise(VO_NVO_read_nb = sum(sum_read_count)))

#Add a column with the total nb of mapped reads on the sample ...

opsins_count_final_df_annot_summary$sample <- 
  paste(opsins_count_final_df_annot_summary$species,
        opsins_count_final_df_annot_summary$tissue,
        sep='_')

number_mapped_reads_df_no <- 
  as.data.frame(number_mapped_reads_df) %>% 
  dplyr::select(sample, number_reads_mapped)

opsins_count_final_df_annot_summary_total <- 
  left_join(opsins_count_final_df_annot_summary,
            number_mapped_reads_df_no,
            by='sample')

## compute the proportion of reads mapped to opsins

opsins_count_final_df_annot_summary_total <- 
  opsins_count_final_df_annot_summary_total %>%
  mutate(prop_reads = VO_NVO_read_nb/number_reads_mapped)


## combine prop_expr and prop_read

prop_read_df <- 
  opsins_count_final_df_annot_summary_total %>%
  dplyr::select(sample, species, tissue, VO_VNO, prop_reads)

VO_NVO_OGG_count_final_df_summary$sample <- 
  paste(VO_NVO_OGG_count_final_df_summary$species,
        VO_NVO_OGG_count_final_df_summary$tissue,
        sep='_')

prop_expr_df <- 
  VO_NVO_OGG_count_final_df_summary %>%
  filter(VO_VNO != "non_opsin") %>%
  dplyr::select(sample, prop_expr, VO_VNO)

prop_expr_read_df <- 
  left_join(prop_read_df, prop_expr_df, by=c("sample", "VO_VNO"))

## Letds check if there is a correlation

## First with NVO + VO



cor_test_result <- 
  cor.test(
    prop_expr_read_df %>% 
      pull(prop_reads),
    prop_expr_read_df %>% 
      pull(prop_expr),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)


lm_read_RE <- 
  lm(data = prop_expr_read_df,
     formula = as.numeric(prop_expr) ~ as.numeric(prop_reads))
function_read_RE <- GLS_function(lm_read_RE)



prop_expr_read_df %>%
  ggplot(., aes(x=prop_reads, y=prop_expr, color=VO_VNO)) +
  geom_jitter() + 
  stat_function(fun = function_read_RE, color="black") +
  #geom_smooth(method = "lm", color = "black") +
  #scale_color_manual(values = opsins_clades_colors) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none") +
  xlab("proportion of reads") +
  ylab("proportion of expression")


## Then with VO only


cor_test_result <- 
  cor.test(
    prop_expr_read_df %>% 
      filter(VO_VNO == "VO") %>%
      pull(prop_reads),
    prop_expr_read_df %>% 
      filter(VO_VNO == "VO") %>%
      pull(prop_expr),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)



prop_expr_read_df %>%
  filter(VO_VNO == "VO") %>%
  ggplot(., aes(x=prop_reads, y=prop_expr, color=VO_VNO)) +
  geom_point() + 
  geom_smooth(method = "lm", color = "black") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none") +
  xlab("proportion of reads") +
  ylab("proportion of expression")


## Finally with NVO only

cor_test_result <- 
  cor.test(
    prop_expr_read_df %>% 
      filter(VO_VNO == "NVO") %>%
      pull(prop_reads),
    prop_expr_read_df %>% 
      filter(VO_VNO == "NVO") %>%
      pull(prop_expr),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)



prop_expr_read_df %>%
  filter(VO_VNO == "NVO") %>%
  ggplot(., aes(x=prop_reads, y=prop_expr, color=VO_VNO)) +
  geom_point() + 
  geom_smooth(method = "lm", color = "black") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none") +
  xlab("proportion of reads") +
  ylab("proportion of expression")



## Now that we have seen it was correlated, extract the slope and intercept !! 


lm_prop <- 
  lm(data = prop_expr_read_df,
     formula = prop_expr ~ prop_reads)
summary(lm_prop)

intercept_val <- summary(lm_prop)$coefficients[1]
coeff_val <- summary(lm_prop)$coefficients[2]

lm_prop_VO <- 
  lm(data = prop_expr_read_df %>% filter(VO_VNO == "VO"),
     formula = prop_expr ~ prop_reads)
summary(lm_prop_VO)

lm_prop_NVO <- 
  lm(data = prop_expr_read_df %>% filter(VO_VNO == "NVO"),
     formula = prop_expr ~ prop_reads)
summary(lm_prop_NVO)


#### Infer relative expression VO/NVO -- All species  ---------------------------------


##Import total reads mapped
number_mapped_reads_df_filter <- 
  number_mapped_reads_df %>%
  filter(! sample %in% all_bad_sample)
number_mapped_reads_df_filter <- as.data.frame(number_mapped_reads_df_filter)


##Import nb of reads on opsins

   
opsins_count_final_df_summary <- 
  opsins_count_final_df %>%
  mutate(VO_VNO = case_when(
    clade %in% vo_clades ~ "VO",
    clade %in% nvo_clades ~ "NVO",
    (! clade %in% nvo_clades) & (! clade %in% vo_clades) ~ "non_opsin",
  ))

opsins_count_final_df_summary <- 
  as.data.frame(opsins_count_final_df_summary %>%
                  group_by(species, tissue, VO_VNO) %>%
                  summarise(VO_NVO_read_nb = sum(sum_read_count)))


opsins_count_final_df_summary$sample <- 
  paste(opsins_count_final_df_summary$species,
        opsins_count_final_df_summary$tissue,
        sep='_')


##Join the two dataframes

number_mapped_reads_df_filter_no <- 
  as.data.frame(number_mapped_reads_df_filter) %>% 
  dplyr::select(sample, number_reads_mapped)

opsins_count_final_df_annot_summary_total <- 
  left_join(opsins_count_final_df_summary,
            number_mapped_reads_df_filter_no,
            by='sample')

##Compute the relative read nb

opsins_count_final_df_annot_summary_total <- 
  opsins_count_final_df_annot_summary_total %>%
  mutate(prop_reads = VO_NVO_read_nb/number_reads_mapped)

opsins_count_final_df_annot_summary_total <- 
  opsins_count_final_df_annot_summary_total %>%
  mutate(relative_expression = (0.00007 + (2.9*prop_reads)))


opsins_count_final_df_annot_summary_total$tissue <- 
  factor(opsins_count_final_df_annot_summary_total$tissue ,
         levels=c("eye", "brain", 
                  "fin", "gill", "skin","bones", 
                  "muscle","heart", "spleen",
                  "gut", "liver", "kidney", "testis",
                  "ovary"
         ))

##Plot the results ! 


opsins_count_final_df_annot_summary_total %>%
  ggplot(., aes(x=tissue, y=log(relative_expression), fill=VO_VNO)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(values = VO_NVO_colors) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  #theme(legend.position="none") +
  xlab("Tissue") +
  ylab("log(relative expression)")


#Relative expression of VO
opsins_count_final_df_annot_summary_total %>% 
  filter(VO_VNO == "VO") %>%
  group_by(tissue) %>%
  summarise(curr_mean = mean(relative_expression)) %>%
  arrange(curr_mean) %>%
  mutate(perc = curr_mean*100)
  
#Range of VO expression in the eye?
eye_VO <- 
  opsins_count_final_df_annot_summary_total %>% 
  filter(VO_VNO == "VO") %>%
  filter(tissue == "eye") %>%
  pull(relative_expression) 
mean(eye_VO)
confidence_interval(eye_VO, 0.95)


#Relative expression of NVO
opsins_count_final_df_annot_summary_total %>% 
  filter(VO_VNO == "NVO") %>%
  group_by(tissue) %>%
  summarise(curr_mean = mean(relative_expression)) %>%
  mutate(curr_mean_100 = curr_mean*100) %>%
  arrange(curr_mean) %>%
  dplyr::select(tissue, curr_mean_100)


#Range of NVO expression in the eye?
eye_NVO <- 
  opsins_count_final_df_annot_summary_total %>% 
  filter(VO_VNO == "NVO") %>%
  filter(tissue == "eye") %>%
  pull(relative_expression) 
mean(eye_NVO)*100
confidence_interval(eye_NVO, 0.95)*100

#Range in other tissues ?

ovary_NVO <- 
  opsins_count_final_df_annot_summary_total %>% 
  filter(VO_VNO == "NVO") %>%
  filter(tissue == "ovary") %>%
  pull(relative_expression) 
mean(ovary_NVO)*100
confidence_interval(ovary_NVO, 0.95)*100



#### Filter to analysis the relative expression subfamilies  ---------------------------------


toolow_reads_samples <- 
  as.data.frame(
    opsins_count_final_df %>%
  group_by(species, tissue) %>%
  summarise(total_reads = sum(sum_read_count)))



toolow_reads_samples %>%
  filter(tissue != 'eye') %>%
  filter(total_reads < 10000) %>%
  ggplot(., aes(x=tissue, y=total_reads)) +
  geom_boxplot() +
  theme_classic()
  



toolow_reads_samples %>%
  filter(total_reads > 1000) %>%
  group_by(tissue) %>%
  summarise(count = n()) 


as.data.frame(
  toolow_reads_samples %>%
    filter(total_reads > 1000) %>%
    group_by(tissue) %>%
    summarise(count = n()) 
) %>% pull(count) %>% sum()







as.data.frame(opsins_count_final_df %>%
  filter(gene_state ==  "Complete") %>%
  group_by(tissue, species) %>%
  summarise(count = n())) %>% 
  dplyr::select(species, count) %>%
  distinct() %>%
  pull(count) %>%
  mean()


toolow_reads_samples$sample <- 
  paste(toolow_reads_samples$species,
        toolow_reads_samples$tissue,
        sep="_")


good_tissues <- 
toolow_reads_samples %>%
  filter(total_reads > 1000) %>%
  group_by(tissue) %>%
  summarise(count = n())  %>%
  filter(count >= 4) %>%
  pull(tissue)

Good_Samples <- 
  toolow_reads_samples %>%
  filter(total_reads > 1000) %>%
  filter(tissue %in% good_tissues) %>%
  pull(sample)
  
 


#### Compute - Relative expression subfamilies  ---------------------------------

#Create a table with tissues of interest
opsins_count_final_df$sample <- 
  paste(opsins_count_final_df$species,
        opsins_count_final_df$tissue,
        sep="_")
Final_opsins_count_df <- 
  opsins_count_final_df %>% filter(sample %in% Good_Samples)


#### Compute the relative expression per sample


#1- Relative NVO/VO
Final_opsins_count_df <- 
  Final_opsins_count_df %>%
  mutate(VO_VNO = case_when(
    clade %in% vo_clades ~ "VO",
    clade %in% nvo_clades ~ "NVO",
    (! clade %in% nvo_clades) & (! clade %in% vo_clades) ~ "non_opsin",
  ))

NVO_VO_RE <-
  as.data.frame(
  Final_opsins_count_df %>%
  group_by(species, tissue, sample, VO_VNO) %>%
  summarise(sum_TPM = sum(TPM))
  ) %>%
  mutate(RE = sum_TPM/1000000)
  

#2- Relative among NVO

Final_opsins_count_df_NVO <- 
  Final_opsins_count_df %>%
  filter(VO_VNO == "NVO")

Among_NVO_RE <-
  as.data.frame(
    Final_opsins_count_df_NVO %>%
      group_by(species, tissue, sample, clade) %>%
      summarise(sum_TPM = sum(TPM))
  ) 

normalize_factor <- 
  as.data.frame(
    Among_NVO_RE %>%
      group_by(sample) %>%
      summarise(Total_TPM = sum(sum_TPM))
  )
  
Among_NVO_RE <- 
  as.data.frame(
  as.data.frame(
  left_join(Among_NVO_RE, normalize_factor, by="sample")) %>%
  rowwise()) %>%
  mutate(RE = sum_TPM/Total_TPM) %>%
  dplyr::select(-Total_TPM)
  
  
#3- Relative among VO - CONE vs RH1

cone_opsin_clades <- c("lws", "sws1", "sws2", "rh2")

Final_opsins_count_df_VO <- 
  Final_opsins_count_df %>%
  filter(VO_VNO == "VO")
  
Final_opsins_count_df_VO <- 
  Final_opsins_count_df_VO %>%
  mutate(Cone_vs_Rho = if_else(
    clade %in% cone_opsin_clades,
    "cone",
    "rho"
  ))


Among_VO_RE <-
  as.data.frame(
    Final_opsins_count_df_VO %>%
      group_by(species, tissue, sample, Cone_vs_Rho) %>%
      summarise(sum_TPM = sum(TPM))
  ) 

normalize_factor <- 
  as.data.frame(
    Among_VO_RE %>%
      group_by(sample) %>%
      summarise(Total_TPM = sum(sum_TPM))
  )

Among_VO_RE <- 
  as.data.frame(
    as.data.frame(
      left_join(Among_VO_RE, normalize_factor, by="sample")) %>%
      rowwise()) %>%
  mutate(RE = sum_TPM/Total_TPM) %>%
  dplyr::select(-Total_TPM)


#4  Relative among VO - cones



Final_opsins_count_df_cones <- 
  Final_opsins_count_df_VO %>%
  filter(Cone_vs_Rho == "cone")



Among_Cone_RE <-
  as.data.frame(
    Final_opsins_count_df_cones %>%
      group_by(species, tissue, sample, clade) %>%
      summarise(sum_TPM = sum(TPM))
  ) 

normalize_factor <- 
  as.data.frame(
    Among_Cone_RE %>%
      group_by(sample) %>%
      summarise(Total_TPM = sum(sum_TPM))
  )

Among_Cone_RE <- 
  as.data.frame(
    as.data.frame(
      left_join(Among_Cone_RE, normalize_factor, by="sample")) %>%
      rowwise()) %>%
  mutate(RE = sum_TPM/Total_TPM) %>%
  dplyr::select(-Total_TPM)







#### Analyse - Relative expression subfamilies -- Eye ---------------------------------

Eye_species <- Final_opsins_count_df %>% filter(tissue == "eye") %>% pull(species) %>% unique()
Eye_tree <- keep.tip(RNA_sp_tree, Eye_species)

NVO_VO_RE_eye <- NVO_VO_RE %>% filter(tissue == "eye") 
Among_VO_RE_eye <- Among_VO_RE %>% filter(tissue == "eye")
Among_Cone_RE_eye <- Among_Cone_RE %>% filter(tissue == "eye")
Among_NVO_RE_eye <- Among_NVO_RE %>% filter(tissue == "eye")

p <- 
  ggtree(Eye_tree) + 
  geom_tiplab(size=3, offset=0.1, fontface=3) + 
  theme(legend.position='none') +
  geom_treescale(fontsize=4, linesize=1, offset=0.2, x=0, y=4)




p1 <- 
  p + 
  geom_facet(panel = 'Opsins', 
             data = NVO_VO_RE_eye, geom = geom_bar, 
             mapping = aes(x = as.numeric(RE), fill = as.factor(VO_VNO)), 
             orientation = 'y', width = 0.8, stat='identity', 
             position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values=VO_NVO_colors) +
  ggnewscale::new_scale_fill() +
  geom_facet(panel = 'Visual opsins', 
             data = Among_VO_RE_eye, geom = geom_bar, 
             mapping = aes(x = as.numeric(RE), fill = as.factor(Cone_vs_Rho)), 
             orientation = 'y', width = 0.8, stat='identity', 
             position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values=Cone_vs_Rho_colors) +
  ggnewscale::new_scale_fill() +
  geom_facet(panel = 'Cone opsins', 
             data = Among_Cone_RE_eye, geom = geom_bar, 
             mapping = aes(x = as.numeric(RE), fill = as.factor(clade)), 
             orientation = 'y', width = 0.8, stat='identity', 
             position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values=Cone_colors) +
  ggnewscale::new_scale_fill() +
  geom_facet(panel = 'Non-visual opsins', 
             data = Among_NVO_RE_eye, geom = geom_bar, 
             mapping = aes(x = as.numeric(RE), fill = as.factor(clade)), 
             orientation = 'y', width = 0.8, stat='identity', 
             position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values=NVO_opsins_clades_colors) +
  ggnewscale::new_scale_fill() +
  theme_tree2(legend.position = 'none') +
  xlim_tree(800)





facet_widths(p1, c(Tree = 1.5))




## Percentage of rh1 expression

Among_VO_RE_eye %>%
  group_by(Cone_vs_Rho) %>%
  summarise(
    min_RE = min(RE),
    max_RE = max(RE),
    mean_RE = mean(RE)
  )
 
Among_VO_RE_eye %>% 
  filter(Cone_vs_Rho == "rho") %>%
  arrange(RE) %>%
  dplyr::select(species, RE)

## Is there a link between the proportion of genes and the expression ?

RH1_copy_ratio <- 
  visual_opsins_table_count_complete_wide %>%
  mutate(rh1_prop = rh1/total_visual) %>%
  mutate(cone_pop = 1-rh1_prop) %>%
  dplyr::select(species, rh1_prop, cone_pop)
  
RH1_expr_ratio <- 
  Among_VO_RE_eye %>% 
  filter(Cone_vs_Rho == "rho") %>%
  dplyr::select(species, RE)

RH1_expr_copy_ratio <-
  left_join(RH1_expr_ratio, RH1_copy_ratio, by="species")




cor_test_result <- 
  cor.test(
    RH1_expr_copy_ratio %>% 
      pull(rh1_prop),
    RH1_expr_copy_ratio %>% 
      pull(RE),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)


RH1_expr_copy_ratio %>%
  ggplot(., aes(x=rh1_prop, y=RE)) +
  geom_point() +
  geom_smooth(method = "lm", color = "black") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  xlab("Relative rh1 copy number") +
  ylab("Relative rh1 expression") +
  theme(legend.position="none")



## Now lets have some stats on cone opsins


Among_Cone_RE_eye %>%
  group_by(clade) %>%
  summarise(
    min_RE = min(RE),
    max_RE = max(RE),
    mean_RE = mean(RE)
  )

Among_Cone_RE_eye %>%
  filter(clade == "sws2") %>%
  dplyr::select(species, RE) %>%
  pull(RE) %>%
  mean()

Among_Cone_RE_eye %>%
  filter(clade == "sws1") %>%
  dplyr::select(species, RE) %>%
  pull(RE) %>%
  mean()


## Correlation with copy number ? 

Cone_copy_ratio <- 
  visual_opsins_table_count_complete_wide %>%
  mutate(rh2_prop = rh2/(lws + rh2 + sws1 + sws2)) %>%
  mutate(lws_prop = lws/(lws + rh2 + sws1 + sws2)) %>%
  mutate(sws1_prop = sws1/(lws + rh2 + sws1 + sws2)) %>%
  mutate(sws2_prop = sws2/(lws + rh2 + sws1 + sws2)) %>%
  dplyr::select(species, rh2_prop, lws_prop, sws1_prop, sws2_prop)

Cone_expr_ratio <- 
  Among_Cone_RE_eye %>%
  dplyr::select(species, clade, RE)
Cone_expr_ratio <- 
  as.data.frame(Cone_expr_ratio %>%
  pivot_wider(names_from = clade, values_from = RE, values_fill = 0))


Cone_expr_copy_ratio <-
  left_join(Cone_expr_ratio, Cone_copy_ratio, by="species")


### LWS
cor_test_result <- 
  cor.test(
    Cone_expr_copy_ratio %>% 
      pull(lws),
    Cone_expr_copy_ratio %>% 
      pull(lws_prop),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)

Cone_expr_copy_ratio$color <- 
  ifelse(Cone_expr_copy_ratio$lws == 0, "not_expressed", "expressed")

Cone_expr_copy_ratio %>%
  ggplot(., aes(x=lws_prop, y=lws, color=color)) +
  geom_jitter() +
  scale_color_manual(values = c("expressed"="black", "not_expressed"="gray")) +
  geom_smooth(method = "lm", color = "black") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  xlab("Relative lws copy number") +
  ylab("Relative lws expression") +
  theme(legend.position="none")


### RH2
cor_test_result <- 
  cor.test(
    Cone_expr_copy_ratio %>% 
      pull(rh2),
    Cone_expr_copy_ratio %>% 
      pull(rh2_prop),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)

Cone_expr_copy_ratio$color <- 
  ifelse(Cone_expr_copy_ratio$rh2 == 0, "not_expressed", "expressed")

Cone_expr_copy_ratio %>%
  ggplot(., aes(x=rh2_prop, y=rh2, color=color)) +
  geom_jitter() +
  scale_color_manual(values = c("expressed"="black", "not_expressed"="gray")) +
  geom_smooth(method = "lm", color = "black") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  xlab("Relative rh2 copy number") +
  ylab("Relative rh2 expression") +
  theme(legend.position="none")



### SWS2
cor_test_result <- 
  cor.test(
    Cone_expr_copy_ratio %>% 
      pull(sws2),
    Cone_expr_copy_ratio %>% 
      pull(sws2_prop),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)

Cone_expr_copy_ratio$color <- 
  ifelse(Cone_expr_copy_ratio$sws2 == 0, "not_expressed", "expressed")

Cone_expr_copy_ratio %>%
  ggplot(., aes(x=sws2_prop, y=sws2, color=color)) +
  scale_color_manual(values = c("expressed"="black", "not_expressed"="gray")) +
  geom_jitter() +
  geom_smooth(method = "lm", color = "black") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  xlab("Relative sws2 copy number") +
  ylab("Relative sws2 expression") +
  theme(legend.position="none")


### SWS1
cor_test_result <- 
  cor.test(
    Cone_expr_copy_ratio %>% 
      filter(sws1 < 0.3) %>%
      pull(sws1),
    Cone_expr_copy_ratio %>% 
      filter(sws1 < 0.3) %>%
      pull(sws1_prop),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)

Cone_expr_copy_ratio$color <- 
  ifelse(Cone_expr_copy_ratio$sws1 == 0, "not_expressed", "expressed")

Cone_expr_copy_ratio %>%
  #filter(sws1 < 0.3) %>%
  ggplot(., aes(x=sws1_prop, y=sws1, color=color)) +
  geom_jitter() +
  geom_smooth(method = "lm", color = "black") +
  scale_color_manual(values = c("expressed"="black", "not_expressed"="gray")) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  xlab("Relative sws1 copy number") +
  ylab("Relative sws1 expression") +
  theme(legend.position="none")



## Now lets have some stats on non visual opsins

as.data.frame(
  Among_NVO_RE_eye %>%
  group_by(clade) %>%
  summarise(mean_RE = mean(RE))) %>%
  arrange(mean_RE)

as.data.frame(
  Among_NVO_RE_eye %>%
  group_by(clade) %>%
  summarise(
    min_RE = min(RE)*100,
    max_RE = max(RE)*100,
    mean_RE = mean(RE)*100
  )) %>%
  arrange(mean_RE)


## Can I see some pinopsin :( ?


Among_NVO_RE_eye %>%
  filter(species == "Polypterus_senegalus") %>%
  dplyr::select(clade, RE) %>%
  mutate(RE_perc = RE*100) %>%
  arrange(RE)

Among_NVO_RE_eye %>%
  filter(species == "Acipenser_ruthenus") %>%
  dplyr::select(clade, RE) %>%
  mutate(RE_perc = RE*100) %>%
  arrange(RE)

Among_NVO_RE_eye %>%
  filter(species == "Amia_calva") %>%
  dplyr::select(clade, RE) %>%
  mutate(RE_perc = RE*100) %>%
  arrange(RE)

##Check if expression depends on the nb of genes ...

NVO_copy_ratio <- 
  non_visual_opsins_table_count_complete_wide %>%
  mutate(exorh_prop = exorh/total_non_visual) %>%
  mutate(opn3_prop = opn3/total_non_visual) %>%
  mutate(opn4m1_3_prop = opn4m1_3/total_non_visual) %>%
  mutate(opn4m2_prop = opn4m2/total_non_visual) %>%
  mutate(opn4x_prop = opn4x/total_non_visual) %>%
  mutate(opn5_prop = opn5/total_non_visual) %>%
  mutate(opn6_prop = opn6/total_non_visual) %>%
  mutate(opn7a_prop = opn7a/total_non_visual) %>%
  mutate(opn7b_prop = opn7b/total_non_visual) %>%
  mutate(opn8a_prop = opn8a/total_non_visual) %>%
  mutate(opn8b_prop = opn8b/total_non_visual) %>%
  mutate(opn8c_prop = opn8c/total_non_visual) %>%
  mutate(opn9_prop = opn9/total_non_visual) %>%
  mutate(parapinopsin_prop = parapinopsin/total_non_visual) %>%
  mutate(parietopsin_prop = parietopsin/total_non_visual) %>%
  mutate(rgr_prop = rgr/total_non_visual) %>%
  mutate(rrh_prop = rrh/total_non_visual) %>%
  mutate(tmt1_prop = tmt1/total_non_visual) %>%
  mutate(tmt2_prop = tmt2/total_non_visual) %>%
  mutate(tmt3_prop = tmt3/total_non_visual) %>%
  mutate(va_prop = va/total_non_visual) %>%
  mutate(pinopsin_prop = pinopsin/total_non_visual) 

test_to_do <- 
  c("exorh_prop","opn3_prop","opn4m1_3_prop","opn4m2_prop","opn4x_prop",
    "opn5_prop","opn6_prop","opn7a_prop","opn7b_prop","opn8a_prop",
    "opn8b_prop","opn8c_prop","opn9_prop","parapinopsin_prop","parietopsin_prop",
    "rgr_prop","rrh_prop","tmt1_prop","tmt2_prop","tmt3_prop","va_prop","pinopsin_prop")

NVO_copy_ratio <- 
  NVO_copy_ratio %>% dplyr::select("species", test_to_do)
  
NVO_expr_ratio <- 
  Among_NVO_RE_eye %>%
  dplyr::select(species, clade, RE)
NVO_expr_ratio <- 
  as.data.frame(NVO_expr_ratio %>%
                  pivot_wider(names_from = clade, values_from = RE, values_fill = 0))


NVO_expr_ratio <-
  left_join(NVO_expr_ratio, NVO_copy_ratio, by="species")



NVO_eye_corr_df <- as.data.frame(NULL)
for(curr_prop in test_to_do){
  
  curr_gene <- gsub("_prop","" ,curr_prop)
  
  cor_test_result <- 
    cor.test(
      NVO_expr_ratio %>% 
        pull(curr_gene),
      NVO_expr_ratio %>% 
        pull(curr_prop),
      method="pearson"
    )
  estimate_p <- round(cor_test_result$estimate, 3)
  pvalue_p <- cor_test_result$p.value
  pvalue_p <- format.pval(pvalue_p, digits = 2)
  
  curr_df <- as.data.frame(cbind(curr_gene, estimate_p, pvalue_p))
  colnames(curr_df) <- c("subfamily","Pearson's R", "p-value")
  
  NVO_eye_corr_df <- rbind(NVO_eye_corr_df, curr_df)
}

colnames(NVO_eye_corr_df) <-  c("subfamily","Pearson's R", "p-value")

gt_NVO_eye <- gt(NVO_eye_corr_df)


#### Analyse - Relative expression subfamilies -- brain ---------------------------------

brain_species <- Final_opsins_count_df %>% filter(tissue == "brain") %>% pull(species) %>% unique()
brain_tree <- keep.tip(RNA_sp_tree, brain_species)

NVO_VO_RE_brain <- NVO_VO_RE %>% filter(tissue == "brain") 
Among_VO_RE_brain <- Among_VO_RE %>% filter(tissue == "brain")
Among_Cone_RE_brain <- Among_Cone_RE %>% filter(tissue == "brain")
Among_NVO_RE_brain <- Among_NVO_RE %>% filter(tissue == "brain")

p <- 
  ggtree(brain_tree) + 
  geom_tiplab(size=3, offset=0.1, fontface=3) + 
  theme(legend.position='none') #+
  #geom_treescale(fontsize=4, linesize=1, offset=0.2, x=0, y=4)




p1 <- 
  p + 
  geom_facet(panel = 'Opsins', 
             data = NVO_VO_RE_brain, geom = geom_bar, 
             mapping = aes(x = as.numeric(RE), fill = as.factor(VO_VNO)), 
             orientation = 'y', width = 0.8, stat='identity', 
             position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values=VO_NVO_colors) +
  ggnewscale::new_scale_fill() +
  geom_facet(panel = 'Non-visual opsins', 
             data = Among_NVO_RE_brain, geom = geom_bar, 
             mapping = aes(x = as.numeric(RE), fill = as.factor(clade)), 
             orientation = 'y', width = 0.8, stat='identity', 
             position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values=NVO_opsins_clades_colors) +
  ggnewscale::new_scale_fill() +
  theme_tree2(legend.position = 'none') +
  xlim_tree(800)




facet_widths(p1, c(Tree = 1.5))



#Lets look at the most expressed NVO

Among_NVO_RE_brain %>% 
  group_by(species) %>% 
  slice(which.max(RE)) 




#What happen to exorh ? 


Among_NVO_RE_brain %>%
  filter(clade == "exorh") %>%
  mutate(RE_100 = RE*100) %>%
  arrange(RE_100) %>%
  dplyr::select(species, RE_100)


##Check if expression depends on the nb of genes ...

NVO_copy_ratio <- 
  non_visual_opsins_table_count_complete_wide %>%
  mutate(exorh_prop = exorh/total_non_visual) %>%
  mutate(opn3_prop = opn3/total_non_visual) %>%
  mutate(opn4m1_3_prop = opn4m1_3/total_non_visual) %>%
  mutate(opn4m2_prop = opn4m2/total_non_visual) %>%
  mutate(opn4x_prop = opn4x/total_non_visual) %>%
  mutate(opn5_prop = opn5/total_non_visual) %>%
  mutate(opn6_prop = opn6/total_non_visual) %>%
  mutate(opn7a_prop = opn7a/total_non_visual) %>%
  mutate(opn7b_prop = opn7b/total_non_visual) %>%
  mutate(opn8a_prop = opn8a/total_non_visual) %>%
  mutate(opn8b_prop = opn8b/total_non_visual) %>%
  mutate(opn8c_prop = opn8c/total_non_visual) %>%
  mutate(opn9_prop = opn9/total_non_visual) %>%
  mutate(parapinopsin_prop = parapinopsin/total_non_visual) %>%
  mutate(parietopsin_prop = parietopsin/total_non_visual) %>%
  mutate(rgr_prop = rgr/total_non_visual) %>%
  mutate(rrh_prop = rrh/total_non_visual) %>%
  mutate(tmt1_prop = tmt1/total_non_visual) %>%
  mutate(tmt2_prop = tmt2/total_non_visual) %>%
  mutate(tmt3_prop = tmt3/total_non_visual) %>%
  mutate(va_prop = va/total_non_visual) %>%
  mutate(pinopsin_prop = pinopsin/total_non_visual) 

test_to_do <- 
  c("exorh_prop","opn3_prop","opn4m1_3_prop","opn4m2_prop","opn4x_prop",
    "opn5_prop","opn6_prop","opn7a_prop","opn7b_prop","opn8a_prop",
    "opn8b_prop","opn8c_prop","opn9_prop","parapinopsin_prop","parietopsin_prop",
    "rgr_prop","rrh_prop","tmt1_prop","tmt2_prop","tmt3_prop","va_prop","pinopsin_prop")

NVO_copy_ratio <- 
  NVO_copy_ratio %>% dplyr::select("species", test_to_do)

NVO_expr_ratio <- 
  Among_NVO_RE_brain %>%
  dplyr::select(species, clade, RE)
NVO_expr_ratio <- 
  as.data.frame(NVO_expr_ratio %>%
                  pivot_wider(names_from = clade, values_from = RE, values_fill = 0))


NVO_expr_ratio <-
  left_join(NVO_expr_ratio, NVO_copy_ratio, by="species")



NVO_brain_corr_df <- as.data.frame(NULL)
for(curr_prop in test_to_do){
  
  curr_gene <- gsub("_prop","" ,curr_prop)
  
  cor_test_result <- 
    cor.test(
      NVO_expr_ratio %>% 
        pull(curr_gene),
      NVO_expr_ratio %>% 
        pull(curr_prop),
      method="pearson"
    )
  estimate_p <- round(cor_test_result$estimate, 3)
  pvalue_p <- cor_test_result$p.value
  pvalue_p <- format.pval(pvalue_p, digits = 2)
  
  curr_df <- as.data.frame(cbind(curr_gene, estimate_p, pvalue_p))
  colnames(curr_df) <- c("subfamily","Pearson's R", "p-value")
  
  NVO_brain_corr_df <- rbind(NVO_brain_corr_df, curr_df)
}

colnames(NVO_brain_corr_df) <-  c("subfamily","Pearson's R", "p-value")

gt_NVO_brain <- gt(NVO_brain_corr_df)


#### Analyse - Relative expression subfamilies -- heart ---------------------------------

heart_species <- Final_opsins_count_df %>% filter(tissue == "heart") %>% pull(species) %>% unique()
heart_tree <- keep.tip(RNA_sp_tree, heart_species)

NVO_VO_RE_heart <- NVO_VO_RE %>% filter(tissue == "heart") 
Among_VO_RE_heart <- Among_VO_RE %>% filter(tissue == "heart")
Among_Cone_RE_heart <- Among_Cone_RE %>% filter(tissue == "heart")
Among_NVO_RE_heart <- Among_NVO_RE %>% filter(tissue == "heart")

p <- 
  ggtree(heart_tree) + 
  geom_tiplab(size=3, offset=0.1, fontface=3) + 
  theme(legend.position='none') #+
#geom_treescale(fontsize=4, linesize=1, offset=0.2, x=0, y=4)




p1 <- 
  p + 
  geom_facet(panel = 'Opsins', 
             data = NVO_VO_RE_heart, geom = geom_bar, 
             mapping = aes(x = as.numeric(RE), fill = as.factor(VO_VNO)), 
             orientation = 'y', width = 0.8, stat='identity', 
             position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values=VO_NVO_colors) +
  ggnewscale::new_scale_fill() +
  geom_facet(panel = 'Non-visual opsins', 
             data = Among_NVO_RE_heart, geom = geom_bar, 
             mapping = aes(x = as.numeric(RE), fill = as.factor(clade)), 
             orientation = 'y', width = 0.8, stat='identity', 
             position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values=NVO_opsins_clades_colors) +
  ggnewscale::new_scale_fill() +
  theme_tree2(legend.position = 'none') +
  xlim_tree(800)




facet_widths(p1, c(Tree = 1.5))



#facet_widths(p2, widths = c(1, 1))


Final_opsins_count_df %>%
  filter(species == "Poecilia_reticulata") %>%
  pull(TPM) %>%
  sum()

as.data.frame(
  Final_opsins_count_df %>%
    group_by(species, tissue, sample, clade) %>%
    summarise(sum_TPM = sum(TPM))
) %>%
  mutate(RE = sum_TPM/1000000)  %>%
  filter(species == "Gadus_morhua") %>%
  filter(tissue == "ovary") %>%
  dplyr::select(clade, RE)


as.data.frame(
  Final_opsins_count_df %>%
    group_by(species, tissue, sample, clade) %>%
    summarise(sum_TPM = sum(TPM))
) %>%
  mutate(RE = sum_TPM/1000000)  %>%
  filter(species == "Amia_calva") %>%
  filter(tissue == "liver") %>%
  dplyr::select(clade, RE)

#### Analyse - Relative expression subfamilies -- liver ---------------------------------

liver_species <- Final_opsins_count_df %>% filter(tissue == "liver") %>% pull(species) %>% unique()
liver_tree <- keep.tip(RNA_sp_tree, liver_species)

NVO_VO_RE_liver <- NVO_VO_RE %>% filter(tissue == "liver") 
Among_VO_RE_liver <- Among_VO_RE %>% filter(tissue == "liver")
Among_Cone_RE_liver <- Among_Cone_RE %>% filter(tissue == "liver")
Among_NVO_RE_liver <- Among_NVO_RE %>% filter(tissue == "liver")

p <- 
  ggtree(liver_tree) + 
  geom_tiplab(size=3, offset=0.1, fontface=3) + 
  theme(legend.position='none') #+
#geom_treescale(fontsize=4, linesize=1, offset=0.2, x=0, y=4)




p1 <- 
  p + 
  geom_facet(panel = 'Opsins', 
             data = NVO_VO_RE_liver, geom = geom_bar, 
             mapping = aes(x = as.numeric(RE), fill = as.factor(VO_VNO)), 
             orientation = 'y', width = 0.8, stat='identity', 
             position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values=VO_NVO_colors) +
  ggnewscale::new_scale_fill() +
  geom_facet(panel = 'Non-visual opsins', 
             data = Among_NVO_RE_liver, geom = geom_bar, 
             mapping = aes(x = as.numeric(RE), fill = as.factor(clade)), 
             orientation = 'y', width = 0.8, stat='identity', 
             position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values=NVO_opsins_clades_colors) +
  ggnewscale::new_scale_fill() +
  theme_tree2(legend.position = 'none') +
  xlim_tree(800)





facet_widths(p1, c(Tree = 1.5))




####

Among_NVO_RE_heart %>% 
  group_by(species) %>% 
  slice(which.max(RE)) 

Among_NVO_RE_liver %>% 
  group_by(species) %>% 
  slice(which.max(RE)) 


Among_NVO_RE_liver %>%
  filter(species == "Gasterosteus_aculeatus_aculeatus") %>%
  filter(RE > 0.3)



#### Analyse - Relative expression subfamilies -- ovary ---------------------------------

ovary_species <- Final_opsins_count_df %>% filter(tissue == "ovary") %>% pull(species) %>% unique()
ovary_tree <- keep.tip(RNA_sp_tree, ovary_species)

NVO_VO_RE_ovary <- NVO_VO_RE %>% filter(tissue == "ovary") 
Among_VO_RE_ovary <- Among_VO_RE %>% filter(tissue == "ovary")
Among_Cone_RE_ovary <- Among_Cone_RE %>% filter(tissue == "ovary")
Among_NVO_RE_ovary <- Among_NVO_RE %>% filter(tissue == "ovary")

p <- 
  ggtree(ovary_tree) + 
  geom_tiplab(size=3, offset=0.1, fontface=3) + 
  theme(legend.position='none') #+
#geom_treescale(fontsize=4, linesize=1, offset=0.2, x=0, y=4)




p1 <- 
  p + 
  geom_facet(panel = 'Opsins', 
             data = NVO_VO_RE_ovary, geom = geom_bar, 
             mapping = aes(x = as.numeric(RE), fill = as.factor(VO_VNO)), 
             orientation = 'y', width = 0.8, stat='identity', 
             position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values=VO_NVO_colors) +
  ggnewscale::new_scale_fill() +
  geom_facet(panel = 'Non-visual opsins', 
             data = Among_NVO_RE_ovary, geom = geom_bar, 
             mapping = aes(x = as.numeric(RE), fill = as.factor(clade)), 
             orientation = 'y', width = 0.8, stat='identity', 
             position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values=NVO_opsins_clades_colors) +
  ggnewscale::new_scale_fill() +
  theme_tree2(legend.position = 'none') +
  xlim_tree(800)



facet_widths(p1, c(Tree = 1.5))





#### Analyse - Relative expression subfamilies -- testis ---------------------------------

testis_species <- Final_opsins_count_df %>% filter(tissue == "testis") %>% pull(species) %>% unique()
testis_tree <- keep.tip(RNA_sp_tree, testis_species)

NVO_VO_RE_testis <- NVO_VO_RE %>% filter(tissue == "testis") 
Among_VO_RE_testis <- Among_VO_RE %>% filter(tissue == "testis")
Among_Cone_RE_testis <- Among_Cone_RE %>% filter(tissue == "testis")
Among_NVO_RE_testis <- Among_NVO_RE %>% filter(tissue == "testis")

p <- 
  ggtree(testis_tree) + 
  geom_tiplab(size=3, offset=0.1, fontface=3) + 
  theme(legend.position='none') #+




p1 <- 
  p + 
  geom_facet(panel = 'Opsins', 
             data = NVO_VO_RE_testis, geom = geom_bar, 
             mapping = aes(x = as.numeric(RE), fill = as.factor(VO_VNO)), 
             orientation = 'y', width = 0.8, stat='identity', 
             position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values=VO_NVO_colors) +
  ggnewscale::new_scale_fill() +
  geom_facet(panel = 'Non-visual opsins', 
             data = Among_NVO_RE_testis, geom = geom_bar, 
             mapping = aes(x = as.numeric(RE), fill = as.factor(clade)), 
             orientation = 'y', width = 0.8, stat='identity', 
             position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values=NVO_opsins_clades_colors) +
  new_scale_fill() +
  theme_tree2(legend.position = 'none') +
  xlim_tree(800)





facet_widths(p1, c(Tree = 1.5))




#### Boxplot NVO expression per tissue per opsin ---------------------------------



Among_NVO_RE_eye %>%
  ggplot(., aes(x=reorder(clade, RE, mean), y=RE)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  xlab("") +
  ylab("Relative expression") +
  theme(legend.position="none")

Among_NVO_RE_brain %>%
  ggplot(., aes(x=reorder(clade, RE, mean), y=RE)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  xlab("") +
  ylab("Relative expression") +
  theme(legend.position="none")


Among_NVO_RE_ovary %>%
  ggplot(., aes(x=reorder(clade, RE, mean), y=RE)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  xlab("") +
  ylab("Relative expression") +
  theme(legend.position="none")


Among_NVO_RE_testis %>%
  ggplot(., aes(x=reorder(clade, RE, mean), y=RE)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  xlab("") +
  ylab("Relative expression") +
  theme(legend.position="none")



as.data.frame(
  Among_NVO_RE_eye %>% group_by(clade) %>% summarise(mean_ops = mean(RE)) %>% arrange(mean_ops)) %>%
  arrange(mean_ops)

as.data.frame(
  Among_NVO_RE_brain %>% group_by(clade) %>% summarise(mean_ops = mean(RE)) %>% arrange(mean_ops)) %>%
  arrange(mean_ops)

as.data.frame(
  Among_NVO_RE_ovary %>% group_by(clade) %>% summarise(mean_ops = mean(RE)) %>% arrange(mean_ops)) %>%
  arrange(mean_ops)

as.data.frame(
  Among_NVO_RE_testis %>% group_by(clade) %>% summarise(mean_ops = mean(RE)) %>% arrange(mean_ops)) %>%
  arrange(mean_ops)






Among_NVO_RE_eye %>%
  filter(species %in% c("Polypterus_senegalus", "Acipenser_ruthenus", "Amia_calva")) %>%
  ggplot(., aes(x=reorder(clade, RE, mean), y=RE)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  xlab("") +
  ylab("Relative expression") +
  theme(legend.position="none")




as.data.frame(
  Among_NVO_RE_eye %>% 
    filter(species %in% c("Polypterus_senegalus", "Acipenser_ruthenus", "Amia_calva")) %>%
    group_by(clade) %>% summarise(mean_ops = mean(RE)) %>% arrange(mean_ops)) %>%
  arrange(mean_ops)


#### Compute Shannon diversity index ---------------------------------

#https://www.omnicalculator.com/ecology/shannon-index
#library(vegan)

all_samples <- Among_NVO_RE %>% pull(sample) %>% unique()

Shannon_index_df <- as.data.frame(NULL)
for(curr_sample in all_samples){
  
  
  curr_sp <- 
    Among_NVO_RE %>%
    filter(sample == curr_sample) %>%
    pull(species) %>%
    unique()
  
  curr_tissue <- 
    Among_NVO_RE %>%
    filter(sample == curr_sample) %>%
    pull(tissue) %>%
    unique()
    
  Shannon_H <- 
    Among_NVO_RE %>%
    filter(sample == curr_sample) %>%
    mutate(ln_RE = log(RE)) %>%
    mutate(ln_RE_x_RE = ln_RE * RE) %>%
    mutate(min_ln_RE_x_RE = -ln_RE_x_RE) %>%
    pull(min_ln_RE_x_RE) %>%
    sum(na.rm=TRUE)
  
  
  curr_df <-
    as.data.frame(
    cbind(
      curr_sample,
      curr_sp,
      curr_tissue,
      Shannon_H
    )
  )
  colnames(curr_df) <- c("sample", "species", "tissue", "Shannon_H")
  
  Shannon_index_df <- rbind(Shannon_index_df, curr_df)
}

Shannon_index_df$Shannon_H <- as.numeric(Shannon_index_df$Shannon_H )


Shannon_index_df %>%
  ggplot(., aes(x=reorder(tissue,Shannon_H, mean) , y=Shannon_H)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  xlab("") +
  ylab("Shannon diversity index") +
  theme(legend.position="none")




Shannon_index_df %>%
  group_by(tissue) %>%
  summarise(mean_H = mean(Shannon_H))
  
  
  
  

#### Megalops cyprinoides -- brain ---------------------------------


Megalops_count_df <- 
  opsins_count_final_df %>% filter(sample %in% "Megalops_cyprinoides_brain")


Megalops_count_df %>% pull(sum_read_count) %>% sum()


NVO_VO_RE_brain <- NVO_VO_RE %>% filter(tissue == "brain") 
Among_VO_RE_brain <- Among_VO_RE %>% filter(tissue == "brain")
Among_Cone_RE_brain <- Among_Cone_RE %>% filter(tissue == "brain")
Among_NVO_RE_brain <- Among_NVO_RE %>% filter(tissue == "brain")

#### Correlation dN/dS per OGG and TPM values ?? gene level ---------------------------------

#Load dN/dS per OGG per species table

omega_per_sp_per_gene_df_wide <- 
  read.table("OGG_omega_table_sp.csv",
            sep=",",
            header=TRUE)




all_genes_count_final_df$sample <- 
  paste(all_genes_count_final_df$species,
        all_genes_count_final_df$tissue,
        sep="_")


#Test for each sample ... 

unique_sample <- OGG_count_final_df_wide %>% pull(sample) %>% unique()
all_plot_names_gene <- c()
plot_list_gene <- list()
Pearson_omega_tpm_df_gene <- as.data.frame(NULL)
for(curr_sample in unique_sample){
  
  curr_species <- 
    all_genes_count_final_df %>%
    filter(sample == curr_sample) %>%
    pull(species) %>% unique()
  
  curr_tissue <- 
    as.character(OGG_count_final_df_wide_shared %>%
                   filter(sample == curr_sample) %>%
                   pull(tissue) %>% unique())
  
  
  curr_all_genes_count_final_df <- 
    all_genes_count_final_df %>%
    filter(sample == curr_sample)
  
  curr_all_genes_count_final_df <- 
    curr_all_genes_count_final_df %>%
    mutate(gene_name_parsed = paste(paste(curr_species, "_rna_", sep=""), gene_name, sep=""))
  
  curr_all_genes_count_final_df$gene_name_parsed <- 
    gsub("\\.", "_", curr_all_genes_count_final_df$gene_name_parsed)
  
  

  
  curr_all_genes_count_final_df <- 
    curr_all_genes_count_final_df %>% 
    dplyr::select(species, gene_name_parsed, tissue, sample, TPM)
  
  colnames(curr_all_genes_count_final_df) <- 
    c("species", "gene_name", "tissue", "sample", "TPM")
  
  
  curr_sample_omegas <- 
    omega_per_sp_per_gene_df_wide %>%
    filter(species == curr_species)
  
  
  if(curr_species == "Albula_goreensis"){
    curr_sample_omegas$gene_name <-
      gsub("_gnl_WGS_JAERUA", "", curr_sample_omegas$gene_name)
  } else if (curr_species == "Aldrovandia_affinis"){
    curr_sample_omegas$gene_name <-
      gsub("_gnl_WGS_JAINUG", "", curr_sample_omegas$gene_name)
  } else if (curr_species == "Atractosteus_spatula"){
    curr_sample_omegas$gene_name <-
      gsub("_gnl_WGS_JAAWVO", "", curr_sample_omegas$gene_name)
  } else if (curr_species == "Collichthys_lucidus"){
    curr_sample_omegas$gene_name <-
      gsub("_gnl_WGS_SCMI", "", curr_sample_omegas$gene_name)
  } else if (curr_species == "Conger_conger"){
    curr_sample_omegas$gene_name <-
      gsub("_gnl_WGS_JAFJMO", "", curr_sample_omegas$gene_name)
  } else if (curr_species == "Megalops_atlanticus"){
    curr_sample_omegas$gene_name <-
      gsub("_gnl_WGS_JAFDVH", "", curr_sample_omegas$gene_name)
  } else if (curr_species == "Nibea_albiflora"){
    curr_sample_omegas$gene_name <-
      gsub("_gnl_WGS_JABGLX", "", curr_sample_omegas$gene_name)
  } else if (curr_species == "Silurus_asotus"){
    curr_sample_omegas$gene_name <-
      gsub("_gnl_WGS_QPAB", "", curr_sample_omegas$gene_name)
  }
  

  curr_TPM_omega <- 
    as.data.frame(
      left_join(curr_all_genes_count_final_df, curr_sample_omegas, by=c("species", "gene_name")))
  
  curr_TPM_omega <- 
    curr_TPM_omega %>% filter(! is.na(MLE)) %>% filter(! is.na(TPM))
  
  print(c(curr_species, nrow(curr_TPM_omega)))
  
  curr_lm <- 
    lm(data = curr_TPM_omega,
       formula = as.numeric(TPM) ~ as.numeric(MLE))
  curr_function <- GLS_function(curr_lm)
  
  
  p = curr_TPM_omega %>%
    ggplot(., aes(x=MLE, y=TPM)) +
    geom_point() +
    theme_classic() +
    stat_function(fun = curr_function, color="black") +
    #labs(subtitle = bquote(R^2 ~ "=" ~ .(curr_R2) ~ ";" ~ FDR.pvalue ~ "=" ~.(curr_pvalue) ~ ";" ~ lambda ~ "=" ~ .(curr_lambda))) +
    xlab("dN/dS") +
    ylab("TPM") +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          plot.subtitle=element_text(size=16),
          legend.position="none") 
  
  
  plot_list_gene[[curr_sample]] = p
  all_plot_names_gene <- c(all_plot_names_gene, curr_sample)
  
  
  
  cor_test_result <- 
    cor.test(
      curr_TPM_omega %>% pull(TPM),
      curr_TPM_omega %>% pull(MLE),
      method="pearson"
    )
  estimate_p <- round(cor_test_result$estimate, 3)
  pvalue_p <- cor_test_result$p.value
  
  
  curr_pearson_df <- 
    as.data.frame(
      cbind(curr_species, curr_tissue, curr_sample, estimate_p, pvalue_p)
    )
  colnames(curr_pearson_df) <- c("species", "tissue","sample", "Pearson_R", "pvalue")
  
  Pearson_omega_tpm_df_gene <- 
    rbind(Pearson_omega_tpm_df_gene, curr_pearson_df)
  
}





Pearson_omega_tpm_df_gene$Pearson_R <- as.numeric(Pearson_omega_tpm_df_gene$Pearson_R)
Pearson_omega_tpm_df_gene$pvalue <- as.numeric(Pearson_omega_tpm_df_gene$pvalue)




Pearson_omega_tpm_df_gene <- 
  read.table("Pearson_omega_tpm_df_gene.csv", 
             sep=",",
             header=TRUE)


nrow(Pearson_omega_tpm_df_gene %>% filter(pvalue < 0.05))

Pearson_omega_tpm_df_gene <- 
  Pearson_omega_tpm_df_gene %>%
  mutate(signif = if_else(
    pvalue < 0.05,
    "Signif",
    "Non_signif"
  ))


order_tissue <- as.data.frame(Pearson_omega_tpm_df_gene %>%
                                group_by(tissue, signif) %>%
                                summarise(count = n())) %>%
  pivot_wider(names_from = signif, values_from = count) %>%
  mutate(Non_signif = coalesce(Non_signif, 0)) %>%
  mutate(ratio = Signif/(Signif + Non_signif)) %>%
  arrange(desc(ratio)) %>%
  pull(tissue)

Pearson_omega_tpm_df_gene$tissue <-
  factor(Pearson_omega_tpm_df_gene$tissue,
         levels=order_tissue)



Pearson_omega_tpm_df_gene %>%
  ggplot(., aes(x=tissue, y=Pearson_R)) +
  geom_violin() +
  geom_jitter(shape=16, position=position_jitter(0.2), aes(color=signif)) +
  scale_color_manual(values = c("Signif" = "#BB9626", "Non_signif" = "#1A71BD")) +
  xlab("") +
  ylab("Pearson's R") +
  theme_classic() + 
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") 



as.data.frame(Pearson_omega_tpm_df_gene %>%
                group_by(tissue, signif) %>%
                summarise(count = n())) %>%
  pivot_wider(names_from = signif, values_from = count) %>%
  mutate(Non_signif = coalesce(Non_signif, 0)) %>%
  mutate(ratio = Signif/(Signif + Non_signif),
         total = Signif + Non_signif) %>%
  dplyr::select(tissue, Signif, total)



#### Correlation dN/dS per OGG and TPM values ?? OGG level ---------------------------------

#Load dN/dS per OGG per species table

mean_omega_per_sp_per_ogg_df_wide <- 
  read.table("mean_omega_per_sp_per_ogg_df_wide.tsv",
             header=TRUE,
             sep=",")


#Test for each sample ... 

unique_sample <- OGG_count_final_df_wide %>% pull(sample) %>% unique()
all_plot_names_all <- c()
plot_list_all <- list()
Pearson_omega_tpm_df_all <- as.data.frame(NULL)
for(curr_sample in unique_sample){
  
  curr_species <- 
    OGG_count_final_df_wide %>%
    filter(sample == curr_sample) %>%
    pull(species) %>% unique()
  
  curr_tissue <- 
    as.character(OGG_count_final_df_wide_shared %>%
    filter(sample == curr_sample) %>%
    pull(tissue) %>% unique())
  
  curr_TPM_df <- 
    OGG_count_final_df_wide %>%
    filter(sample == curr_sample)
  
  curr_TPM_df_long <- 
    as.data.frame(
      curr_TPM_df %>%
        pivot_longer(!c(species, tissue, sample), names_to = "OGG", values_to = "TPM"))
  
  
  curr_sample_omegas <- 
    mean_omega_per_sp_per_ogg_df_wide %>%
    filter(species == curr_species)
  
  curr_sample_omegas_long <- 
    as.data.frame(
      curr_sample_omegas %>%
        pivot_longer(!c(species), names_to = "OGG", values_to = "omega"))
  
  
  curr_TPM_omega <- 
    as.data.frame(
      left_join(curr_TPM_df_long, curr_sample_omegas_long, by=c("species", "OGG")))
  
  
  curr_lm <- 
    lm(data = curr_TPM_omega %>% filter(! is.na(omega)) %>% filter(! is.na(TPM)),
       formula = as.numeric(TPM) ~ as.numeric(omega))
  curr_function <- GLS_function(curr_lm)
  
  
  p = curr_TPM_omega %>%
    ggplot(., aes(x=omega, y=TPM)) +
    geom_point() +
    theme_classic() +
    stat_function(fun = curr_function, color="black") +
    #labs(subtitle = bquote(R^2 ~ "=" ~ .(curr_R2) ~ ";" ~ FDR.pvalue ~ "=" ~.(curr_pvalue) ~ ";" ~ lambda ~ "=" ~ .(curr_lambda))) +
    xlab("dN/dS") +
    ylab("TPM") +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          plot.subtitle=element_text(size=16),
          legend.position="none") 
  
  
  plot_list_all[[curr_sample]] = p
  all_plot_names_all <- c(all_plot_names_all, curr_sample)
  
  
  
  cor_test_result <- 
    cor.test(
      curr_TPM_omega %>% filter(! is.na(omega)) %>% filter(! is.na(TPM)) %>% pull(TPM),
      curr_TPM_omega %>% filter(! is.na(omega)) %>% filter(! is.na(TPM)) %>% pull(omega),
      method="pearson"
    )
  estimate_p <- round(cor_test_result$estimate, 3)
  pvalue_p <- cor_test_result$p.value

  
  curr_pearson_df <- 
    as.data.frame(
      cbind(curr_species, curr_tissue, curr_sample, estimate_p, pvalue_p)
    )
  colnames(curr_pearson_df) <- c("species", "tissue","sample", "Pearson_R", "pvalue")
  
  Pearson_omega_tpm_df_all <- 
    rbind(Pearson_omega_tpm_df_all, curr_pearson_df)
  
}





Pearson_omega_tpm_df_all$Pearson_R <- as.numeric(Pearson_omega_tpm_df_all$Pearson_R)
Pearson_omega_tpm_df_all$pvalue <- as.numeric(Pearson_omega_tpm_df_all$pvalue)




Pearson_omega_tpm_df_all <- 
  read.table("Pearson_omega_tpm_df_all.csv", 
             sep=",",
             header=TRUE)


nrow(Pearson_omega_tpm_df_all %>% filter(pvalue < 0.05))

Pearson_omega_tpm_df_all <- 
  Pearson_omega_tpm_df_all %>%
  mutate(signif = if_else(
    pvalue < 0.05,
    "Signif",
    "Non_signif"
  ))


order_tissue <- as.data.frame(Pearson_omega_tpm_df_all %>%
  group_by(tissue, signif) %>%
  summarise(count = n())) %>%
  pivot_wider(names_from = signif, values_from = count) %>%
  mutate(ratio = Signif/(Signif + Non_signif)) %>%
  arrange(desc(ratio)) %>%
  pull(tissue)

Pearson_omega_tpm_df_all$tissue <-
  factor(Pearson_omega_tpm_df_all$tissue,
         levels=order_tissue)



Pearson_omega_tpm_df_all %>%
  ggplot(., aes(x=tissue, y=Pearson_R)) +
  geom_violin() +
  geom_jitter(shape=16, position=position_jitter(0.2), aes(color=signif)) +
  scale_color_manual(values = c("Signif" = "#BB9626", "Non_signif" = "#1A71BD")) +
  xlab("") +
  ylab("Pearson's R") +
  theme_classic() + 
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") 



  
#### Comparison dN/dS gene level --  OGG level ---------------------------------

Pearson_omega_tpm_df_gene <- 
  read.table("Pearson_omega_tpm_df_gene.csv", 
             sep=",",
             header=TRUE)

Pearson_omega_tpm_df_all <- 
  read.table("Pearson_omega_tpm_df_all.csv", 
             sep=",",
             header=TRUE)



Pearson_omega_tpm_df_gene_all <- 
  left_join(Pearson_omega_tpm_df_gene, Pearson_omega_tpm_df_all,
          by=c("species", "tissue", "sample"),
          suffix = c(".gene", ".OGG"))




Pearson_omega_tpm_df_gene_all <- 
  Pearson_omega_tpm_df_gene_all %>%
  mutate(Signif = case_when(
    (pvalue.OGG < 0.05) & (pvalue.gene < 0.05) ~ "signif_both", 
    (pvalue.OGG >= 0.05) & (pvalue.gene < 0.05) ~ "signif_gene_only", 
    (pvalue.OGG < 0.05) & (pvalue.gene >= 0.05) ~ "signif_OGG_only", 
    (pvalue.OGG >= 0.05) & (pvalue.gene >= 0.05) ~ "signif_none"
  ))





cor_test_result <- 
  cor.test(
    Pearson_omega_tpm_df_gene_all %>% pull(Pearson_R.OGG),
    Pearson_omega_tpm_df_gene_all %>% pull(Pearson_R.gene),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value


Pearson_omega_tpm_df_gene_all %>%
  ggplot(., aes(x=Pearson_R.gene, y=Pearson_R.OGG, color=Signif)) +
  geom_point() +
  theme_classic() +
  geom_abline(slope=1, linetype="dashed", color = "gray") +
  geom_smooth(method = "lm", color = "black") +
  scale_color_manual(values = 
                       c("signif_both" = "#BB9626", "signif_gene_only" = "#DC267F",
                         "signif_OGG_only" = "black", "signif_none" = "#1A71BD")) +
  xlab("Pearson's R - Gene level") +
  ylab("Pearson's R - OGG level") +
  theme_classic() + 
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") 





