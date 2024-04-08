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
library(grid)
library(phangorn)

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



##### BUSCO species filter and Figure  ---------------------------------


teleost_B90 <- 
  species_table %>%
  filter(infraclass == "Teleostei") %>%
  filter(perc_BUSCO >= 0.9) %>%
  pull(genome_name)

nonteleost_B80 <- 
  species_table %>%
  filter(infraclass != "Teleostei") %>%
  filter(perc_BUSCO >= 0.8) %>%
  pull(genome_name)


species_table <- 
  species_table %>%
  mutate(Teleost_or_not = if_else(
    infraclass == "Teleostei",
    "Teleost",
    "Non-teleost"
  ))

final_dataset <- c(teleost_B90, nonteleost_B80)



#put species infraclass in a dataframe 
infraclass_df <- as.data.frame(enframe(c(species_table %>% pull(Teleost_or_not) %>% unique(), "Total")) %>% dplyr::select(value))
colnames(infraclass_df) <- "Teleost_or_not"

#Loop over every BUSCO score from 0 to 1 and compute the number of species filtered at each score
BUSCO_erosion_df <- c()
for (i in seq(0, 1, 0.01)){
  
  current_df <- 
    species_table %>% 
    dplyr::filter(perc_BUSCO >= i) %>% 
    group_by(Teleost_or_not) %>% 
    summarise(number = n()) %>%
    adorn_totals("row") %>%
    mutate(Busco_percentage = i)
  
  
  merged_df <- left_join(infraclass_df, current_df, by="Teleost_or_not")
  merged_df$number[is.na(merged_df$number)] <- 0
  merged_df$Busco_percentage[is.na(merged_df$Busco_percentage)] <- i
  
  BUSCO_erosion_df <- rbind(BUSCO_erosion_df, merged_df)
  
  
}



#Define the number of species filtered 
BUSCO_90_nb_teleost <- BUSCO_erosion_df %>% filter(Busco_percentage == 0.8) %>% filter(Teleost_or_not == "Teleost") %>% pull(number)
BUSCO_80_nb_nonteleost <- BUSCO_erosion_df %>% filter(Busco_percentage == 0.85) %>% filter(Teleost_or_not == "Non-teleost") %>% pull(number)


BUSCO_erosion_df %>%
  ggplot(., aes(x=Busco_percentage, y=number, color=Teleost_or_not)) +
  geom_point() +
  theme_classic() +
  geom_line() +
  scale_color_manual(values = c("#FFB000", "#DC267F", "gray")) + 
  ylab("Number of genomes") +
  xlab("Proportion of complete BUSCO genes")+
  #geom_segment(aes(x = 0.8, y = 0, xend = 0.8, yend = BUSCO_90_nb_teleost),color="#DC267F", linetype="dashed") +
  #geom_segment(aes(x = 0, y = BUSCO_90_nb_teleost, xend = 0.8, yend = BUSCO_90_nb_teleost),color="#DC267F", linetype="dashed") +
  #geom_segment(aes(x = 0.8, y = 0, xend = 0.8, yend = BUSCO_80_nb_nonteleost),color="#FFB000", linetype="dashed") +
  #geom_segment(aes(x = 0, y = BUSCO_80_nb_nonteleost, xend = 0.8, yend = BUSCO_80_nb_nonteleost),color="#FFB000", linetype="dashed") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none")





##### Data load - Sequencing technology -----------------

sequencing_techno_df <-
  read.table("Sequencing_tech_simplified.csv",
             header=FALSE,
             sep=",")
colnames(sequencing_techno_df) <- c("assembly_accession", "techno")

sequencing_techno_df[sequencing_techno_df == ''] <- NA



##### Data load - Sequencing date -----------------

date_assembly_df <-
  read.table("Assembly_dates.csv",
             header=FALSE,
             sep=",")
colnames(date_assembly_df) <- c("assembly_accession", "date")


##### Data load - FishBase table -----------------

fishbase_info <-
  read.table("FishBaseTable.tsv",
             header=TRUE,
             sep="\t", 
             quote = "")


### Redo some columns for further pgls analysis .. 

#Fresh or Salt Water
fishbase_info <-
  fishbase_info %>%
  mutate(Fresh_salt = case_when(
    (Fresh == 1  & Brack == 1 & Saltwater == 1) ~ "Fresh_Brack_Salt",
    (Fresh == 1  & Brack == 1 & Saltwater == 0) ~ "Fresh_Brack",
    (Fresh == 0  & Brack == 1 & Saltwater == 1) ~ "Salt_Brack",
    (Fresh == 1  & Brack == 0 & Saltwater == 1) ~ "Fresh_Salt",
    (Fresh == 1  & Brack == 0 & Saltwater == 0) ~ "Fresh",
    (Fresh == 0  & Brack == 1 & Saltwater == 0) ~ "Brack",
    (Fresh == 0  & Brack == 0 & Saltwater == 1) ~ "Salt",
  ))

#Cave species. Remove poecilia mexicana

fishbase_info <-
  fishbase_info %>%
  mutate(Cave_Surface = if_else(
    Cave == -1,
    "Cave",
    "Surface"
  ))

fishbase_info %>%
  filter(Cave_Surface == "Cave") %>%
  pull(my_species_name)

fishbase_info[(fishbase_info$my_species_name == "Poecilia mexicana"),"Cave_Surface"] <- "Surface"


fishbase_info <- 
  fishbase_info %>%
  mutate(Cave_Surface_filter = case_when(
    Cave_Surface == "Surface" ~ "Surface", 
    Cave_Surface == "Cave" ~ "Cave", 
    is.na(Cave_Surface) ~ "Surface"
  )) 


#Relative eye size (ED / SL) + absolute eye size 

fishbase_info %>%
  filter(my_species_name == "Amia calva") %>%
  pull(CommonLength) 

fishbase_info$TL <- as.numeric(fishbase_info$TL) 

fishbase_info <- 
  fishbase_info %>%
  mutate(ED_SL = ED/SL) %>%
  mutate(Absolute_SL = SL/TL*CommonLength) %>%
  mutate(Absolute_ED = ED_SL * Absolute_SL)



#Diet categories


fishbase_info <-
  fishbase_info %>%
  mutate(Diet_cat = case_when(
    DietTroph <= 2.19  ~ "Herbivore",
    (DietTroph >= 2.2 & DietTroph <= 2.79) ~ "Omnivore",
    DietTroph >= 2.8 ~ "Carnivore"
  ))


fishbase_info <-
  fishbase_info %>%
  mutate(Food_cat = case_when(
    FoodTroph <= 2.19  ~ "Herbivore",
    (FoodTroph >= 2.2 & FoodTroph <= 2.79) ~ "Omnivore",
    FoodTroph >= 2.8 ~ "Carnivore"
  ))



#Simplify Anacat as migratory or non migratory.. 

AnaCat_vector <- fishbase_info$AnaCat
AnaCat_vector[AnaCat_vector==""]<-NA
AnaCat_vector[AnaCat_vector=="oceano-estuarine"]<-NA
fishbase_info$AnaCat <- AnaCat_vector



fishbase_info <- fishbase_info %>%
  mutate(Diadromous_cat = case_when(
    AnaCat %in%  c("non-migratory", "oceanodromous", "potamodromous") ~ "Non-diadromous",
    AnaCat %in% c("amphidromous", "anadromous", "catadromous") ~ "Diadromous",
    
  ))
  


#Extract columns which are of a particular interest


fishbase_info_filtered <- 
  fishbase_info %>%
  dplyr::select(my_species_name, Species, 
                ED_SL, Absolute_SL, Absolute_ED,
                CommonLength, Cave_Surface_filter, Fresh_salt,
                BodyShapeI, DemersPelag,
                Weight, Demersal,
                DietTroph, FoodTroph, EnvTemp, DietTroph,
                FoodTroph, Food_cat, Diet_cat, TempMin, TempMax,
                TempPreferred, Pelagic, Oceanic, Lakes, Epipelagic, 
                SubLittoral, LittoralZone, Neritic, Sand, Mud,
                Rocky, SoftBottom, HardBottom, AnaCat, Diadromous_cat)




fishbase_info_filtered[(fishbase_info_filtered$my_species_name == "Seriola lalandi dorsalis"),"DemersPelag"] <- "benthopelagic"
fishbase_info_filtered[(fishbase_info_filtered$my_species_name == "Oostethus manadensis"),"DemersPelag"] <- "demersal"
fishbase_info_filtered[(fishbase_info_filtered$my_species_name == "Heteromormyrus longilateralis"),"DemersPelag"] <- "benthopelagic"
fishbase_info_filtered[(fishbase_info_filtered$my_species_name == "Tenebrosternarchus preto"),"DemersPelag"] <- "pelagic"
fishbase_info_filtered[(fishbase_info_filtered$my_species_name == "Gymnocypris eckloni scoliostomus"),"DemersPelag"] <- "benthopelagic"
fishbase_info_filtered[(fishbase_info_filtered$my_species_name == "Rhinichthys klamathensis goyatoka"),"DemersPelag"] <- "demersal"



##### Data load - Depth -----------------

#Musilova and Coresti table #https://www.sciencedirect.com/science/article/pii/S0042698923000287
cortesi_depth_table <- 
  read.table("Depth_table_Cortesi_2023.tsv",
             sep="\t",
             header=TRUE)
cortesi_depth_table <- 
  cortesi_depth_table %>%
  dplyr::select(species, Depth_shallow, Depth_deep, Depth_mean)

colnames(cortesi_depth_table) <- c("species", "cortesi_depth_shallow", "cortesi_depth_deep",
                                   "cortesi_depth_mean")

#Duhamet table #https://onlinelibrary.wiley.com/doi/10.1002/ece3.9672
Duhamet_depth_table <- 
  read.table("Duhamet_2023.tsv",
             sep="\t",
             header=TRUE)

Duhamet_depth_table <- 
  Duhamet_depth_table %>%
  dplyr::select(Species, Depth_min, Depth_max)
colnames(Duhamet_depth_table) <- c("species", "duhamet_depth_shallow", "duhamet_depth_deep")

Duhamet_depth_table <- 
  as.data.frame(
    Duhamet_depth_table %>%
      rowwise() %>%
      mutate(duhamet_depth_mean = mean(c(duhamet_depth_shallow, duhamet_depth_deep), na.rm=TRUE))
  )

Duhamet_depth_table$species <- sub(" ", "_", Duhamet_depth_table$species)


#Fishbase table 
Fishbase_depth_table <- 
  fishbase_info %>%
  dplyr::select(my_species_name, 
                DepthRangeShallow,DepthRangeDeep)
colnames(Fishbase_depth_table) <- c("species", "fishbase_depth_shallow", "fishbase_depth_deep")

Fishbase_depth_table <- 
  as.data.frame(
    Fishbase_depth_table %>%
      rowwise() %>%
      mutate(fishbase_depth_mean = mean(c(fishbase_depth_shallow, fishbase_depth_deep), na.rm=TRUE))
  )


Fishbase_depth_table$species <- sub(" ", "_", Fishbase_depth_table$species)

#OBIS table

OBIS_depth_table <- 
  read.table("OBIS_depth_table.tsv",
             sep="\t",
             header=TRUE)

OBIS_depth_table <- 
  OBIS_depth_table %>%
  dplyr::select(species, mean_depth, mean_depth_min, mean_depth_max)
colnames(OBIS_depth_table) <- c("species", "OBIS_depth_mean","OBIS_depth_shallow", "OBIS_depth_deep")

OBIS_depth_table$species <- sub(" ", "_", OBIS_depth_table$species)



#Merge tables

Depth_Cortesi_Duhamet <- 
  full_join(cortesi_depth_table, Duhamet_depth_table, by="species")

Depth_Cortesi_Duhamet_Fishbase <- 
  full_join(Depth_Cortesi_Duhamet, Fishbase_depth_table, by="species")

Depth_Cortesi_Duhamet_Fishbase_OBIS <- 
  full_join(Depth_Cortesi_Duhamet_Fishbase, OBIS_depth_table, by="species")



Depth_Cortesi_Duhamet_Fishbase_OBIS[Depth_Cortesi_Duhamet_Fishbase_OBIS == "NA"] <- NA
Depth_Cortesi_Duhamet_Fishbase_OBIS[Depth_Cortesi_Duhamet_Fishbase_OBIS == "NaN"] <- NA


## DEEP FILTER
#Filter in this order of importance = Coresti > OBIS > Duhamet > Fishbase

cortesi_sp <- 
  Depth_Cortesi_Duhamet_Fishbase_OBIS %>%
  filter(! is.na(cortesi_depth_deep)) %>%
  pull(species)

OBIS_sp <- 
  Depth_Cortesi_Duhamet_Fishbase_OBIS %>%
  filter(is.na(cortesi_depth_deep)) %>%
  filter(! is.na(OBIS_depth_deep)) %>%
  pull(species)


duhamet_sp <- 
  Depth_Cortesi_Duhamet_Fishbase_OBIS %>%
  filter(is.na(cortesi_depth_deep)) %>%
  filter(is.na(OBIS_depth_deep)) %>%
  filter(! is.na(duhamet_depth_deep)) %>%
  pull(species)


fishbase_sp <- 
  Depth_Cortesi_Duhamet_Fishbase_OBIS %>%
  filter(is.na(cortesi_depth_deep)) %>%
  filter(is.na(OBIS_depth_deep)) %>%
  filter(is.na(duhamet_depth_deep)) %>%
  filter(! is.na(fishbase_depth_deep)) %>%
  pull(species)



first_filter <- 
  Depth_Cortesi_Duhamet_Fishbase_OBIS %>%
  mutate(depth_deep = case_when(
    species %in% cortesi_sp  ~ cortesi_depth_deep,
    species %in% OBIS_sp  ~ OBIS_depth_deep,
    species %in% duhamet_sp  ~ duhamet_depth_deep)) %>%
  mutate(depth_deep_source = case_when(
    species %in% cortesi_sp  ~ "Cortesi_2023",
    species %in% OBIS_sp  ~ "OBIS",
    species %in% duhamet_sp  ~ "Duhamet_2023")) %>%
  filter(! is.na(depth_deep)) %>%
  dplyr::select(species, depth_deep, depth_deep_source)


second_filter <- 
  Depth_Cortesi_Duhamet_Fishbase_OBIS %>%
  mutate(depth_deep = case_when(
    species %in% fishbase_sp  ~ fishbase_depth_deep)) %>%
  mutate(depth_deep_source = case_when(
    species %in% fishbase_sp  ~ "Fishbase")) %>%
  filter(! is.na(depth_deep)) %>%
  dplyr::select(species, depth_deep, depth_deep_source)


Deep_depth_df <- rbind(first_filter, second_filter)


## SHALLOW FILTER
#Filter in this order of importance = Coresti > OBIS > Duhamet > Fishbase

cortesi_sp <- 
  Depth_Cortesi_Duhamet_Fishbase_OBIS %>%
  filter(! is.na(cortesi_depth_shallow)) %>%
  pull(species)

OBIS_sp <- 
  Depth_Cortesi_Duhamet_Fishbase_OBIS %>%
  filter(is.na(cortesi_depth_shallow)) %>%
  filter(! is.na(OBIS_depth_shallow)) %>%
  pull(species)


duhamet_sp <- 
  Depth_Cortesi_Duhamet_Fishbase_OBIS %>%
  filter(is.na(cortesi_depth_shallow)) %>%
  filter(is.na(OBIS_depth_shallow)) %>%
  filter(! is.na(duhamet_depth_shallow)) %>%
  pull(species)


fishbase_sp <- 
  Depth_Cortesi_Duhamet_Fishbase_OBIS %>%
  filter(is.na(cortesi_depth_shallow)) %>%
  filter(is.na(OBIS_depth_shallow)) %>%
  filter(is.na(duhamet_depth_shallow)) %>%
  filter(! is.na(fishbase_depth_shallow)) %>%
  pull(species)




first_filter <- 
  Depth_Cortesi_Duhamet_Fishbase_OBIS %>%
  mutate(depth_shallow = case_when(
    species %in% cortesi_sp  ~ cortesi_depth_shallow,
    species %in% OBIS_sp  ~ OBIS_depth_shallow,
    species %in% duhamet_sp  ~ duhamet_depth_shallow)) %>%
  mutate(depth_shallow_source = case_when(
    species %in% cortesi_sp  ~ "Cortesi_2023",
    species %in% OBIS_sp  ~ "OBIS",
    species %in% duhamet_sp  ~ "Duhamet_2023")) %>%
  filter(! is.na(depth_shallow)) %>%
  dplyr::select(species, depth_shallow, depth_shallow_source)


second_filter <- 
  Depth_Cortesi_Duhamet_Fishbase_OBIS %>%
  mutate(depth_shallow = case_when(
    species %in% fishbase_sp  ~ fishbase_depth_shallow)) %>%
  mutate(depth_shallow_source = case_when(
    species %in% fishbase_sp  ~ "Fishbase")) %>%
  filter(! is.na(depth_shallow)) %>%
  dplyr::select(species, depth_shallow, depth_shallow_source)


Shallow_depth_df <- rbind(first_filter, second_filter)


## MEAN FILTER
#Filter in this order of importance = Coresti > OBIS > Duhamet > Fishbase



cortesi_sp <- 
  Depth_Cortesi_Duhamet_Fishbase_OBIS %>%
  filter(! is.na(cortesi_depth_mean)) %>%
  pull(species)

OBIS_sp <- 
  Depth_Cortesi_Duhamet_Fishbase_OBIS %>%
  filter(is.na(cortesi_depth_mean)) %>%
  filter(! is.na(OBIS_depth_mean)) %>%
  pull(species)


duhamet_sp <- 
  Depth_Cortesi_Duhamet_Fishbase_OBIS %>%
  filter(is.na(cortesi_depth_mean)) %>%
  filter(is.na(OBIS_depth_mean)) %>%
  filter(! is.na(duhamet_depth_mean)) %>%
  pull(species)


fishbase_sp <- 
  Depth_Cortesi_Duhamet_Fishbase_OBIS %>%
  filter(is.na(cortesi_depth_mean)) %>%
  filter(is.na(OBIS_depth_mean)) %>%
  filter(is.na(duhamet_depth_mean)) %>%
  filter(! is.na(fishbase_depth_mean)) %>%
  pull(species)



first_filter <- 
  Depth_Cortesi_Duhamet_Fishbase_OBIS %>%
  mutate(depth_mean = case_when(
    species %in% cortesi_sp  ~ cortesi_depth_mean,
    species %in% OBIS_sp  ~ OBIS_depth_mean,
    species %in% duhamet_sp  ~ duhamet_depth_mean)) %>%
  mutate(depth_mean_source = case_when(
    species %in% cortesi_sp  ~ "Cortesi_2023",
    species %in% OBIS_sp  ~ "OBIS",
    species %in% duhamet_sp  ~ "Duhamet_2023")) %>%
  filter(! is.na(depth_mean)) %>%
  dplyr::select(species, depth_mean, depth_mean_source)


second_filter <- 
  Depth_Cortesi_Duhamet_Fishbase_OBIS %>%
  mutate(depth_mean = case_when(
    species %in% fishbase_sp  ~ fishbase_depth_mean)) %>%
  mutate(depth_mean_source = case_when(
    species %in% fishbase_sp  ~ "Fishbase")) %>%
  filter(! is.na(depth_mean)) %>%
  dplyr::select(species, depth_mean, depth_mean_source)


Mean_depth_df <- rbind(first_filter, second_filter)

## Merge dataframes


mean_shallow <- full_join(Mean_depth_df, Shallow_depth_df, by="species")
Depth_dataframe <- full_join(mean_shallow, Deep_depth_df, by="species")



species_df_filter <- 
  species_table %>%
  filter(genome_name %in% final_dataset) %>%
  dplyr::select(genome_name)

colnames(species_df_filter) <- "species"



Depth_dataframe_filter <- 
  left_join(species_df_filter, Depth_dataframe, by="species")


#Set Riverine/pond/spring living species to 2 m (Following Musilova and Cortesi 2023)


River_sp <- scan("River_sp.txt", what="character")


lacking_sp <- 
  Depth_dataframe_filter %>%
  filter(is.na(depth_mean)) %>%
  pull(species) 


river_to_add <- 
  intersect(River_sp, lacking_sp)



for(curr_species in river_to_add){
  
  Depth_dataframe_filter[(Depth_dataframe_filter$species ==  curr_species),"depth_mean"] <- 2
  Depth_dataframe_filter[(Depth_dataframe_filter$species ==  curr_species),"depth_shallow"] <- 2
  Depth_dataframe_filter[(Depth_dataframe_filter$species ==  curr_species),"depth_deep"] <- 2
  
  Depth_dataframe_filter[(Depth_dataframe_filter$species ==  curr_species),"depth_shallow_source"] <- "River"
  Depth_dataframe_filter[(Depth_dataframe_filter$species ==  curr_species),"depth_mean_source"] <- "River"
  Depth_dataframe_filter[(Depth_dataframe_filter$species ==  curr_species),"depth_deep_source"] <- "River"
  
}


nrow(Depth_dataframe_filter %>%
       filter(is.na(depth_mean)))





##### Data load - Opsins table  ---------------------------------

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



## same with pseudogenes


visual_opsins_table_count_pseudogene <- 
  visual_opsins_table_count %>%
  filter(gene_state == "Pseudogene")

non_visual_opsins_table_count_pseudogene <- 
  non_visual_opsins_table_count %>%
  filter(gene_state == "Pseudogene")


all_opsins_table_count_pseudogene <- 
  rbind(visual_opsins_table_count_pseudogene,
        non_visual_opsins_table_count_pseudogene)

visual_opsins_table_count_pseudogene_wide <- 
  as.data.frame(
    visual_opsins_table_count_pseudogene %>%
      dplyr::select(! gene_state) %>%
      pivot_wider(names_from = clade, values_from = count)
  )

visual_opsins_table_count_pseudogene_wide <- 
  replace(visual_opsins_table_count_pseudogene_wide, is.na(visual_opsins_table_count_pseudogene_wide), 0)



non_visual_opsins_table_count_pseudogene_wide <- 
  as.data.frame(
    non_visual_opsins_table_count_pseudogene %>%
      dplyr::select(! gene_state) %>%
      pivot_wider(names_from = clade, values_from = count)
  )

non_visual_opsins_table_count_pseudogene_wide <- 
  replace(non_visual_opsins_table_count_pseudogene_wide, 
          is.na(non_visual_opsins_table_count_pseudogene_wide), 0)


visual_opsins_table_count_pseudogene_wide <- 
  visual_opsins_table_count_pseudogene_wide %>%
  mutate(total_visual = rowSums(across(where(is.numeric))))


non_visual_opsins_table_count_pseudogene_wide <- 
  non_visual_opsins_table_count_pseudogene_wide %>%
  mutate(total_non_visual = rowSums(across(where(is.numeric))))

all_ospins_pseudogene_wide <- 
  full_join(non_visual_opsins_table_count_pseudogene_wide, 
            visual_opsins_table_count_pseudogene_wide, by="species")

all_ospins_pseudogene_wide <- 
  replace(all_ospins_pseudogene_wide, is.na(all_ospins_pseudogene_wide), 0)


all_ospins_pseudogene_wide <- 
  all_ospins_pseudogene_wide %>%
  mutate(total_opsins = 
           total_visual + total_non_visual)


sp_list_df <- as.data.frame(visual_opsins_table_count %>% pull(species) %>% unique())
colnames(sp_list_df) <- c("species")



all_ospins_pseudogene_wide <- left_join(sp_list_df, all_ospins_pseudogene_wide, by='species')
all_ospins_pseudogene_wide <- 
  replace(all_ospins_pseudogene_wide, is.na(all_ospins_pseudogene_wide), 0)



all_ospins_pseudogene_wide %>% pull(total_opsins) %>% sum()





## same with incomplete


visual_opsins_table_count_incomplete <- 
  visual_opsins_table_count %>%
  filter(gene_state == "Incomplete")

non_visual_opsins_table_count_incomplete <- 
  non_visual_opsins_table_count %>%
  filter(gene_state == "Incomplete")


all_opsins_table_count_incomplete <- 
  rbind(visual_opsins_table_count_incomplete,
        non_visual_opsins_table_count_incomplete)

visual_opsins_table_count_incomplete_wide <- 
  as.data.frame(
    visual_opsins_table_count_incomplete %>%
      dplyr::select(! gene_state) %>%
      pivot_wider(names_from = clade, values_from = count)
  )

visual_opsins_table_count_incomplete_wide <- 
  replace(visual_opsins_table_count_incomplete_wide, is.na(visual_opsins_table_count_incomplete_wide), 0)



non_visual_opsins_table_count_incomplete_wide <- 
  as.data.frame(
    non_visual_opsins_table_count_incomplete %>%
      dplyr::select(! gene_state) %>%
      pivot_wider(names_from = clade, values_from = count)
  )

non_visual_opsins_table_count_incomplete_wide <- 
  replace(non_visual_opsins_table_count_incomplete_wide, 
          is.na(non_visual_opsins_table_count_incomplete_wide), 0)


visual_opsins_table_count_incomplete_wide <- 
  visual_opsins_table_count_incomplete_wide %>%
  mutate(total_visual = rowSums(across(where(is.numeric))))


non_visual_opsins_table_count_incomplete_wide <- 
  non_visual_opsins_table_count_incomplete_wide %>%
  mutate(total_non_visual = rowSums(across(where(is.numeric))))

all_ospins_incomplete_wide <- 
  full_join(non_visual_opsins_table_count_incomplete_wide, 
            visual_opsins_table_count_incomplete_wide, by="species")

all_ospins_incomplete_wide <- 
  replace(all_ospins_incomplete_wide, is.na(all_ospins_incomplete_wide), 0)


all_ospins_incomplete_wide <- 
  all_ospins_incomplete_wide %>%
  mutate(total_opsins = 
           total_visual + total_non_visual)


sp_list_df <- as.data.frame(visual_opsins_table_count %>% pull(species) %>% unique())
colnames(sp_list_df) <- c("species")



all_ospins_incomplete_wide <- left_join(sp_list_df, all_ospins_incomplete_wide, by='species')
all_ospins_incomplete_wide <- 
  replace(all_ospins_incomplete_wide, is.na(all_ospins_incomplete_wide), 0)



all_ospins_incomplete_wide %>% pull(total_opsins) %>% sum()


##### Data load - Nocturnality   ---------------------------------

Diur_Noctu_table <-
  read.table("Nocturnal_Diurnal_DB_woref.csv",
             header=TRUE,
             sep=",", 
             quote = "")

colnames(Diur_Noctu_table) <- c("species", "Diel_pattern")




Diur_Noctu_table <-
  Diur_Noctu_table %>%
  mutate(Diel_pattern_simplified = case_when(
    Diel_pattern == "Cathemeral/Nocturnal" ~ "Nocturnal",
    Diel_pattern == "Nocturnal" ~ "Nocturnal",
    Diel_pattern == "Crepuscular/Nocturnal" ~ "Nocturnal",
    Diel_pattern == "Diurnal" ~ "Diurnal",
    Diel_pattern == "Crepuscular/Diurnal" ~ "Diurnal",
    Diel_pattern == "Crepuscular" ~ "Crepuscular",
    Diel_pattern == "Cathemeral" ~ "Cathemeral"
  ))



Diur_Noctu_table$species <- sub(" ", "_", Diur_Noctu_table$species)



##### Data load - Electric   ---------------------------------

electric_table <-
  read.table("Electrogenic_table.tsv",
             header=TRUE,
             sep="\t")

##### Data load - Amphibious   ---------------------------------


amphibious_table <-
  read.table("Amphibious_table.tsv",
             header=TRUE,
             sep="\t")


curr_vector <- amphibious_table$Amphibious
curr_vector[curr_vector == 1]<- "Yes"
curr_vector[curr_vector == 0 ]<- "No"
amphibious_table$Amphibious <- curr_vector

curr_vector <- amphibious_table$Terrestrial_exploration
curr_vector[curr_vector == 1]<- "Yes"
curr_vector[curr_vector == 0 ]<- "No"
amphibious_table$Terrestrial_exploration <- curr_vector


##### Data load - Chemoreceptor genes --------------------

chemoreceptor_table <-
  read.table("Chemoreceptor_table_policarpo.tsv",
             header=TRUE,
             sep="\t")

##### Data load - Genome sizes --------------------


genome_size_table <- 
  read.table("Genomes_sizes.tsv",
             header=FALSE,
             sep=",")

colnames(genome_size_table) <- c("species", "genome_size_bp")

genome_size_table <- 
  genome_size_table %>%
  mutate(genome_size_mbp = genome_size_bp/1000000)

##### Data load - OGG --------------------

OGG_table <- 
  read.table("Orthogroups.GeneCount.tsv",
             header=TRUE,
             sep="\t")

#Turn the df

rownames(OGG_table)<- OGG_table$Orthogroup
OGG_table$Orthogroup<-NULL
OGG_table_df <- as.data.frame(t(OGG_table))
OGG_table_df <- tibble::rownames_to_column(OGG_table_df, "species")


#Import chemoreceptors and opsins OGG 

opsins_OGG <- scan("list_opsin_OGG.txt",
                   what="character")

OLR_OGG <- scan("list_OLR_OGG.txt",
                   what="character")

#Remove chemoreceptor and opsins OGG from the OGG table

OGG_table_df_filtered <- 
  OGG_table_df %>%
  dplyr::select(! c(opsins_OGG))


##### Data load - Species trees   ---------------------------------


ASTRAL_species_tree_nodelabel <- 
  read.tree("AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel")
ASTRAL_species_tree_nodelabel$edge.length <- ASTRAL_species_tree_nodelabel$edge.length * 1000



##### Order tree + Circle for each subfamily  ---------------------------------


colnames(species_table) <- c("species", colnames(species_table)[c(2:length(colnames(species_table)))])

all_ospins_complete_wide_final_dataset <- 
  all_ospins_complete_wide %>%
  filter(species %in% final_dataset)

all_ospins_complete_wide_final_dataset <- 
  left_join(all_ospins_complete_wide_final_dataset, species_table, by="species")

order_list <- 
  species_table %>%
  filter(species %in% final_dataset) %>%
  pull(order) %>%
  unique()

order_tree <- ASTRAL_species_tree_nodelabel
all_species_list <- order_tree$tip.label


all_ospins_complete_wide_final_dataset %>%
  filter(order == "Batrachoidiformes")




Cyprin_main_WGD <- 
  all_ospins_complete_wide_final_dataset %>% 
  filter(family == "Cyprinidae") %>% 
  filter(WGD == "Yes") %>%
  filter(! species %in% c("Barbus_barbus", "Oxygymnocypris_stewartii")) %>%
  pull(species)


Cyprin_no_WGD <- 
  all_ospins_complete_wide_final_dataset %>% 
  filter(order == "Cypriniformes") %>% 
  filter(WGD == "No") %>%
  pull(species)

order_list <- order_list[! order_list %in% c('Acipenseriformes', 'Cypriniformes')]

for(curr_order in order_list){
  
  curr_species_list <- 
    species_table %>%
    filter(species %in% final_dataset) %>%
    filter(order == curr_order) %>%
    pull(species)
  
  curr_species <- curr_species_list[1]
  
  all_species_list[all_species_list == curr_species] <- curr_order
  
  
}

all_species_list[all_species_list == "Xyrauchen_texanus"] <- "Catostomidae"
all_species_list[all_species_list == "Carassius_auratus"] <- "Cypriniforme_and_WGD"
all_species_list[all_species_list == "Danio_rerio"] <- "Cypriniforme_no_WGD"



order_tree$tip.label <- all_species_list

order_list <- c(order_list, "Catostomidae",
                "Oxygymnocypris_stewartii", "Barbus_barbus",
                "Cypriniforme_no_WGD", "Polyodon_spathula", "Acipenser_ruthenus", "Cypriniforme_and_WGD")


order_tree <- 
  keep.tip(order_tree, 
           c(order_list))



ggtree(order_tree) + 
  xlim(0, 600) +
  geom_tiplab()


#Prepare species/opsin table

all_ospins_complete_wide_final_dataset <- 
  all_ospins_complete_wide %>%
  filter(species %in% final_dataset)

all_ospins_complete_wide_final_dataset <- 
  left_join(all_ospins_complete_wide_final_dataset, species_table, by="species")





all_ospins_complete_wide_final_dataset <- 
  all_ospins_complete_wide_final_dataset %>%
  mutate(order_reprez_modif = case_when(
    species == "Acipenser_ruthenus" ~ "Acipenser_ruthenus",
    species == "Polyodon_spathula" ~ "Polyodon_spathula",
    species %in% c("Myxocyprinus_asiaticus", "Xyrauchen_texanus") ~ "Catostomidae",
    species == "Barbus_barbus" ~ "Barbus_barbus",
    species == "Oxygymnocypris_stewartii" ~ "Oxygymnocypris_stewartii",
    species %in% Cyprin_main_WGD ~ "Cypriniforme_and_WGD",
    species %in% Cyprin_no_WGD ~ "Cypriniforme_no_WGD"
  ))


all_ospins_complete_wide_final_dataset <- 
  all_ospins_complete_wide_final_dataset %>%
  mutate(order_reprez = if_else(
    is.na(order_reprez_modif),
    order,
    order_reprez_modif
  ))



all_ospins_complete_wide_final_dataset %>%
  filter(infraclass != "Teleostei")
#Prepare tables of presence absence and mean for each order ... 


colnames(all_ospins_complete_wide_final_dataset) <- 
  sub("-", "_", colnames(all_ospins_complete_wide_final_dataset))

exorh_prez_abs <- 
  as.data.frame(
    all_ospins_complete_wide_final_dataset %>%
      dplyr::select(order_reprez, exorh) %>%
      group_by(order_reprez) %>%
      summarise(mean_comp = mean(exorh)) %>%
      mutate(prez_abs = if_else(
        mean_comp > 0,
        "Presence", 
        "Absence"
      )) %>%
      dplyr::select(order_reprez, prez_abs, mean_comp)
  )

lws_prez_abs <- 
  as.data.frame(
    all_ospins_complete_wide_final_dataset %>%
      dplyr::select(order_reprez, lws) %>%
      group_by(order_reprez) %>%
      summarise(mean_comp = mean(lws)) %>%
      mutate(prez_abs = if_else(
        mean_comp > 0,
        "Presence", 
        "Absence"
      )) %>%
      dplyr::select(order_reprez, prez_abs, mean_comp)
  )
opn3_prez_abs <- 
  as.data.frame(
    all_ospins_complete_wide_final_dataset %>%
      dplyr::select(order_reprez, opn3) %>%
      group_by(order_reprez) %>%
      summarise(mean_comp = mean(opn3)) %>%
      mutate(prez_abs = if_else(
        mean_comp > 0,
        "Presence", 
        "Absence"
      )) %>%
      dplyr::select(order_reprez, prez_abs, mean_comp)
  )
opn4m1_3_prez_abs <- 
  as.data.frame(
    all_ospins_complete_wide_final_dataset %>%
      dplyr::select(order_reprez, opn4m1_3) %>%
      group_by(order_reprez) %>%
      summarise(mean_comp = mean(opn4m1_3)) %>%
      mutate(prez_abs = if_else(
        mean_comp > 0,
        "Presence", 
        "Absence"
      )) %>%
      dplyr::select(order_reprez, prez_abs, mean_comp)
  )
opn4m2_prez_abs <- 
  as.data.frame(
    all_ospins_complete_wide_final_dataset %>%
      dplyr::select(order_reprez, opn4m2) %>%
      group_by(order_reprez) %>%
      summarise(mean_comp = mean(opn4m2)) %>%
      mutate(prez_abs = if_else(
        mean_comp > 0,
        "Presence", 
        "Absence"
      )) %>%
      dplyr::select(order_reprez, prez_abs, mean_comp)
  )
opn4x_prez_abs <- 
  as.data.frame(
    all_ospins_complete_wide_final_dataset %>%
      dplyr::select(order_reprez, opn4x) %>%
      group_by(order_reprez) %>%
      summarise(mean_comp = mean(opn4x)) %>%
      mutate(prez_abs = if_else(
        mean_comp > 0,
        "Presence", 
        "Absence"
      )) %>%
      dplyr::select(order_reprez, prez_abs, mean_comp)
  )
opn5_prez_abs <- 
  as.data.frame(
    all_ospins_complete_wide_final_dataset %>%
      dplyr::select(order_reprez, opn5) %>%
      group_by(order_reprez) %>%
      summarise(mean_comp = mean(opn5)) %>%
      mutate(prez_abs = if_else(
        mean_comp > 0,
        "Presence", 
        "Absence"
      )) %>%
      dplyr::select(order_reprez, prez_abs, mean_comp)
  )
opn6_prez_abs <- 
  as.data.frame(
    all_ospins_complete_wide_final_dataset %>%
      dplyr::select(order_reprez, opn6) %>%
      group_by(order_reprez) %>%
      summarise(mean_comp = mean(opn6)) %>%
      mutate(prez_abs = if_else(
        mean_comp > 0,
        "Presence", 
        "Absence"
      )) %>%
      dplyr::select(order_reprez, prez_abs, mean_comp)
  )
opn7a_prez_abs <- 
  as.data.frame(
    all_ospins_complete_wide_final_dataset %>%
      dplyr::select(order_reprez, opn7a) %>%
      group_by(order_reprez) %>%
      summarise(mean_comp = mean(opn7a)) %>%
      mutate(prez_abs = if_else(
        mean_comp > 0,
        "Presence", 
        "Absence"
      )) %>%
      dplyr::select(order_reprez, prez_abs, mean_comp)
  )
opn7b_prez_abs <- 
  as.data.frame(
    all_ospins_complete_wide_final_dataset %>%
      dplyr::select(order_reprez, opn7b) %>%
      group_by(order_reprez) %>%
      summarise(mean_comp = mean(opn7b)) %>%
      mutate(prez_abs = if_else(
        mean_comp > 0,
        "Presence", 
        "Absence"
      )) %>%
      dplyr::select(order_reprez, prez_abs, mean_comp)
  )
opn8a_prez_abs <- 
  as.data.frame(
    all_ospins_complete_wide_final_dataset %>%
      dplyr::select(order_reprez, opn8a) %>%
      group_by(order_reprez) %>%
      summarise(mean_comp = mean(opn8a)) %>%
      mutate(prez_abs = if_else(
        mean_comp > 0,
        "Presence", 
        "Absence"
      )) %>%
      dplyr::select(order_reprez, prez_abs, mean_comp)
  )
opn8b_prez_abs <- 
  as.data.frame(
    all_ospins_complete_wide_final_dataset %>%
      dplyr::select(order_reprez, opn8b) %>%
      group_by(order_reprez) %>%
      summarise(mean_comp = mean(opn8b)) %>%
      mutate(prez_abs = if_else(
        mean_comp > 0,
        "Presence", 
        "Absence"
      )) %>%
      dplyr::select(order_reprez, prez_abs, mean_comp)
  )
opn8c_prez_abs <- 
  as.data.frame(
    all_ospins_complete_wide_final_dataset %>%
      dplyr::select(order_reprez, opn8c) %>%
      group_by(order_reprez) %>%
      summarise(mean_comp = mean(opn8c)) %>%
      mutate(prez_abs = if_else(
        mean_comp > 0,
        "Presence", 
        "Absence"
      )) %>%
      dplyr::select(order_reprez, prez_abs, mean_comp)
  )
opn9_prez_abs <- 
  as.data.frame(
    all_ospins_complete_wide_final_dataset %>%
      dplyr::select(order_reprez, opn9) %>%
      group_by(order_reprez) %>%
      summarise(mean_comp = mean(opn9)) %>%
      mutate(prez_abs = if_else(
        mean_comp > 0,
        "Presence", 
        "Absence"
      )) %>%
      dplyr::select(order_reprez, prez_abs, mean_comp)
  )
parapinopsin_prez_abs <- 
  as.data.frame(
    all_ospins_complete_wide_final_dataset %>%
      dplyr::select(order_reprez, parapinopsin) %>%
      group_by(order_reprez) %>%
      summarise(mean_comp = mean(parapinopsin)) %>%
      mutate(prez_abs = if_else(
        mean_comp > 0,
        "Presence", 
        "Absence"
      )) %>%
      dplyr::select(order_reprez, prez_abs, mean_comp)
  )
parietopsin_prez_abs <- 
  as.data.frame(
    all_ospins_complete_wide_final_dataset %>%
      dplyr::select(order_reprez, parietopsin) %>%
      group_by(order_reprez) %>%
      summarise(mean_comp = mean(parietopsin)) %>%
      mutate(prez_abs = if_else(
        mean_comp > 0,
        "Presence", 
        "Absence"
      )) %>%
      dplyr::select(order_reprez, prez_abs, mean_comp)
  )
pinopsin_prez_abs <- 
  as.data.frame(
    all_ospins_complete_wide_final_dataset %>%
      dplyr::select(order_reprez, pinopsin) %>%
      group_by(order_reprez) %>%
      summarise(mean_comp = mean(pinopsin)) %>%
      mutate(prez_abs = if_else(
        mean_comp > 0,
        "Presence", 
        "Absence"
      )) %>%
      dplyr::select(order_reprez, prez_abs, mean_comp)
  )
rgr_prez_abs <- 
  as.data.frame(
    all_ospins_complete_wide_final_dataset %>%
      dplyr::select(order_reprez, rgr) %>%
      group_by(order_reprez) %>%
      summarise(mean_comp = mean(rgr)) %>%
      mutate(prez_abs = if_else(
        mean_comp > 0,
        "Presence", 
        "Absence"
      )) %>%
      dplyr::select(order_reprez, prez_abs, mean_comp)
  )
rh1_prez_abs <- 
  as.data.frame(
    all_ospins_complete_wide_final_dataset %>%
      dplyr::select(order_reprez, rh1) %>%
      group_by(order_reprez) %>%
      summarise(mean_comp = mean(rh1)) %>%
      mutate(prez_abs = if_else(
        mean_comp > 0,
        "Presence", 
        "Absence"
      )) %>%
      dplyr::select(order_reprez, prez_abs, mean_comp)
  )
rh2_prez_abs <- 
  as.data.frame(
    all_ospins_complete_wide_final_dataset %>%
      dplyr::select(order_reprez, rh2) %>%
      group_by(order_reprez) %>%
      summarise(mean_comp = mean(rh2)) %>%
      mutate(prez_abs = if_else(
        mean_comp > 0,
        "Presence", 
        "Absence"
      )) %>%
      dplyr::select(order_reprez, prez_abs, mean_comp)
  )
rrh_prez_abs <- 
  as.data.frame(
    all_ospins_complete_wide_final_dataset %>%
      dplyr::select(order_reprez, rrh) %>%
      group_by(order_reprez) %>%
      summarise(mean_comp = mean(rrh)) %>%
      mutate(prez_abs = if_else(
        mean_comp > 0,
        "Presence", 
        "Absence"
      )) %>%
      dplyr::select(order_reprez, prez_abs, mean_comp)
  )
sws1_prez_abs <- 
  as.data.frame(
    all_ospins_complete_wide_final_dataset %>%
      dplyr::select(order_reprez, sws1) %>%
      group_by(order_reprez) %>%
      summarise(mean_comp = mean(sws1)) %>%
      mutate(prez_abs = if_else(
        mean_comp > 0,
        "Presence", 
        "Absence"
      )) %>%
      dplyr::select(order_reprez, prez_abs, mean_comp)
  )
sws2_prez_abs <- 
  as.data.frame(
    all_ospins_complete_wide_final_dataset %>%
      dplyr::select(order_reprez, sws2) %>%
      group_by(order_reprez) %>%
      summarise(mean_comp = mean(sws2)) %>%
      mutate(prez_abs = if_else(
        mean_comp > 0,
        "Presence", 
        "Absence"
      )) %>%
      dplyr::select(order_reprez, prez_abs, mean_comp)
  )
tmt1_prez_abs <- 
  as.data.frame(
    all_ospins_complete_wide_final_dataset %>%
      dplyr::select(order_reprez, tmt1) %>%
      group_by(order_reprez) %>%
      summarise(mean_comp = mean(tmt1)) %>%
      mutate(prez_abs = if_else(
        mean_comp > 0,
        "Presence", 
        "Absence"
      )) %>%
      dplyr::select(order_reprez, prez_abs, mean_comp)
  )
tmt2_prez_abs <- 
  as.data.frame(
    all_ospins_complete_wide_final_dataset %>%
      dplyr::select(order_reprez, tmt2) %>%
      group_by(order_reprez) %>%
      summarise(mean_comp = mean(tmt2)) %>%
      mutate(prez_abs = if_else(
        mean_comp > 0,
        "Presence", 
        "Absence"
      )) %>%
      dplyr::select(order_reprez, prez_abs, mean_comp)
  )
tmt3_prez_abs <- 
  as.data.frame(
    all_ospins_complete_wide_final_dataset %>%
      dplyr::select(order_reprez, tmt3) %>%
      group_by(order_reprez) %>%
      summarise(mean_comp = mean(tmt3)) %>%
      mutate(prez_abs = if_else(
        mean_comp > 0,
        "Presence", 
        "Absence"
      )) %>%
      dplyr::select(order_reprez, prez_abs, mean_comp)
  )
va_prez_abs <- 
  as.data.frame(
    all_ospins_complete_wide_final_dataset %>%
      dplyr::select(order_reprez, va) %>%
      group_by(order_reprez) %>%
      summarise(mean_comp = mean(va)) %>%
      mutate(prez_abs = if_else(
        mean_comp > 0,
        "Presence", 
        "Absence"
      )) %>%
      dplyr::select(order_reprez, prez_abs, mean_comp)
  )




ggtree(order_tree, size=2,
       branch.length="none") +
  geom_tiplab() +
  xlim(0, 100)


p1 <- ggtree(order_tree, size=2,
             branch.length="none") +
  geom_tiplab() + 
  geom_facet(
    panel="rh1",
    data=rh1_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="black",
    width = 0.5
  ) +
  geom_facet(
    panel="exorh",
    data=exorh_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#E6F5D0",
    width = 0.5
  ) +
  geom_facet(
    panel="rh2",
    data=rh2_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#276419",
    width = 0.5
  ) +
  geom_facet(
    panel="sws2",
    data=sws2_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#41B6C4",
    width = 0.5
  ) +
  geom_facet(
    panel="sws1",
    data=sws1_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#FEE391",
    width = 0.5
  ) +
  geom_facet(
    panel="lws",
    data=lws_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#67001F",
    width = 0.5
  ) +
  geom_facet(
    panel="pinopsin",
    data=pinopsin_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#3690C0",
    width = 0.5
  ) +
  geom_facet(
    panel="va",
    data=va_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#8C510A",
    width = 0.5
  ) +
  geom_facet(
    panel="parapinopsin",
    data=parapinopsin_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#A8DDB5",
    width = 0.5
  ) +
  geom_facet(
    panel="parietopsin",
    data=parietopsin_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#FDAE61",
    width = 0.5
  ) +
  geom_facet(
    panel="tmt1",
    data=tmt1_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#FC4E2A",
    width = 0.5
  ) +
  geom_facet(
    panel="tmt2",
    data=tmt2_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#045A8D",
    width = 0.5
  ) +
  geom_facet(
    panel="tmt3",
    data=tmt3_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#D1E5F0",
    width = 0.5
  ) +
  geom_facet(
    panel="opn3",
    data=opn3_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#F4A582",
    width = 0.5
  ) +
  geom_facet(
    panel="opn8a",
    data=opn8a_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#B3B3B3",
    width = 0.5
  ) +
  geom_facet(
    panel="opn8b",
    data=opn8b_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#2D004B",
    width = 0.5
  ) +
  geom_facet(
    panel="opn8c",
    data=opn8c_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="goldenrod",
    width = 0.5
  ) +
  geom_facet(
    panel="opn7a",
    data=opn7a_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#08519C",
    width = 0.5
  ) +
  geom_facet(
    panel="opn7b",
    data=opn7b_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#EF6548",
    width = 0.5
  ) +
  geom_facet(
    panel="opn5",
    data=opn5_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#9970AB",
    width = 0.5
  ) +
  geom_facet(
    panel="opn9",
    data=opn9_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#88419D",
    width = 0.5
  ) +
  geom_facet(
    panel="opn6",
    data=opn6_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#B15928",
    width = 0.5
  ) +
  geom_facet(
    panel="rgr",
    data=rgr_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#C2A5CF",
    width = 0.5
  ) +
  geom_facet(
    panel="rrh",
    data=rrh_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#ABD9E9",
    width = 0.5
  ) +
  geom_facet(
    panel="opn4m1_3",
    data=opn4m1_3_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#35978F",
    width = 0.5
  ) +
  geom_facet(
    panel="opn4m2",
    data=opn4m2_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#238443",
    width = 0.5
  ) +
  geom_facet(
    panel="opn4x",
    data=opn4x_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#525252",
    width = 0.5
  ) +
  scale_alpha_manual(values = c(
    Presence=1,
    Absence =0
  )) +
  scale_size(range = c(2, 10)) +
  theme(legend.position = "none") 




gt = ggplot_gtable(ggplot_build(p1))
gt$layout$l[grep('panel-1', gt$layout$name)] 
gt$widths[5] = 6*gt$widths[5]
gt$widths[7] = 0.5*gt$widths[7]
gt$widths[9] = 0.5*gt$widths[9]
gt$widths[11] = 0.5*gt$widths[11]
gt$widths[13] = 0.5*gt$widths[13]
gt$widths[15] = 0.5*gt$widths[15]
gt$widths[17] = 0.5*gt$widths[17]
gt$widths[19] = 0.5*gt$widths[19]
gt$widths[21] = 0.5*gt$widths[21]
gt$widths[23] = 0.5*gt$widths[23]
gt$widths[25] = 0.5*gt$widths[25]
gt$widths[27] = 0.5*gt$widths[27]
gt$widths[29] = 0.5*gt$widths[29]
gt$widths[31] = 0.5*gt$widths[31]
gt$widths[33] = 0.5*gt$widths[33]
gt$widths[35] = 0.5*gt$widths[35]
gt$widths[37] = 0.5*gt$widths[37]
gt$widths[39] = 0.5*gt$widths[39]
gt$widths[41] = 0.5*gt$widths[41]
gt$widths[43] = 0.5*gt$widths[43]
gt$widths[45] = 0.5*gt$widths[45]
gt$widths[47] = 0.5*gt$widths[47]
gt$widths[49] = 0.5*gt$widths[49]
gt$widths[51] = 0.5*gt$widths[51]
gt$widths[53] = 0.5*gt$widths[53]
gt$widths[55] = 0.5*gt$widths[55]
gt$widths[57] = 0.5*gt$widths[57]
gt$widths[59] = 0.5*gt$widths[59]


grid.draw(gt) # plot with grid draw


all_ospins_complete_wide_final_dataset %>%
  filter(infraclass != "Teleostei")

ggtree(order_tree, size=2,
       branch.length="none") +
  geom_tiplab() + 
  geom_facet(
    panel="rh1",
    data=rh1_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="black",
    width = 0.5
  ) +
  geom_facet(
    panel="exorh",
    data=exorh_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#E6F5D0",
    width = 0.5
  ) +
  geom_facet(
    panel="rh2",
    data=rh2_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#276419",
    width = 0.5
  ) +
  geom_facet(
    panel="sws2",
    data=sws2_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#41B6C4",
    width = 0.5
  ) +
  geom_facet(
    panel="sws1",
    data=sws1_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#FEE391",
    width = 0.5
  ) +
  geom_facet(
    panel="lws",
    data=lws_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#67001F",
    width = 0.5
  ) +
  geom_facet(
    panel="pinopsin",
    data=pinopsin_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#3690C0",
    width = 0.5
  ) +
  geom_facet(
    panel="va",
    data=va_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#8C510A",
    width = 0.5
  ) +
  geom_facet(
    panel="parapinopsin",
    data=parapinopsin_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#A8DDB5",
    width = 0.5
  ) +
  geom_facet(
    panel="parietopsin",
    data=parietopsin_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#FDAE61",
    width = 0.5
  ) +
  geom_facet(
    panel="tmt1",
    data=tmt1_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#FC4E2A",
    width = 0.5
  ) +
  geom_facet(
    panel="tmt2",
    data=tmt2_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#045A8D",
    width = 0.5
  ) +
  geom_facet(
    panel="tmt3",
    data=tmt3_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#D1E5F0",
    width = 0.5
  ) +
  geom_facet(
    panel="opn3",
    data=opn3_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#F4A582",
    width = 0.5
  ) +
  geom_facet(
    panel="opn8a",
    data=opn8a_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#B3B3B3",
    width = 0.5
  ) +
  geom_facet(
    panel="opn8b",
    data=opn8b_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#2D004B",
    width = 0.5
  ) +
  geom_facet(
    panel="opn8c",
    data=opn8c_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="goldenrod",
    width = 0.5
  ) +
  geom_facet(
    panel="opn7a",
    data=opn7a_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#08519C",
    width = 0.5
  ) +
  geom_facet(
    panel="opn7b",
    data=opn7b_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#EF6548",
    width = 0.5
  ) +
  geom_facet(
    panel="opn5",
    data=opn5_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#9970AB",
    width = 0.5
  ) +
  geom_facet(
    panel="opn9",
    data=opn9_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#88419D",
    width = 0.5
  ) +
  geom_facet(
    panel="opn6",
    data=opn6_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#B15928",
    width = 0.5
  ) +
  geom_facet(
    panel="rgr",
    data=rgr_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#C2A5CF",
    width = 0.5
  ) +
  geom_facet(
    panel="rrh",
    data=rrh_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#ABD9E9",
    width = 0.5
  ) +
  geom_facet(
    panel="opn4m1_3",
    data=opn4m1_3_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#35978F",
    width = 0.5
  ) +
  geom_facet(
    panel="opn4m2",
    data=opn4m2_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#238443",
    width = 0.5
  ) +
  geom_facet(
    panel="opn4x",
    data=opn4x_prez_abs,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_comp),
    color="#525252",
    width = 0.5
  ) +
  scale_alpha_manual(values = c(
    Presence=1,
    Absence =0
  )) +
  scale_size(range = c(2, 10))




ggtree(order_tree, size=2,
       branch.length="none") +
  geom_tiplab() +
  xlim(0, 30)



##### Boxplot per subfamily  ---------------------------------

all_opsins_table_count_complete$clade <-
  factor(all_opsins_table_count_complete$clade ,
         levels=c("rh1", "exorh", "rh2", "sws2", "sws1", "lws", "pinopsin", "va", "parapinopsin",
                  "parietopsin", "tmt1", "tmt2", "tmt3","opn3", "opn8a", "opn8b", 
                  "opn8c", "opn7a", "opn7b", "opn5", "opn9", "opn6", "rgr", "rrh", "opn4m1_3", 
                  "opn4m2", "opn4x"))


all_opsins_table_count_complete %>%
  ggplot(., aes(x=clade, y=count, color=clade)) +
  geom_boxplot() +
  scale_color_manual(values = opsins_clades_colors) +
  xlab("Subfamily") +
  ylab("Copy number") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none") +
  ylim(0, 10)



all_opsins_table_count_complete %>%
  ggplot(., aes(x=clade, y=count, color=clade)) +
  geom_violin() +
  scale_color_manual(values = opsins_clades_colors) +
  xlab("Subfamily") +
  ylab("Copy number") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none") +
  ylim(1, 10)



sp_wgd <- species_table %>% dplyr::select(species, WGD)
colnames(sp_wgd) <- c("species", "WGD")
all_ospins_complete_wide_wgd <- left_join(all_ospins_complete_wide, sp_wgd, by="species")

test <- 
  all_ospins_complete_wide_wgd %>%
  dplyr::select(species, WGD, total_visual, total_non_visual)
test_long <- 
  as.data.frame(
    test %>% pivot_longer(!c(species, WGD), names_to = "category", values_to = "count")
  )


test_long %>%
  ggplot(., aes(x=category, y=count)) +
  geom_boxplot() +
  #xlab("Subfamily") +
  ylab("Copy number") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none")


test_long %>%
  filter(category == "total_visual") %>%
  ggplot(., aes(x=category,y=count)) +
  geom_boxplot() +
  ylab("Copy number") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none") 

test_long %>%
  filter(category == "total_visual") %>%
  ggplot(., aes(x=category, y=count)) +
  geom_violin() +
  geom_jitter(shape=16, color="gray", alpha=0.8, position=position_jitter(0.2)) +
  ylab("Copy number") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none")

test_long %>%
  filter(category == "total_visual") %>%
  ggplot(., aes(x=category, y=count)) +
  geom_violin() +
  geom_jitter(shape=16, alpha=0.8, position=position_jitter(0.2), aes(color = WGD)) +
  ylab("Copy number") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none") +
  scale_color_manual(values = c("Yes" = "#fdae61", "No" = "#2c7bb6"))

test_long %>%
  filter(category == "total_non_visual") %>%
  ggplot(., aes(x=category, y=count)) +
  geom_boxplot() +
  ylab("Copy number") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none") 


test_long %>%
  filter(category == "total_non_visual") %>%
  ggplot(., aes(x=category, y=count)) +
  geom_violin() +
  geom_jitter(shape=16, color="gray", alpha=0.8, position=position_jitter(0.2)) +
  ylab("Copy number") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none") 



test_long %>%
  filter(category == "total_non_visual") %>%
  ggplot(., aes(x=category, y=count)) +
  geom_violin() +
  geom_jitter(shape=16, alpha=0.8, position=position_jitter(0.2), aes(color = WGD)) +
  ylab("Copy number") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none") +
  scale_color_manual(values = c("Yes" = "#fdae61", "No" = "#2c7bb6"))



test_long %>%
  filter(category == "total_non_visual") %>%
  filter(WGD == "Yes") %>%
  arrange()


tail(
  all_ospins_complete_wide %>%
  dplyr::select(species, total_visual) %>%
  arrange(total_visual),
  3)
head(
  all_ospins_complete_wide %>%
    dplyr::select(species, total_visual) %>%
    arrange(total_visual),
  3)



tail(
  all_ospins_complete_wide %>%
    dplyr::select(species, total_non_visual) %>%
    arrange(total_non_visual),
  3)
head(
  all_ospins_complete_wide %>%
    dplyr::select(species, total_non_visual) %>%
    arrange(total_non_visual),
  3)

all_ospins_complete_wide %>%
  filter(species == "Polypterus_senegalus")


all_ospins_complete_wide_ecology %>%
  filter(order == "Ophidiiformes")

all_ospins_complete_wide_ecology %>%
  filter(species == "Ilyophis_sp_1_JC-2022") %>%
  pull(order)

##### Plot species tree with the nb of complete opsin --   ---------------------------------

WGD_teleost <- 
  species_table %>%
  filter(species %in% final_dataset) %>%
  filter(infraclass == "Teleostei") %>%
  pull(species)

MRCA_WGD_teleost <- 
  getMRCA(ASTRAL_species_tree_nodelabel, WGD_teleost)



WGD_1_list <- 
  species_table %>%
  filter(species %in% final_dataset) %>%
  filter(family == "Cyprinidae") %>%
  filter(WGD == "Yes") %>%
  filter(! species %in% c("Barbus_barbus", "Oxygymnocypris_stewartii")) %>%
  pull(species)

MRCA_WGD_1_astral <- 
  getMRCA(ASTRAL_species_tree_nodelabel, WGD_1_list)


WGD_2_list <- 
  species_table %>%
  filter(species %in% final_dataset) %>%
  filter(family == "Salmonidae") %>%
  filter(WGD == "Yes") %>%
  pull(species)

MRCA_WGD_2_astral <- 
  getMRCA(ASTRAL_species_tree_nodelabel, WGD_2_list)



WGD_3_list <- 
  species_table %>%
  filter(species %in% final_dataset) %>%
  filter(family == "Acipenseridae") %>%
  filter(WGD == "Yes") %>%
  pull(species)

WGD_4_list <- 
  species_table %>%
  filter(species %in% final_dataset) %>%
  filter(family == "Polyodontidae") %>%
  filter(WGD == "Yes") %>%
  pull(species)


WGD_5_list <- 
  species_table %>%
  filter(species %in% final_dataset) %>%
  filter(family == "Catostomidae") %>%
  filter(WGD == "Yes") %>%
  pull(species)
MRCA_WGD_5_astral <- 
  getMRCA(ASTRAL_species_tree_nodelabel, WGD_5_list)



species_alone <- c(WGD_4_list, WGD_3_list, "Barbus_barbus", "Oxygymnocypris_stewartii")
misc_table <- as.data.frame(species_alone) %>% mutate(value = 1)
colnames(misc_table) <- c("node","value")
misc_table$node <- as.numeric(misc_table$node)
node_label_corresp <- 
  left_join(ASTRAL_species_tree_nodelabel, misc_table, by = 'node') %>%
  dplyr::select(node, label)
WGD_sp_alone <- 
  node_label_corresp %>%
  filter(label %in% species_alone) %>%
  pull(node)


ASTRAL_species_tree_nodelabel_space <- ASTRAL_species_tree_nodelabel
ASTRAL_species_tree_nodelabel_space$tip.label <- 
  gsub("_", " ", ASTRAL_species_tree_nodelabel_space$tip.label)

species_table_space <- species_table
species_table_space$genome_name <- gsub("_", " ", species_table_space$species)

visual_opsins_table_count_complete_space <- visual_opsins_table_count_complete
visual_opsins_table_count_complete_space$species <- gsub("_", " ", visual_opsins_table_count_complete_space$species)

non_visual_opsins_table_count_complete_space <- non_visual_opsins_table_count_complete
non_visual_opsins_table_count_complete_space$species <- gsub("_", " ", non_visual_opsins_table_count_complete_space$species)



non_visual_opsins_table_count_complete_space$clade <-
  factor(non_visual_opsins_table_count_complete_space$clade ,
         levels=c("exorh", "pinopsin", "va", "parapinopsin",
                  "parietopsin", "tmt1", "tmt2", "tmt3","opn3", "opn8a", "opn8b", 
                  "opn8c", "opn7a", "opn7b", "opn5", "opn9", "opn6", "rgr", "rrh", "opn4m1_3", 
                  "opn4m2", "opn4x"))


visual_opsins_table_count_complete_space$clade <-
  factor(visual_opsins_table_count_complete_space$clade ,
         levels=c("rh1", "rh2", "sws2", "sws1", "lws"))


rownames(species_table_space) <- species_table_space$genome_name
species_table_space <- species_table_space %>% dplyr::select(genome_name, order)


ggtree(ASTRAL_species_tree_nodelabel_space, layout="circular", size=0.2) %<+% species_table_space  +
  geom_tiplab(size=0.3, offset=0.07, fontface=3) +
  aes(color=as.character(order)) +
  geom_point2(aes(subset = node == 582,  x = x - branch.length * 0.5), size = 0.3, colour = 'red', shape=18) +
  geom_point2(aes(subset = node == 997,  x = x - branch.length * 0.5), size = 0.3, colour = 'red', shape=18) +
  geom_point2(aes(subset = node == 604,  x = x - branch.length * 0.5), size = 0.3, colour = 'red', shape=18) +
  geom_point2(aes(subset = node == 32,  x = x - branch.length * 0.5), size = 0.3, colour = 'red', shape=18) +
  geom_point2(aes(subset = node == 33,  x = x - branch.length * 0.5), size = 0.3, colour = 'red', shape=18) +
  geom_point2(aes(subset = node == 532,  x = x - branch.length * 0.5), size = 0.3, colour = 'red', shape=18) +
  geom_point2(aes(subset = node == 533,  x = x - branch.length * 0.5), size = 0.3, colour = 'red', shape=18) +
  geom_point2(aes(subset = node == 539,  x = x - branch.length * 0.5), size = 0.3, colour = 'red', shape=18) +
  #scale_color_manual(values = tribes.colors) +
  #new_scale_colour() + 
  geom_fruit(
    data=visual_opsins_table_count_complete_space,
    geom=geom_bar,
    mapping=aes(y=species, x=count, fill=clade),
    pwidth=0.5,
    stat="identity",
    orientation="y", 
    position="auto",
    offset=0.35,  #distance depuis les tips
    axis.params=list(
      axis="x", 
      text.angle=-45, 
      hjust=0  
    ),
    grid.params=list() 
  ) +
  geom_fruit(
    data=non_visual_opsins_table_count_complete_space,
    geom=geom_bar,
    mapping=aes(y=species, x=count, fill=clade),
    pwidth=0.5,
    stat="identity",
    orientation="y", 
    position="auto",
    offset=0.1,  #distance depuis les tips
    axis.params=list(
      axis="x", 
      text.angle=-45, 
      hjust=0  
    ),
    grid.params=list() 
  ) +
  scale_fill_manual(values = opsins_clades_colors) +
  theme(legend.position = "none") 





#Compute the number of order which form monophyletic clades

species_test <- "Neomul"
misc_table <- as.data.frame(species_test) %>% mutate(value = 1)
colnames(misc_table) <- c("node","value")
misc_table$node <- as.numeric(misc_table$node)
node_label_corresp <- 
  left_join(ASTRAL_species_tree_nodelabel, misc_table, by = 'node') %>%
  dplyr::select(node, label)



order_list <- 
  species_table %>%
  filter(species %in% final_dataset) %>%
  pull(order) %>%
  unique()

Order_mono_df <- as.data.frame(NULL)
for(curr_order in order_list){
  
  curr_sp_list <- 
    species_table %>%
    filter(species %in% final_dataset) %>%
    filter(order == curr_order) %>%
    pull(species)
  
  if(length(curr_sp_list) > 1){
    
    curr_MRCA <- getMRCA(ASTRAL_species_tree_nodelabel, curr_sp_list)
    curr_desc <- Descendants(ASTRAL_species_tree_nodelabel, curr_MRCA)[[1]]
    curr_desc_label <- node_label_corresp %>% filter(node %in% curr_desc) %>% pull(label)
    curr_desc_orders <- species_table %>% filter(species %in% curr_desc_label) %>% pull(order) %>% unique()
      
    if(length(curr_desc_label) > length(curr_sp_list)) {
      curr_state="paraphyletic"
      
    } else {
      
      curr_state="monophyletic"
      
    }
    
  } else {
    
    curr_state="monophyletic"
    
  }
  
  
  curr_df <- 
    as.data.frame(cbind(curr_order, curr_state))
  colnames(curr_df) <- c("order", "state")
  
  Order_mono_df <- rbind(Order_mono_df, curr_df)
  
}


Order_mono_df %>%
  group_by(state) %>%
  summarise(n())


Order_mono_df %>%
  filter(state == "paraphyletic")


##### Plot species tree with the nb of pseudogenes  --   ---------------------------------

visual_opsins_table_count_pseudo <- 
  visual_opsins_table_count %>%
  filter(gene_state == "Pseudogene")

non_visual_opsins_table_count_pseudo <- 
  non_visual_opsins_table_count %>%
  filter(gene_state == "Pseudogene")


visual_opsins_table_count_pseudo_space <- visual_opsins_table_count_pseudo
visual_opsins_table_count_pseudo_space$species <- gsub("_", " ", visual_opsins_table_count_pseudo_space$species)

non_visual_opsins_table_count_pseudo_space <- non_visual_opsins_table_count_pseudo
non_visual_opsins_table_count_pseudo_space$species <- gsub("_", " ", non_visual_opsins_table_count_pseudo_space$species)

ASTRAL_species_tree_nodelabel_space <- ASTRAL_species_tree_nodelabel
ASTRAL_species_tree_nodelabel_space$tip.label <- 
  gsub("_", " ", ASTRAL_species_tree_nodelabel_space$tip.label)

species_table_space <- species_table
species_table_space$genome_name <- gsub("_", " ", species_table_space$species)


non_visual_opsins_table_count_pseudo_space$clade <-
  factor(non_visual_opsins_table_count_pseudo_space$clade ,
         levels=c("exorh", "pinopsin", "va", "parapinopsin",
                  "parietopsin", "tmt1", "tmt2", "tmt3","opn3", "opn8a", "opn8b", 
                  "opn8c", "opn7a", "opn7b", "opn5", "opn9", "opn6", "rgr", "rrh", "opn4m1_3", 
                  "opn4m2", "opn4x"))


visual_opsins_table_count_pseudo_space$clade <-
  factor(visual_opsins_table_count_pseudo_space$clade ,
         levels=c("rh1", "rh2", "sws2", "sws1", "lws"))


rownames(species_table_space) <- species_table_space$species
species_table_space <- species_table_space %>% dplyr::select(species, order)



ggtree(ASTRAL_species_tree_nodelabel_space, layout="circular", size=0.2) %<+% species_table_space  +
  geom_tiplab(size=0.3, offset=0.07, fontface=3) +
  aes(color=order) +
  geom_point2(aes(subset = node == 582,  x = x - branch.length * 0.5), size = 0.3, colour = 'red', shape=18) +
  geom_point2(aes(subset = node == 997,  x = x - branch.length * 0.5), size = 0.3, colour = 'red', shape=18) +
  geom_point2(aes(subset = node == 604,  x = x - branch.length * 0.5), size = 0.3, colour = 'red', shape=18) +
  geom_point2(aes(subset = node == 32,  x = x - branch.length * 0.5), size = 0.3, colour = 'red', shape=18) +
  geom_point2(aes(subset = node == 33,  x = x - branch.length * 0.5), size = 0.3, colour = 'red', shape=18) +
  geom_point2(aes(subset = node == 532,  x = x - branch.length * 0.5), size = 0.3, colour = 'red', shape=18) +
  geom_point2(aes(subset = node == 533,  x = x - branch.length * 0.5), size = 0.3, colour = 'red', shape=18) +
  geom_point2(aes(subset = node == 539,  x = x - branch.length * 0.5), size = 0.3, colour = 'red', shape=18) +
  #scale_color_manual(values = tribes.colors) +
  #new_scale_colour() + 
  geom_fruit(
    data=visual_opsins_table_count_pseudo_space,
    geom=geom_bar,
    mapping=aes(y=species, x=count, fill=clade),
    pwidth=0.5,
    stat="identity",
    orientation="y", 
    position="auto",
    offset=0.35,  #distance depuis les tips
    axis.params=list(
      axis="x", 
      text.angle=-45, 
      hjust=0  
    ),
    grid.params=list() 
  ) +
  geom_fruit(
    data=non_visual_opsins_table_count_pseudo_space,
    geom=geom_bar,
    mapping=aes(y=species, x=count, fill=clade),
    pwidth=0.5,
    stat="identity",
    orientation="y", 
    position="auto",
    offset=0.1,  #distance depuis les tips
    axis.params=list(
      axis="x", 
      text.angle=-45, 
      hjust=0  
    ),
    grid.params=list() 
  ) +
  scale_fill_manual(values = opsins_clades_colors) +
  theme(legend.position = "none") 





##### Prepare pGLS data -------------------


#First combine the different tables
colnames(fishbase_info_filtered) <- 
  c("species", colnames(fishbase_info_filtered)[c(2:length(colnames(fishbase_info_filtered)))])

fishbase_info_filtered$species <- sub(" ", "_", fishbase_info_filtered$species)

all_ospins_complete_wide_ecology <- 
  left_join(all_ospins_complete_wide, fishbase_info_filtered, by="species")

all_ospins_complete_wide_ecology <- 
  left_join(all_ospins_complete_wide_ecology, Diur_Noctu_table, by="species")

all_ospins_complete_wide_ecology <- 
  left_join(all_ospins_complete_wide_ecology, chemoreceptor_table, by="species")

all_ospins_complete_wide_ecology <- 
  left_join(all_ospins_complete_wide_ecology, electric_table, by="species")

all_ospins_complete_wide_ecology <- 
  left_join(all_ospins_complete_wide_ecology, amphibious_table, by="species")

all_ospins_complete_wide_ecology <- 
  left_join(all_ospins_complete_wide_ecology, Depth_dataframe_filter, by="species")


all_ospins_complete_wide_ecology <- 
  left_join(all_ospins_complete_wide_ecology, genome_size_table, by="species")


all_ospins_complete_wide_ecology <- 
  left_join(all_ospins_complete_wide_ecology, OGG_table_df_filtered, by="species")


colnames(species_table) <- c("species", colnames(species_table)[c(2:length(colnames(species_table)))])

all_ospins_complete_wide_ecology <- 
  left_join(all_ospins_complete_wide_ecology, species_table, by="species")



all_ospins_complete_wide_ecology <- 
  all_ospins_complete_wide_ecology %>%
  mutate(Cave_Surface_filtered = case_when(
    Cave_Surface_filter == "Surface" ~ "Surface", 
    Cave_Surface_filter == "Cave" ~ "Cave", 
    is.na(Cave_Surface_filter) ~ "Surface"
  )) %>%
  dplyr::select(-Cave_Surface_filter)


#Remove the only species in the "Brack" category 
all_ospins_complete_wide_ecology[(all_ospins_complete_wide_ecology$species == "Albula_glossodonta"),"Fresh_salt"] <- NA

#Remove the categorie "unknwon" of DemersPelag
all_ospins_complete_wide_ecology[(all_ospins_complete_wide_ecology$species == "Lampris_incognitus"),"DemersPelag"] <- NA
all_ospins_complete_wide_ecology[(all_ospins_complete_wide_ecology$species == "Lampris_megalopsis"),"DemersPelag"] <- NA
all_ospins_complete_wide_ecology[(all_ospins_complete_wide_ecology$species == "Oreonectes_daqikongensis"),"DemersPelag"] <- NA
all_ospins_complete_wide_ecology[(all_ospins_complete_wide_ecology$species == "Sebastes_cheni"),"DemersPelag"] <- NA



all_ospins_complete_wide_ecology <- 
  all_ospins_complete_wide_ecology %>%
  distinct()


#Transform depth to categories

#Set Pseudoliparis sp yap trench as 7,000m

all_ospins_complete_wide_ecology[(all_ospins_complete_wide_ecology$species == "Pseudoliparis_sp_Yap_Trench"),"depth_mean"] <- 7000


all_ospins_complete_wide_ecology <- 
  all_ospins_complete_wide_ecology %>% 
  mutate(depth_category = case_when(
    depth_mean <= 30 ~ "shallow",
    depth_mean > 30 & depth_mean <= 150 ~ "mesophotic",
    depth_mean > 150 & depth_mean <= 300 ~ "rariphotic",
    depth_mean > 300 & depth_mean <= 1000 ~ "mesopelagic",
    depth_mean > 1000  ~ "bathypelagic",
  )) 


all_ospins_complete_wide_ecology <- 
  all_ospins_complete_wide_ecology %>% 
  mutate(depth_category_lily = case_when(
    depth_mean <= 50 ~ "shallow",
    depth_mean > 50 & depth_mean <= 200 ~ "middle",
    depth_mean > 200 ~ "deep"
  )) 


#Prepare Caper data for the two species tree

caper_data_opsins_ecology_ASTRAL <- 
  comparative.data(phy = ASTRAL_species_tree_nodelabel, data = all_ospins_complete_wide_ecology,
                   names.col = species, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)



##### Sequencing techno or date impact -------------------


#Prepare dataframe

all_ospins_complete_wide_ecology_techno <-
  left_join(all_ospins_complete_wide_ecology, sequencing_techno_df, by="assembly_accession")

all_ospins_complete_wide_ecology_techno <-
  left_join(all_ospins_complete_wide_ecology_techno, date_assembly_df, by="assembly_accession")



all_ospins_complete_wide_ecology_techno <- 
  all_ospins_complete_wide_ecology_techno %>%
  dplyr::select(species, assembly_accession, techno, total_visual, total_non_visual, assembly_level, WGD, date)

all_ospins_complete_long_ecology_techno <- 
  all_ospins_complete_wide_ecology_techno %>%
  pivot_longer(!c(species, assembly_accession, techno, assembly_level, WGD, date), names_to = "VO_NVO", values_to = "count")

## Remove contig as it can be either long or short contigs ..


all_ospins_complete_long_ecology_techno <- 
  all_ospins_complete_long_ecology_techno %>%
  mutate(assembly_level_s = if_else(
    assembly_level == "Contig",
    NA,
    assembly_level
  ))


all_ospins_complete_wide_ecology_techno <- 
  all_ospins_complete_wide_ecology_techno %>%
  mutate(assembly_level_s = if_else(
    assembly_level == "Contig",
    NA,
    assembly_level
  ))

#Remove species with additional WGD

all_ospins_complete_wide_ecology_techno <- 
  all_ospins_complete_wide_ecology_techno %>%
  filter(WGD == "No")
all_ospins_complete_wide_ecology_techno <- 
  as.data.frame(all_ospins_complete_wide_ecology_techno)



all_ospins_complete_long_ecology_techno <- 
  all_ospins_complete_long_ecology_techno %>%
  filter(WGD == "No")
all_ospins_complete_long_ecology_techno <- 
  as.data.frame(all_ospins_complete_long_ecology_techno)


all_ospins_complete_wide_ecology_techno <- 
  all_ospins_complete_wide_ecology_techno %>%
  mutate(total_opsins = total_visual + total_non_visual)

#Prepare Caper data

curr_astral_caper <- 
  comparative.data(phy = ASTRAL_species_tree_nodelabel, data = all_ospins_complete_wide_ecology_techno,
                   names.col = species, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)


#Run pGLS between opsin number and assembly level

visual_vs_level <-
    pgls(total_visual ~ assembly_level_s, 
         data = curr_astral_caper, 
         lambda = "ML")
summary(visual_vs_level)

non_visual_vs_level <-
  pgls(total_non_visual ~ assembly_level_s, 
       data = curr_astral_caper, 
       lambda = "ML")
summary(non_visual_vs_level)
  

#Run pGLS between opsin number and sequencing technology

visual_vs_techno <-
  pgls(total_visual ~ techno, 
       data = curr_astral_caper, 
       lambda = "ML")
summary(visual_vs_techno)
non_visual_vs_techno <-
  pgls(total_non_visual ~ techno, 
       data = curr_astral_caper, 
       lambda = "ML")
summary(non_visual_vs_techno)




all_ospins_complete_long_ecology_techno %>%
  filter(! is.na(assembly_level_s)) %>%
  ggplot(., aes(x=VO_NVO, y=count, fill=assembly_level_s)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Chromosome" = "black", "Scaffold" = "gray")) +
  xlab("") +
  ylab("# of genes") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none") 


all_ospins_complete_long_ecology_techno %>%
  filter(! is.na(techno)) %>%
  ggplot(., aes(x=VO_NVO, y=count, fill=techno)) +
  geom_boxplot() +
  scale_fill_manual(values = c("long" = "black", "short" = "gray")) +
  xlab("") +
  ylab("# of genes") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none") 





all_ospins_complete_wide_ecology_techno %>%
  group_by(assembly_level_s) %>%
  summarise(mean_non_visual = mean(total_non_visual),
            mean_visual = mean(total_visual))




#Run pGLS between opsin number and date


visual_vs_date <-
  pgls(total_visual ~ as.factor(date), 
       data = curr_astral_caper, 
       lambda = "ML")
summary(visual_vs_date)
non_visual_vs_date <-
  pgls(total_non_visual ~ as.factor(date), 
       data = curr_astral_caper, 
       lambda = "ML")
summary(non_visual_vs_date)


all_ospins_complete_long_ecology_techno$date <-
  factor(all_ospins_complete_long_ecology_techno$date ,
         levels=c("2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019",
                  "2020", "2021", "2022", "2023"))



all_ospins_complete_long_ecology_techno %>%
  filter(! is.na(date)) %>%
  ggplot(., aes(x=VO_NVO, y=count, fill=as.factor(date))) +
  geom_boxplot() +
  scale_fill_manual(values = c("long" = "black", "short" = "gray")) +
  xlab("") +
  ylab("# of genes") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none") 


##### pGLS between opsins and WGD -------------------


all_ospins_complete_wide_ecology_teleost <- 
  all_ospins_complete_wide_ecology %>%
  filter(infraclass == "Teleostei")

all_ospins_complete_wide_ecology_teleost %>%
  group_by(WGD) %>%
  summarise(n())

curr_astral_caper <- 
  comparative.data(phy = ASTRAL_species_tree_nodelabel, data = all_ospins_complete_wide_ecology_teleost,
                   names.col = species, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)


pgls_visual_complete_WGD <-
  pgls(total_visual ~ WGD, 
       data = curr_astral_caper, 
       lambda = "ML",
       bounds=list(lambda=c(0,1)))
summary(pgls_visual_complete_WGD)

pgls_non_visual_complete_WGD <-
  pgls(total_non_visual ~ WGD, 
       data = curr_astral_caper, 
       lambda = "ML",
       bounds=list(lambda=c(0,1)))
summary(pgls_non_visual_complete_WGD)

test <- all_ospins_complete_wide_ecology_teleost %>%
  dplyr::select(species, WGD, total_visual, total_non_visual)
test_long <- 
  test %>%
  pivot_longer(!c(species, WGD), names_to = "VO_NVO", values_to = "count")




test_long %>%
  ggplot(., aes(x=VO_NVO, y=count, fill=WGD)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") +
  ylab("# of genes") +
  scale_fill_manual(values = c("Yes" = "gray", "No" = "black"))


all_ospins_complete_wide_ecology_teleost %>%
  group_by(WGD) %>%
  summarise(mean_visual = mean(total_visual),
            mean_non_visual = mean(total_non_visual))


## Now with pseudogenes


visual_opsins_table_count_pseudo <- 
  visual_opsins_table_count %>%
  filter(gene_state == "Pseudogene")

non_visual_opsins_table_count_pseudo <- 
  non_visual_opsins_table_count %>%
  filter(gene_state == "Pseudogene")

sp_wgd <- all_ospins_complete_wide_ecology %>% dplyr::select(species, WGD)

visual_pseudo_df <- 
  visual_opsins_table_count_pseudo %>%
  group_by(species) %>%
  summarise(total_visual_opsin = sum(count))
non_visual_pseudo_df <- 
  non_visual_opsins_table_count_pseudo %>%
  group_by(species) %>%
  summarise(total_nonvisual_opsin = sum(count))

all_pseudo_df <- as.data.frame(full_join(visual_pseudo_df, non_visual_pseudo_df, by="species"))

all_ospins_pseudogene_wide_ecology_long <- as.data.frame(left_join(sp_wgd, all_pseudo_df, by="species"))
all_ospins_pseudogene_wide_ecology_long <- mutate_all(all_ospins_pseudogene_wide_ecology_long, ~replace_na(.,0))

curr_astral_caper <- 
  comparative.data(phy = ASTRAL_species_tree_nodelabel, data = all_ospins_pseudogene_wide_ecology_long,
                   names.col = species, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)


all_ospins_pseudogene_wide_ecology_long %>%
  filter(WGD == "Yes")
pgls_visual_pseudo_WGD <-
  pgls(total_visual_opsin ~ WGD, 
       data = curr_astral_caper, 
       lambda = "ML",
       bounds=list(lambda=c(0,1)))
summary(pgls_visual_pseudo_WGD)

pgls_non_visual_pseudo_WGD <-
  pgls(total_nonvisual_opsin ~ WGD, 
       data = curr_astral_caper, 
       lambda = "ML",
       bounds=list(lambda=c(0,1)))
summary(pgls_non_visual_pseudo_WGD)


test_long <- 
  all_ospins_pseudogene_wide_ecology_long %>%
  pivot_longer(!c(species, WGD), names_to = "VO_NVO", values_to = "count")



test_long %>%
  ggplot(., aes(x=VO_NVO, y=count, fill=WGD)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") +
  ylab("# of pseudogenes") +
  scale_fill_manual(values = c("Yes" = "gray", "No" = "black"))


all_ospins_pseudogene_wide_ecology_long %>%
  group_by(WGD) %>%
  summarise(mean_visual = mean(total_visual_opsin),
            mean_non_visual = mean(total_nonvisual_opsin))




##### pGLS between opsins and depth -------------------

responses_variables <- 
  all_ospins_complete_wide %>% 
  dplyr::select(-species) 
responses_variables <- colnames(responses_variables)


Depth_pgls_df <- as.data.frame(NULL)
for(curr_response in responses_variables){
  
  print(curr_response)
  
  curr_data_opsins <- 
    all_ospins_complete_wide_ecology %>%
    dplyr::select(species, curr_response, depth_category) %>%
    filter(! is.na(depth_category))
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
  
  
  curr_astral_caper <- 
    comparative.data(phy = ASTRAL_species_tree_nodelabel, data = curr_data_opsins,
                     names.col = species, vcv = TRUE,
                     na.omit = FALSE, warn.dropped = TRUE)
  
  
  
  if(curr_response %in% c("rh1")){ #reduce the lambda search boundary for failed computations
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
        "depth_category"
      )
    )
  
  colnames(curr_df) <- c("Response", "R2", "pvalue", "lambda","predictor")
  
  Depth_pgls_df <- 
    rbind(Depth_pgls_df, curr_df)
  
  
}






##### pGLS between opsins and depth -- Teleost NO WGD -------------------

responses_variables <- 
  all_ospins_complete_wide %>% 
  dplyr::select(-species) 
responses_variables <- colnames(responses_variables)


Depth_pgls_df <- as.data.frame(NULL)
for(curr_response in responses_variables){
  
  print(curr_response)
  
  curr_data_opsins <- 
    all_ospins_complete_wide_ecology %>%
    filter(infraclass == "Teleostei") %>%
    filter(WGD == "No") %>%
    dplyr::select(species, curr_response, depth_category) %>%
    filter(! is.na(depth_category))
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
  
  
  curr_astral_caper <- 
    comparative.data(phy = ASTRAL_species_tree_nodelabel, data = curr_data_opsins,
                     names.col = species, vcv = TRUE,
                     na.omit = FALSE, warn.dropped = TRUE)
  
  
  
  if(curr_response %in% c("opn4m2", "opn8a")){ #reduce the lambda search boundary for failed computations
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
        "depth_category"
      )
    )
  
  colnames(curr_df) <- c("Response", "R2", "pvalue", "lambda","predictor")
  
  Depth_pgls_df <- 
    rbind(Depth_pgls_df, curr_df)
  
  
}







##### pGLS between opsins and depth -- Teleost Chromosome -------------------

responses_variables <- 
  all_ospins_complete_wide %>% 
  dplyr::select(-species) 
responses_variables <- colnames(responses_variables)



Depth_pgls_df <- as.data.frame(NULL)
for(curr_response in responses_variables){
  
  print(curr_response)
  
  curr_data_opsins <- 
    all_ospins_complete_wide_ecology %>%
    filter(assembly_level == "Chromosome") %>%
    filter(infraclass == "Teleostei") %>%
    filter(WGD == "No") %>%
    dplyr::select(species, curr_response, depth_category) %>%
    filter(! is.na(depth_category))
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
  
  
  curr_astral_caper <- 
    comparative.data(phy = ASTRAL_species_tree_nodelabel, data = curr_data_opsins,
                     names.col = species, vcv = TRUE,
                     na.omit = FALSE, warn.dropped = TRUE)
  
  
  
  if(curr_response %in% c("todetermine")){ #reduce the lambda search boundary for failed computations
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
        "depth_category"
      )
    )
  
  colnames(curr_df) <- c("Response", "R2", "pvalue", "lambda","predictor")
  
  Depth_pgls_df <- 
    rbind(Depth_pgls_df, curr_df)
  
  
}




##### pGLS between opsins and ecology (categorical variables)  -------------------

responses_variables <- 
  all_ospins_complete_wide %>% 
  dplyr::select(-species) 
responses_variables <- colnames(responses_variables)


predictor_variables <- 
  c("Fresh_salt", "DemersPelag", "Food_cat", "Diel_pattern_simplified",
    "Electrogenic_simplified", "EnvTemp", "AnaCat", "Diadromous_cat",
    "Amphibious", "Terrestrial_exploration")

all_ospins_complete_wide_ecology %>% group_by(Fresh_salt) %>% summarise(n())
all_ospins_complete_wide_ecology %>% group_by(Food_cat) %>% summarise(n())
all_ospins_complete_wide_ecology %>% group_by(EnvTemp) %>% summarise(n())
all_ospins_complete_wide_ecology %>% group_by(DemersPelag) %>% summarise(n())
all_ospins_complete_wide_ecology %>% group_by(Diel_pattern_simplified) %>% summarise(n())
all_ospins_complete_wide_ecology %>% group_by(Electrogenic_simplified) %>% summarise(n())
all_ospins_complete_wide_ecology %>% group_by(AnaCat) %>% summarise(n())
all_ospins_complete_wide_ecology %>% group_by(Diadromous_cat) %>% summarise(n())
all_ospins_complete_wide_ecology %>% group_by(Amphibious) %>% summarise(n())
all_ospins_complete_wide_ecology %>% group_by(Terrestrial_exploration) %>% summarise(n())

#list_summaries <- list()
#names_test_v <- c()

ASTRAL_quali_pgls <- as.data.frame(NULL)
for(curr_response in responses_variables){
  for(curr_predictor in predictor_variables){
    
    
    current_test_name <- paste(curr_response, curr_predictor, sep="_")
    print(current_test_name)
    
    curr_data_opsins <- 
      all_ospins_complete_wide_ecology %>%
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
    
    
    curr_astral_caper <- 
      comparative.data(phy = ASTRAL_species_tree_nodelabel, data = curr_data_opsins,
                       names.col = species, vcv = TRUE,
                       na.omit = FALSE, warn.dropped = TRUE)
    
    
    
    if(current_test_name %in% c("parapinopsin_EnvTemp","opn5_DemersPelag","lws_AnaCat","rgr_Terrestrial_exploration","opn8c_Diadromous_cat","total_non_visual_Fresh_salt","total_non_visual_Diel_pattern_simplified","parapinopsin-a_WGD", "va_Food_cat", "total_non_visual_Diel_pattern", "rh2_WGD")){ #reduce the lambda search boundary for failed computations
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
      
      
    ASTRAL_quali_pgls <- 
      rbind(ASTRAL_quali_pgls,
            curr_df)

    
  }
  
}



ASTRAL_quali_pgls <- read.table("ASTRAL_quali_pgls", 
                                header=TRUE, sep="\t")





##### pGLS between opsins and ecology (categorical variables) -- Teleost NO WGD  -------------------

responses_variables <- 
  all_ospins_complete_wide %>% 
  dplyr::select(-species) 
responses_variables <- colnames(responses_variables)



predictor_variables <- 
  c("Fresh_salt", "DemersPelag", "Food_cat", "Diel_pattern_simplified",
    "Electrogenic_simplified", "EnvTemp", "AnaCat", "Diadromous_cat",
    "Amphibious", "Terrestrial_exploration")


all_ospins_complete_wide_ecology$Food_cat %>% unique()
all_ospins_complete_wide_ecology$EnvTemp %>% unique()



#list_summaries <- list()
#names_test_v <- c()


ASTRAL_quali_pgls_teleost <- as.data.frame(NULL)
for(curr_response in responses_variables){
  for(curr_predictor in predictor_variables){
    
    
    current_test_name <- paste(curr_response, curr_predictor, sep="_")
    print(current_test_name)
    
    curr_data_opsins <- 
      all_ospins_complete_wide_ecology %>%
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
    
    curr_astral_caper <- 
      comparative.data(phy = ASTRAL_species_tree_nodelabel, data = curr_data_opsins,
                       names.col = species, vcv = TRUE,
                       na.omit = FALSE, warn.dropped = TRUE)
    
    
    
    if(current_test_name %in% c("lws_DemersPelag","opn9_Terrestrial_exploration","total_opsins_Food_cat","total_visual_DemersPelag","tmt2_DemersPelag","rgr_Electrogenic_simplified","total_non_visual_Diel_pattern_simplified","parapinopsin-a_WGD", "va_Food_cat", "total_non_visual_Diel_pattern", "rh2_WGD", "tmt2_Diel_pattern_simplified")){ #reduce the lambda search boundary for failed computations
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
    
    
    ASTRAL_quali_pgls_teleost <- 
      rbind(ASTRAL_quali_pgls_teleost,
            curr_df)
    
    
  }
  
}


ASTRAL_quali_pgls_teleost <- read.table("ASTRAL_quali_pgls_teleost", 
                                header=TRUE, sep="\t")




##### pGLS between opsins and ecology (categorical variables) -- Teleost Chromosome  -------------------

responses_variables <- 
  all_ospins_complete_wide %>% 
  dplyr::select(-species) 
responses_variables <- colnames(responses_variables)



predictor_variables <- 
  c("Fresh_salt", "DemersPelag", "Food_cat", "Diel_pattern_simplified",
    "Electrogenic_simplified", "EnvTemp", "AnaCat", "Diadromous_cat",
    "Amphibious", "Terrestrial_exploration")


all_ospins_complete_wide_ecology$Food_cat %>% unique()
all_ospins_complete_wide_ecology$EnvTemp %>% unique()

#list_summaries <- list()
#names_test_v <- c()


ASTRAL_quali_pgls_teleost <- as.data.frame(NULL)
for(curr_response in responses_variables){
  for(curr_predictor in predictor_variables){
    
    
    current_test_name <- paste(curr_response, curr_predictor, sep="_")
    print(current_test_name)
    
    curr_data_opsins <- 
      all_ospins_complete_wide_ecology %>%
      filter(assembly_level == "Chromosome") %>%
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
    
    curr_astral_caper <- 
      comparative.data(phy = ASTRAL_species_tree_nodelabel, data = curr_data_opsins,
                       names.col = species, vcv = TRUE,
                       na.omit = FALSE, warn.dropped = TRUE)
    
    
    
    if(current_test_name %in% c("rh2_Diadromous_cat","rrh_Terrestrial_exploration", "opn8a_DemersPelag", "opn4x_Amphibious","lws_DemersPelag","opn9_Terrestrial_exploration","total_opsins_Food_cat","total_visual_DemersPelag","tmt2_DemersPelag","rgr_Electrogenic_simplified","total_non_visual_Diel_pattern_simplified","parapinopsin-a_WGD", "va_Food_cat", "total_non_visual_Diel_pattern", "rh2_WGD", "tmt2_Diel_pattern_simplified")){ #reduce the lambda search boundary for failed computations
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
    
    
    ASTRAL_quali_pgls_teleost <- 
      rbind(ASTRAL_quali_pgls_teleost,
            curr_df)
    
    
  }
  
}



ASTRAL_quali_pgls_teleost <- read.table("ASTRAL_quali_pgls_teleost", 
                                        header=TRUE, sep="\t")





##### pGLS between opsins and ecology (numeric variables)  -------------------

responses_variables <- 
  all_ospins_complete_wide %>% 
  dplyr::select(-species) 
responses_variables <- colnames(responses_variables)



predictor_variables <- 
  c("ED_SL",
    "Absolute_ED", "Absolute_SL")



all_ospins_complete_wide_ecology$nb_taste_receptor <- as.numeric(all_ospins_complete_wide_ecology$nb_taste_receptor)
all_ospins_complete_wide_ecology$nb_olfactory_receptor <- as.numeric(all_ospins_complete_wide_ecology$nb_olfactory_receptor)
all_ospins_complete_wide_ecology$TempPreferred <- as.numeric(all_ospins_complete_wide_ecology$TempPreferred)
all_ospins_complete_wide_ecology$ED_SL <- as.numeric(all_ospins_complete_wide_ecology$ED_SL)
all_ospins_complete_wide_ecology$Absolute_ED <- as.numeric(all_ospins_complete_wide_ecology$Absolute_ED)
all_ospins_complete_wide_ecology$Absolute_SL <- as.numeric(all_ospins_complete_wide_ecology$Absolute_SL)
all_ospins_complete_wide_ecology$Weight <- as.numeric(all_ospins_complete_wide_ecology$Weight)


Numerical_pgls_df <- as.data.frame(NULL)
for(curr_response in responses_variables){
  for(curr_predictor in predictor_variables){
    
    current_test_name <- paste(curr_response, curr_predictor, sep="_")
    print(current_test_name)
    
    curr_data_opsins <- 
      all_ospins_complete_wide_ecology %>%
      dplyr::select(species, curr_response, curr_predictor)
    
    
    colnames(curr_data_opsins) <- c("species", "reponse", "predictor")
    
    curr_astral_caper <- 
      comparative.data(phy = ASTRAL_species_tree_nodelabel, data = curr_data_opsins,
                       names.col = species, vcv = TRUE,
                       na.omit = FALSE, warn.dropped = TRUE)
    
    
    if(current_test_name %in% c("sws2_Weight","sws1_Weight","opn5_TempPreferred","lws_Weight","total_cone_visual_nb_taste_receptor","rh2_nb_olfactory_receptor","opn9_Absolute_SL", "lws_nb_taste_receptor", "lws_Absolute_SL")){ #reduce the lambda search boundary for failed computations
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
          curr_predictor
        )
      )
    
    colnames(curr_df) <- c("Response", "R2", "pvalue", "lambda", "slope","predictor")
    
    Numerical_pgls_df <- 
      rbind(Numerical_pgls_df, curr_df)
  }
  
}



##### pGLS between opsins and ecology (numeric variables) -- Teleost NO WGD  -------------------

responses_variables <- 
  all_ospins_complete_wide %>% 
  dplyr::select(-species) 
responses_variables <- colnames(responses_variables)


predictor_variables <- 
  c("ED_SL",
    "Absolute_ED", "Absolute_SL")



all_ospins_complete_wide_ecology$ED_SL <- as.numeric(all_ospins_complete_wide_ecology$ED_SL)
all_ospins_complete_wide_ecology$Absolute_ED <- as.numeric(all_ospins_complete_wide_ecology$Absolute_ED)
all_ospins_complete_wide_ecology$Absolute_SL <- as.numeric(all_ospins_complete_wide_ecology$Absolute_SL)


Numerical_pgls_df_teleost <- as.data.frame(NULL)
for(curr_response in responses_variables){
  for(curr_predictor in predictor_variables){
    
    current_test_name <- paste(curr_response, curr_predictor, sep="_")
    print(current_test_name)
    
    curr_data_opsins <- 
      all_ospins_complete_wide_ecology %>%
      filter(infraclass == "Teleostei") %>%
      filter(WGD == "No") %>%
      dplyr::select(species, curr_response, curr_predictor)
    
    
    colnames(curr_data_opsins) <- c("species", "reponse", "predictor")
    
    curr_astral_caper <- 
      comparative.data(phy = ASTRAL_species_tree_nodelabel, data = curr_data_opsins,
                       names.col = species, vcv = TRUE,
                       na.omit = FALSE, warn.dropped = TRUE)
    
    
    if(current_test_name %in% c("opn4m2_Absolute_SL","sws2_Weight","sws1_Weight","opn5_TempPreferred","lws_Weight","total_cone_visual_nb_taste_receptor","rh2_nb_olfactory_receptor","opn9_Absolute_SL", "lws_nb_taste_receptor", "lws_Absolute_SL", "opn4m2_ED_SL")){ #reduce the lambda search boundary for failed computations
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
          curr_predictor
        )
      )
    
    colnames(curr_df) <- c("Response", "R2", "pvalue", "lambda", "slope","predictor")
    
    Numerical_pgls_df_teleost <- 
      rbind(Numerical_pgls_df_teleost, curr_df)
  }
  
}



##### pGLS between opsins and ecology (numeric variables) -- Teleost Chromosome -------------------

responses_variables <- 
  all_ospins_complete_wide %>% 
  dplyr::select(-species) 
responses_variables <- colnames(responses_variables)



predictor_variables <- 
  c("ED_SL",
    "Absolute_ED", "Absolute_SL")



all_ospins_complete_wide_ecology$ED_SL <- as.numeric(all_ospins_complete_wide_ecology$ED_SL)
all_ospins_complete_wide_ecology$Absolute_ED <- as.numeric(all_ospins_complete_wide_ecology$Absolute_ED)
all_ospins_complete_wide_ecology$Absolute_SL <- as.numeric(all_ospins_complete_wide_ecology$Absolute_SL)


Numerical_pgls_df_teleost <- as.data.frame(NULL)
for(curr_response in responses_variables){
  for(curr_predictor in predictor_variables){
    
    current_test_name <- paste(curr_response, curr_predictor, sep="_")
    print(current_test_name)
    
    curr_data_opsins <- 
      all_ospins_complete_wide_ecology %>%
      filter(assembly_level == "Chromosome") %>%
      filter(infraclass == "Teleostei") %>%
      filter(WGD == "No") %>%
      dplyr::select(species, curr_response, curr_predictor)
    
    
    colnames(curr_data_opsins) <- c("species", "reponse", "predictor")
    
    curr_astral_caper <- 
      comparative.data(phy = ASTRAL_species_tree_nodelabel, data = curr_data_opsins,
                       names.col = species, vcv = TRUE,
                       na.omit = FALSE, warn.dropped = TRUE)
    
    
    if(current_test_name %in% c("opn4m2_Absolute_SL","sws2_Weight","sws1_Weight","opn5_TempPreferred","lws_Weight","total_cone_visual_nb_taste_receptor","rh2_nb_olfactory_receptor","opn9_Absolute_SL", "lws_nb_taste_receptor", "lws_Absolute_SL", "opn4m2_ED_SL")){ #reduce the lambda search boundary for failed computations
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
          curr_predictor
        )
      )
    
    colnames(curr_df) <- c("Response", "R2", "pvalue", "lambda", "slope","predictor")
    
    Numerical_pgls_df_teleost <- 
      rbind(Numerical_pgls_df_teleost, curr_df)
  }
  
}





##### Ecomorphology - pGLS Summary and analysis  -------------------

#Import tables

astral_pgls_depth <-
  read.table("PGLS_astral_depthCat.tsv",
             sep="\t",
             header=TRUE)

astral_pgls_num <- 
  read.table("PGLS_numerical_Astral.tsv",
             sep="\t",
             header=TRUE)


ASTRAL_quali_pgls <- 
  read.table("ASTRAL_quali_pgls", 
             sep="\t",
             header=TRUE)

ASTRAL_quali_pgls <- 
  ASTRAL_quali_pgls %>% 
  filter(! predictor %in% c("DemersPelag", "Amphibious", "Terrestrial_exploration", "Diadromous_cat"))




ASTRAL_quali_pgls <- 
  ASTRAL_quali_pgls %>%
  mutate(slope = NA)
astral_pgls_depth <- 
  astral_pgls_depth %>%
  mutate(slope = NA)


ASTRAL_quali_pgls <- 
  ASTRAL_quali_pgls %>%
  dplyr::select(Response,R2, pvalue, lambda, slope, predictor)
astral_pgls_depth <- 
  astral_pgls_depth %>%
  dplyr::select(Response,R2, pvalue, lambda, slope, predictor)


colnames(ASTRAL_quali_pgls) <- colnames(astral_pgls_depth)
colnames(astral_pgls_depth) <- colnames(astral_pgls_depth)


#Merge tables

astral_pgls_all_num <- 
  rbind(astral_pgls_depth, astral_pgls_num)

astral_ecology_pgls <- 
  rbind(astral_pgls_all_num, ASTRAL_quali_pgls)


as.data.frame(
  astral_ecology_pgls %>%
    group_by(Response) %>%
    summarise(n()))

as.data.frame(
  astral_ecology_pgls %>%
    group_by(predictor) %>%
    summarise(n()))

astral_ecology_pgls <- 
  astral_ecology_pgls %>%
  filter(Response != "total_cone_visual")


#Add the FDR and Bonferroni p-value 

unique_response <- astral_ecology_pgls$Response %>% unique()
astral_ecology_pgls_corr <- as.data.frame(NULL)
for(curr_response in unique_response){
  
  curr_df <- 
    astral_ecology_pgls %>%
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
  
  astral_ecology_pgls_corr <- 
    rbind(astral_ecology_pgls_corr, curr_df)
  
  
}


## Define significance levels with the difference corrections

astral_ecology_pgls_corr <- 
  astral_ecology_pgls_corr %>%
  mutate(FDR_significant = if_else(
    as.numeric(FDR_pvalue) < 0.05,
    "Significant", 
    "N.S"
  ))

astral_ecology_pgls_corr <- 
  astral_ecology_pgls_corr %>%
  mutate(Bonferroni_significant = if_else(
    as.numeric(Bonferroni_pvalue) < 0.05,
    "Significant", 
    "N.S"
  ))

astral_ecology_pgls_corr <- 
  astral_ecology_pgls_corr %>%
  mutate(significant = if_else(
    as.numeric(pvalue) < 0.05,
    "Significant", 
    "N.S"
  ))


## Define the slope ...
astral_ecology_pgls_corr <- 
  astral_ecology_pgls_corr %>%
  mutate(slope_sign = case_when(
    slope <= 0 ~ "minus",
    slope > 0 ~ "plus",
    is.na(slope) ~ "cat",
  ))


## Define a shape for slopes for future graphics
slope_sign_shapes <- 
  c("minus"="\u25BC",
    "plus"="\u25B2",
    "cat"="\u25FE")

#Lets make a summary graphic 

astral_ecology_pgls_corr$Response <-
  factor(astral_ecology_pgls_corr$Response ,
         levels=c("rh1", "exorh", "rh2", "sws2", "sws1", "lws", "pinopsin", "va", "parapinopsin",
                  "parietopsin", "tmt1", "tmt2", "tmt3","opn3", "opn8b", "opn8a", 
                  "opn8c", "opn7b", "opn7a", "opn5", "opn9", "opn6", "rgr", "rrh", "opn4m1_3", 
                  "opn4m2", "opn4x", "total_non_visual", "total_visual", "total_opsins"))


astral_ecology_pgls_corr$predictor <-
  factor(astral_ecology_pgls_corr$predictor ,
         levels=c("Absolute_SL", "Absolute_ED", 
                  "ED_SL", "Electrogenic_simplified", "Food_cat","AnaCat", 
                  "Diel_pattern_simplified",
                  "EnvTemp", "Fresh_salt",
                  "depth_category"
         ))



Signif_colors <-
  c("N.S" = 0,
    "Significant" = 1)




ggplot(astral_ecology_pgls_corr, 
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

ggplot(astral_ecology_pgls_corr, 
       aes(y=predictor, x=Response, fill= R2, alpha=Bonferroni_significant, color="black")) + 
  geom_tile(color = "black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_minimal() +
  scale_fill_viridis_c(option = "plasma") +
  scale_alpha_manual(values = Signif_colors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("")






##### Ecomorphology - pGLS Summary and analysis -- Teleost NO WGD  -------------------

#Import tables

astral_pgls_depth_teleost <-
  read.table("PGLS_astral_depthCat_teleost.tsv",
             sep="\t",
             header=TRUE)

astral_pgls_num_teleost <- 
  read.table("PGLS_numerical_Astral_teleost.tsv",
             sep="\t",
             header=TRUE)


ASTRAL_quali_pgls_teleost <- 
  read.table("ASTRAL_quali_pgls_teleost", 
             sep="\t",
             header=TRUE)

ASTRAL_quali_pgls_teleost <- 
  ASTRAL_quali_pgls_teleost %>% 
  filter(! predictor %in% c("DemersPelag", "Amphibious", "Terrestrial_exploration", "Diadromous_cat"))




ASTRAL_quali_pgls_teleost <- 
  ASTRAL_quali_pgls_teleost %>%
  mutate(slope = NA)
astral_pgls_depth_teleost <- 
  astral_pgls_depth_teleost %>%
  mutate(slope = NA)


ASTRAL_quali_pgls_teleost <- 
  ASTRAL_quali_pgls_teleost %>%
  dplyr::select(Response,R2, pvalue, lambda, slope, predictor)
astral_pgls_depth_teleost <- 
  astral_pgls_depth_teleost %>%
  dplyr::select(Response,R2, pvalue, lambda, slope, predictor)


colnames(ASTRAL_quali_pgls_teleost) <- colnames(astral_pgls_depth_teleost)
colnames(astral_pgls_depth_teleost) <- colnames(astral_pgls_depth_teleost)


#Merge tables

astral_pgls_all_num_teleost <- 
  rbind(astral_pgls_depth_teleost, astral_pgls_num_teleost)

astral_ecology_pgls_teleost <- 
  rbind(astral_pgls_all_num_teleost, ASTRAL_quali_pgls_teleost)


as.data.frame(
  astral_ecology_pgls_teleost %>%
    group_by(Response) %>%
    summarise(n()))

as.data.frame(
  astral_ecology_pgls_teleost %>%
    group_by(predictor) %>%
    summarise(n()))

astral_ecology_pgls_teleost <- 
  astral_ecology_pgls_teleost %>%
  filter(Response != "total_cone_visual")


#Add the FDR and Bonferroni p-value 

unique_response <- astral_ecology_pgls_teleost$Response %>% unique()
astral_ecology_pgls_corr_teleost <- as.data.frame(NULL)
for(curr_response in unique_response){
  
  curr_df <- 
    astral_ecology_pgls_teleost %>%
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
  
  astral_ecology_pgls_corr_teleost <- 
    rbind(astral_ecology_pgls_corr_teleost, curr_df)
  
  
}


## Define significance levels with the difference corrections

astral_ecology_pgls_corr_teleost <- 
  astral_ecology_pgls_corr_teleost %>%
  mutate(FDR_significant = if_else(
    as.numeric(FDR_pvalue) < 0.05,
    "Significant", 
    "N.S"
  ))

astral_ecology_pgls_corr_teleost <- 
  astral_ecology_pgls_corr_teleost %>%
  mutate(Bonferroni_significant = if_else(
    as.numeric(Bonferroni_pvalue) < 0.05,
    "Significant", 
    "N.S"
  ))

astral_ecology_pgls_corr_teleost <- 
  astral_ecology_pgls_corr_teleost %>%
  mutate(significant = if_else(
    as.numeric(pvalue) < 0.05,
    "Significant", 
    "N.S"
  ))


## Define the slope ...
astral_ecology_pgls_corr_teleost <- 
  astral_ecology_pgls_corr_teleost %>%
  mutate(slope_sign = case_when(
    slope <= 0 ~ "minus",
    slope > 0 ~ "plus",
    is.na(slope) ~ "cat",
  ))


## Define a shape for slopes for future graphics
slope_sign_shapes <- 
  c("minus"="\u25BC",
    "plus"="\u25B2",
    "cat"="\u25FE")

#Lets make a summary graphic 

astral_ecology_pgls_corr_teleost$Response <-
  factor(astral_ecology_pgls_corr_teleost$Response ,
         levels=c("rh1", "exorh", "rh2", "sws2", "sws1", "lws", "pinopsin", "va", "parapinopsin",
                  "parietopsin", "tmt1", "tmt2", "tmt3","opn3", "opn8b", "opn8a", 
                  "opn8c", "opn7b", "opn7a", "opn5", "opn9", "opn6", "rgr", "rrh", "opn4m1_3", 
                  "opn4m2", "opn4x", "total_non_visual", "total_visual", "total_opsins"))


astral_ecology_pgls_corr_teleost$predictor <-
  factor(astral_ecology_pgls_corr_teleost$predictor ,
         levels=c("Absolute_SL", "Absolute_ED", 
                  "ED_SL", "Electrogenic_simplified", "Food_cat","AnaCat", 
                  "Diel_pattern_simplified",
                  "EnvTemp", "Fresh_salt",
                  "depth_category"
                  ))



Signif_colors <-
  c("N.S" = 0,
    "Significant" = 1)




ggplot(astral_ecology_pgls_corr_teleost, 
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

ggplot(astral_ecology_pgls_corr_teleost, 
       aes(y=predictor, x=Response, fill= R2, alpha=Bonferroni_significant, color="black")) + 
  geom_tile(color = "black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_minimal() +
  scale_fill_viridis_c(option = "plasma") +
  scale_alpha_manual(values = Signif_colors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("")





ggplot(astral_ecology_pgls_corr_teleost, 
       aes(y=predictor, x=Response, fill= R2, alpha=significant, color="black")) + 
  geom_tile(color = "black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_minimal() +
  scale_fill_viridis_c(option = "plasma") +
  scale_alpha_manual(values = Signif_colors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("")






##### Ecology graphics  -------------------



all_ospins_complete_wide_ecology %>%
  pull(EnvTemp) %>%
  unique()


significant_correlations <- 
  astral_ecology_pgls_corr %>%
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
  
  curr_all_ospins_complete_wide_ecology <- 
    all_ospins_complete_wide_ecology %>%
    dplyr::select(species, curr_predictor, curr_response)
  colnames(curr_all_ospins_complete_wide_ecology) <- c("species", "predictor", "response")
  
  
  curr_list_name <- paste(curr_response, curr_predictor, sep="_")
  
  if(is.na(slope)){
    
    good_preds <- 
      curr_all_ospins_complete_wide_ecology %>%
      group_by(predictor) %>%
      summarise(count = n()) %>%
      filter(count > 2) %>%
      pull(predictor)
    
    curr_all_ospins_complete_wide_ecology <-
      curr_all_ospins_complete_wide_ecology %>%
      filter(predictor %in% good_preds)
    
    
    p = curr_all_ospins_complete_wide_ecology %>%
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
      lm(data = curr_all_ospins_complete_wide_ecology %>% filter(! is.na(predictor)),
         formula = response ~ as.numeric(predictor))
    curr_function <- GLS_function(curr_lm)
    
    p = curr_all_ospins_complete_wide_ecology %>%
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



##### Ecology graphics -- Teleost no WGD  -------------------



significant_correlations <- 
  astral_ecology_pgls_corr_teleost %>%
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
  
  curr_all_ospins_complete_wide_ecology <- 
    all_ospins_complete_wide_ecology %>%
    filter(infraclass == "Teleostei") %>%
    filter(WGD == "No") %>%
    dplyr::select(species, curr_predictor, curr_response)
  colnames(curr_all_ospins_complete_wide_ecology) <- c("species", "predictor", "response")
  
  
  curr_list_name <- paste(curr_response, curr_predictor, sep="_")
  
  if(is.na(slope)){
    
    
    good_preds <- 
      curr_all_ospins_complete_wide_ecology %>%
      group_by(predictor) %>%
      summarise(count = n()) %>%
      filter(count > 2) %>%
      pull(predictor)
    
    curr_all_ospins_complete_wide_ecology <-
      curr_all_ospins_complete_wide_ecology %>%
      filter(predictor %in% good_preds)
    
    
    p = curr_all_ospins_complete_wide_ecology %>%
      filter(! is.na(predictor)) %>%
      ggplot(., aes(x=predictor, y=response)) +
      geom_violin() + 
      stat_summary(fun=mean, geom="point", shape=18, size=4) +
      theme_classic() +
      xlab(curr_predictor) +
      ylab(curr_response) +
      theme(axis.text=element_text(size=16),
            axis.title=element_text(size=18),
            plot.subtitle=element_text(size=16),
            legend.position="none") 
    
  } else {
    

    
    curr_lm <- 
      lm(data = curr_all_ospins_complete_wide_ecology %>% filter(! is.na(predictor)),
         formula = response ~ as.numeric(predictor))
    curr_function <- GLS_function(curr_lm)
    
    p = curr_all_ospins_complete_wide_ecology %>%
      filter(! is.na(predictor)) %>%
      ggplot(., aes(x=predictor, y=response)) +
      geom_point() +
      theme_classic() +
      stat_function(fun = curr_function, color="black") +
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


all_ospins_complete_wide_ecology %>%
  filter(! is.na(Fresh_salt)) %>%
  dplyr::select(Fresh_salt, sws1) %>%
  group_by(Fresh_salt) %>%
  summarise(mean_sws1 = mean(sws1))

all_ospins_complete_wide_ecology %>%
  filter(! is.na(Fresh_salt)) %>%
  dplyr::select(Fresh_salt, opn9) %>%
  group_by(Fresh_salt) %>%
  summarise(mean_opn9 = mean(opn9))


all_ospins_complete_wide_ecology %>%
  filter(! is.na(Fresh_salt)) %>%
  dplyr::select(Fresh_salt, rh2) %>%
  group_by(Fresh_salt) %>%
  summarise(mean_rh2 = mean(rh2))



##### Ecology graphics -- Teleost no WGD - Better  -------------------


significant_correlations <- 
  astral_ecology_pgls_corr_teleost %>%
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
    all_ospins_complete_wide_ecology %>%
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
    all_ospins_complete_wide_ecology %>%
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


### Salinity 

significant_correlation_salinity <- 
  significant_correlations %>%
  filter(predictor == "Fresh_salt") 


plot_list <- list()
all_plot_names <- c()
for (row in 1:nrow(significant_correlation_salinity)){
  
  curr_lambda <- significant_correlation_salinity[row, "lambda"]
  curr_R2  <- significant_correlation_salinity[row, "R2"]
  curr_response  <- significant_correlation_salinity[row, "Response"]
  curr_predictor  <- significant_correlation_salinity[row, "predictor"]
  
  
  curr_pvalue <- significant_correlation_salinity[row, "FDR_pvalue"]
  slope <- significant_correlation_salinity[row, "slope"]
  
  curr_all_ospins_complete_wide_ecology <- 
    all_ospins_complete_wide_ecology %>%
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
           levels=c("Fresh", "Fresh_Brack", "Fresh_Brack_Salt", "Salt_Brack", "Salt"
           ))
  
  p = curr_all_ospins_complete_wide_ecology %>%
    filter(! is.na(predictor)) %>%
    ggplot(., aes(x=predictor, y=response)) +
    geom_violin() + 
    geom_jitter(shape=16, color="gray", alpha=0.8, position=position_jitter(0.2)) +
    stat_summary(fun=mean, geom="point", shape=18, size=2, color="red") +
    theme_classic() +
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
    all_ospins_complete_wide_ecology %>%
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
    all_ospins_complete_wide_ecology %>%
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
    all_ospins_complete_wide_ecology %>%
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

##### pGLS between non visual and visual opsins  -------------------

pgls_visual_nonvisual <-
  pgls(total_non_visual ~ total_visual, 
       data = caper_data_opsins_ecology_ASTRAL, 
       lambda = "ML",
       bounds=list(lambda=c(0,1)))
summary(pgls_visual_nonvisual)

sum_fit_phy <- summary(pgls_visual_nonvisual)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC(sum_fit_phy$adj.r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
lambda <- as.numeric(sum_fit_phy$param[2])
slope <- formatC(sum_fit_phy$coefficients[2], digits = 3)


lm_nonvisual_visual <- 
  lm(data = all_ospins_complete_wide_ecology,
     formula = as.numeric(total_non_visual) ~ as.numeric(total_visual))
function_nonvisual_visual <- GLS_function(lm_nonvisual_visual)



#Perform the same test with only teleost and w/o WGD please ! 

caper_data_opsins_ecology_ASTRAL_teleost <- 
  comparative.data(phy = ASTRAL_species_tree_nodelabel, data = 
                     all_ospins_complete_wide_ecology %>%
                     filter(infraclass == "Teleostei") %>%
                     filter(WGD == "No"),
                   names.col = species, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)


pgls_visual_nonvisual_teleost <-
  pgls(total_non_visual ~ total_visual, 
       data = caper_data_opsins_ecology_ASTRAL_teleost, 
       lambda = "ML",
       bounds=list(lambda=c(0,1)))
summary(pgls_visual_nonvisual_teleost)


sum_fit_phy <- summary(pgls_visual_nonvisual_teleost)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC(sum_fit_phy$adj.r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
lambda <- as.numeric(sum_fit_phy$param[2])
slope <- formatC(sum_fit_phy$coefficients[2], digits = 3)


lm_nonvisual_visual_teleost <- 
  lm(data = all_ospins_complete_wide_ecology %>%
       filter(infraclass == "Teleostei") %>%
       filter(WGD == "No"),
     formula = as.numeric(total_non_visual) ~ as.numeric(total_visual))
function_nonvisual_visual_teleost <- GLS_function(lm_nonvisual_visual_teleost)




all_ospins_complete_wide_ecology <- 
  all_ospins_complete_wide_ecology %>%
  mutate(teleost_WGD = if_else(
    (infraclass == "Teleostei") & (WGD == "Yes"),
    "Teleost_WGD",
    "Other"
  ))



all_ospins_complete_wide_ecology %>%
  ggplot(., aes(x=total_visual, y=total_non_visual, color=teleost_WGD)) +
  geom_point() +
  theme_classic() +
  stat_function(fun = function_nonvisual_visual, color="gray") +
  stat_function(fun = function_nonvisual_visual_teleost, color="black") +
  #labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ pvalue ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(lambda))) +
  xlab("# of visual opsins") +
  ylab("# of non-visual opsins") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") +
  scale_color_manual(values = c("Teleost_WGD" = "gray", "Other" = "black"))



##### pGLS between opsins subfamilies   -------------------


responses_variables <- 
  all_ospins_complete_wide %>% 
  dplyr::select(-c(species, total_non_visual, total_visual, total_cone_visual, total_opsins))
responses_variables <- colnames(responses_variables)


between_opsins_pgls_df <- as.data.frame(NULL)
between_opsins_pgls_df_genomesize <- as.data.frame(NULL)
between_opsins_pgls_df_genomesize_int <- as.data.frame(NULL)
between_ospins_AIC_df <- as.data.frame(NULL)
best_models_df <- as.data.frame(NULL)

for(curr_response in responses_variables){
  
  predictor_variables <- 
    responses_variables[-which(responses_variables %in% c(curr_response))]
  
  for(curr_predictor in predictor_variables){
    
    current_test_name <- paste(curr_response, curr_predictor, sep="_")
    print(current_test_name)
    
    curr_data_opsins <- 
      all_ospins_complete_wide_ecology %>%
      dplyr::select(species, curr_response, curr_predictor, genome_size_mbp)
    
    
    colnames(curr_data_opsins) <- c("species", "reponse", "predictor", "genome_size")
    
    curr_astral_caper <- 
      comparative.data(phy = ASTRAL_species_tree_nodelabel, data = curr_data_opsins,
                       names.col = species, vcv = TRUE,
                       na.omit = FALSE, warn.dropped = TRUE)
    
    
    #Perform a test without the genome size into account
      #Use trycatch to prevent loop for breaking if pgls return an error
      tryCatch({
        curr_pgls_result <-
          pgls(reponse ~ predictor, 
               data = curr_astral_caper, 
               lambda = "ML",
               bounds=list(lambda=c(0,1)))
      }, error=function(e){})
  
      #Check if there was an error on the pgls computation. If yes, reduce lambda boundaries untill
      #there is no error anymore. 
    
      if(nrow(between_opsins_pgls_df) >= 1){
        sum_fit_phy <- summary(curr_pgls_result)
        current_R2 <- formatC(sum_fit_phy$adj.r.squared, digits = 2)
        current_slope <- formatC(sum_fit_phy$coefficients[2], digits = 3)
        current_lambda <- as.numeric(sum_fit_phy$param[2])
        previous_R2 <- tail(between_opsins_pgls_df, 1) %>% pull(R2)
        previous_slope <- tail(between_opsins_pgls_df, 1) %>% pull(slope)
        previous_lambda <- tail(between_opsins_pgls_df, 1) %>% pull(lambda)
      
        if(current_R2 == previous_R2 & current_slope == previous_slope & current_lambda == previous_lambda){
          tryCatch({
            curr_pgls_result <-
              pgls(reponse ~ predictor, 
                   data = curr_astral_caper, 
                   lambda = "ML",
                   bounds=list(lambda=c(0.7,1)))
          }, error=function(e){})
  
        }
        
        
        sum_fit_phy <- summary(curr_pgls_result)
        current_R2 <- formatC(sum_fit_phy$adj.r.squared, digits = 2)
        current_slope <- formatC(sum_fit_phy$coefficients[2], digits = 3)
        current_lambda <- as.numeric(sum_fit_phy$param[2])
        previous_R2 <- tail(between_opsins_pgls_df, 1) %>% pull(R2)
        previous_slope <- tail(between_opsins_pgls_df, 1) %>% pull(slope)
        previous_lambda <- tail(between_opsins_pgls_df, 1) %>% pull(lambda)
        
        if(current_R2 == previous_R2 & current_slope == previous_slope & current_lambda == previous_lambda){
          tryCatch({
            curr_pgls_result <-
              pgls(reponse ~ predictor, 
                   data = curr_astral_caper, 
                   lambda = "ML",
                   bounds=list(lambda=c(0.9,1)))
          }, error=function(e){})
          
        }
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
          curr_predictor
        )
      )
    
    colnames(curr_df) <- c("Response", "R2", "pvalue", "lambda", "slope","predictor")
    
    between_opsins_pgls_df <- 
      rbind(between_opsins_pgls_df, curr_df)
    
    rm(sum_fit_phy)
    #Perform a test with  genome size into account and no interaction
    
 
    tryCatch({
      curr_pgls_result_gs <-
        pgls(reponse ~ predictor + genome_size, 
             data = curr_astral_caper, 
             lambda = "ML",
             bounds=list(lambda=c(0,1)))
    }, error=function(e){})
    
    
    if(nrow(between_opsins_pgls_df_genomesize) >= 1){
      sum_fit_phy_gs <- summary(curr_pgls_result_gs)
      current_R2 <- formatC(sum_fit_phy_gs$adj.r.squared, digits = 2)
      current_slope <- formatC(sum_fit_phy_gs$coefficients[2], digits = 3)
      current_lambda <- as.numeric(sum_fit_phy_gs$param[2])
      previous_R2 <- tail(between_opsins_pgls_df_genomesize, 1) %>% pull(R2)
      previous_slope <- tail(between_opsins_pgls_df_genomesize, 1) %>% pull(slope_pred1)
      previous_lambda <- tail(between_opsins_pgls_df_genomesize, 1) %>% pull(lambda)
      
      if(current_R2 == previous_R2 & current_slope == previous_slope & current_lambda == previous_lambda){
        tryCatch({
          curr_pgls_result_gs <-
            pgls(reponse ~ predictor + genome_size, 
                 data = curr_astral_caper, 
                 lambda = "ML",
                 bounds=list(lambda=c(0.7,1)))
        }, error=function(e){})
        
      }
      
      
      sum_fit_phy_gs <- summary(curr_pgls_result_gs)
      current_R2 <- formatC(sum_fit_phy_gs$adj.r.squared, digits = 2)
      current_slope <- formatC(sum_fit_phy_gs$coefficients[2], digits = 3)
      current_lambda <- as.numeric(sum_fit_phy_gs$param[2])
      previous_R2 <- tail(between_opsins_pgls_df_genomesize, 1) %>% pull(R2)
      previous_slope <- tail(between_opsins_pgls_df_genomesize, 1) %>% pull(slope_pred1)
      previous_lambda <- tail(between_opsins_pgls_df_genomesize, 1) %>% pull(lambda)
      
      if(current_R2 == previous_R2 & current_slope == previous_slope & current_lambda == previous_lambda){
        tryCatch({
          curr_pgls_result_gs <-
            pgls(reponse ~ predictor + genome_size, 
                 data = curr_astral_caper, 
                 lambda = "ML",
                 bounds=list(lambda=c(0.9,1)))
        }, error=function(e){})
        
      }
    }
    
    
    sum_fit_phy_gs <- summary(curr_pgls_result_gs)
    PGLS_r2_gs =      formatC(sum_fit_phy_gs$adj.r.squared, digits = 2)
    slope_gs <- formatC(sum_fit_phy_gs$coefficients[2], digits = 3)
    lambda_gs <- as.numeric(sum_fit_phy_gs$param[2])
    PGLS_pvalue_gs_pred =  formatC(sum_fit_phy_gs$coefficients[11], digits = 3)
    PGLS_pvalue_gs_gs =  formatC(sum_fit_phy_gs$coefficients[12], digits = 3)
    PGLS_pvalue_gs_model = pf(sum_fit_phy_gs$fstatistic[1], 
                              sum_fit_phy_gs$fstatistic[2], 
                              sum_fit_phy_gs$fstatistic[3], 
                              lower.tail = FALSE)
    
    
    
    if (PGLS_pvalue_gs_pred == "   0"){ PGLS_pvalue_gs_pred = "< 2.2e-16"}
    if (PGLS_pvalue_gs_gs == "   0"){ PGLS_pvalue_gs_gs = "< 2.2e-16"}
    if (PGLS_pvalue_gs_model == "   0"){ PGLS_pvalue_gs_model = "< 2.2e-16"}
    
    print(PGLS_pvalue_gs_model)
    
    curr_df <- 
      as.data.frame(
        cbind(
          curr_response,
          PGLS_r2_gs,
          PGLS_pvalue_gs_pred,
          PGLS_pvalue_gs_gs,
          PGLS_pvalue_gs_model,
          lambda_gs,
          slope_gs,
          curr_predictor
        )
      )
    
    colnames(curr_df) <- c("Response", "R2", "pvalue_pred1", "pvalue_genomesize",
                           "pvalue_model", "lambda", "slope_pred1","predictor")
    
    between_opsins_pgls_df_genomesize <- 
      rbind(between_opsins_pgls_df_genomesize, curr_df)
    
    rm(sum_fit_phy_gs)

    #Perform a test with  genome size into account and interaction ! 
    
 
    tryCatch({
      curr_pgls_result_gs_int <-
        pgls(reponse ~ predictor * genome_size, 
             data = curr_astral_caper, 
             lambda = "ML",
             bounds=list(lambda=c(0,1)))
    }, error=function(e){})
    
    
    if(nrow(between_opsins_pgls_df_genomesize_int) >= 1){
      sum_fit_phy_gs_int <- summary(curr_pgls_result_gs_int)
      current_R2 <- formatC(sum_fit_phy_gs_int$adj.r.squared, digits = 2)
      current_slope <- formatC(sum_fit_phy_gs_int$coefficients[2], digits = 3)
      current_lambda <- as.numeric(sum_fit_phy_gs_int$param[2])
      previous_R2 <- tail(between_opsins_pgls_df_genomesize_int, 1) %>% pull(R2)
      previous_slope <- tail(between_opsins_pgls_df_genomesize_int, 1) %>% pull(slope_pred1)
      previous_lambda <- tail(between_opsins_pgls_df_genomesize_int, 1) %>% pull(lambda)
      
      if(current_R2 == previous_R2 & current_slope == previous_slope & current_lambda == previous_lambda){
        tryCatch({
          curr_pgls_result_gs_int <-
            pgls(reponse ~ predictor * genome_size, 
                 data = curr_astral_caper, 
                 lambda = "ML",
                 bounds=list(lambda=c(0.7,1)))
        }, error=function(e){})
        
      }
      
      
      sum_fit_phy_gs_int <- summary(curr_pgls_result_gs_int)
      current_R2 <- formatC(sum_fit_phy_gs_int$adj.r.squared, digits = 2)
      current_slope <- formatC(sum_fit_phy_gs_int$coefficients[2], digits = 3)
      current_lambda <- as.numeric(sum_fit_phy_gs_int$param[2])
      previous_R2 <- tail(between_opsins_pgls_df_genomesize_int, 1) %>% pull(R2)
      previous_slope <- tail(between_opsins_pgls_df_genomesize_int, 1) %>% pull(slope_pred1)
      previous_lambda <- tail(between_opsins_pgls_df_genomesize_int, 1) %>% pull(lambda)
      
      if(current_R2 == previous_R2 & current_slope == previous_slope & current_lambda == previous_lambda){
        tryCatch({
          curr_pgls_result_gs_int <-
            pgls(reponse ~ predictor * genome_size, 
                 data = curr_astral_caper, 
                 lambda = "ML",
                 bounds=list(lambda=c(0.9,1)))
        }, error=function(e){})
        
      }
    }
    
    
    sum_fit_phy_gs_int <- summary(curr_pgls_result_gs_int)
    PGLS_r2_gs_int =      formatC(sum_fit_phy_gs_int$adj.r.squared, digits = 2)
    slope_gs_int <- formatC(sum_fit_phy_gs_int$coefficients[2], digits = 3)
    lambda_gs_int <- as.numeric(sum_fit_phy_gs_int$param[2])
    PGLS_pvalue_gs_pred_int =  formatC(sum_fit_phy_gs_int$coefficients[14], digits = 3)
    PGLS_pvalue_gs_gs_int =  formatC(sum_fit_phy_gs_int$coefficients[15], digits = 3)
    PGLS_pvalue_interaction =  formatC(sum_fit_phy_gs_int$coefficients[16], digits = 3)
    PGLS_pvalue_gs_model_int = pf(sum_fit_phy_gs_int$fstatistic[1], 
                                  sum_fit_phy_gs_int$fstatistic[2], 
                                  sum_fit_phy_gs_int$fstatistic[3], 
                                  lower.tail = FALSE)
    
    
    
    if (PGLS_pvalue_gs_pred_int == "   0"){ PGLS_pvalue_gs_pred_int = "< 2.2e-16"}
    if (PGLS_pvalue_gs_gs_int == "   0"){ PGLS_pvalue_gs_gs_int = "< 2.2e-16"}
    if (PGLS_pvalue_interaction == "   0"){ PGLS_pvalue_interaction = "< 2.2e-16"}
    if (PGLS_pvalue_gs_model_int == "   0"){ PGLS_pvalue_gs_model_int = "< 2.2e-16"}
    
    print(PGLS_pvalue_gs_model_int)
    
    curr_df <- 
      as.data.frame(
        cbind(
          curr_response,
          PGLS_r2_gs_int,
          PGLS_pvalue_gs_pred_int,
          PGLS_pvalue_gs_gs_int,
          PGLS_pvalue_interaction,
          PGLS_pvalue_gs_model_int,
          lambda_gs_int,
          slope_gs_int,
          curr_predictor
        )
      )
    
    colnames(curr_df) <- c("Response", "R2", "pvalue_pred1", "pvalue_genomesize",
                           "pvalue_interaction",
                           "pvalue_model", "lambda", "slope_pred1","predictor")
    
    between_opsins_pgls_df_genomesize_int <- 
      rbind(between_opsins_pgls_df_genomesize_int, curr_df)
    
    rm(sum_fit_phy_gs_int)

    #Compare the AIC of the three models
    
    curr_AIC_models <- AIC(curr_pgls_result, 
                           curr_pgls_result_gs, 
                           curr_pgls_result_gs_int)
    curr_AIC_models <- tibble::rownames_to_column(curr_AIC_models, "Test")
    
    curr_AIC_models <- 
      curr_AIC_models %>%
      mutate(model = case_when(
        Test == "curr_pgls_result" ~ "Model simple",
        Test == "curr_pgls_result_gs" ~ "Genome size without interaction",
        Test == "curr_pgls_result_gs_int" ~ "Genome size with interaction"
      )) %>%
      dplyr::select(!Test)
    
    
    curr_AIC_models <- 
      curr_AIC_models %>%
      mutate(Response = curr_response,
             predictor = curr_predictor)
    
    curr_AIC_simple <- 
      curr_AIC_models %>%
      filter(model == "Model simple") %>%
      pull(AIC)
    
    
    curr_AIC_models <- 
      curr_AIC_models %>%
      mutate(Delta_AIC = AIC - curr_AIC_simple) #An absolute difference superior to 2 in AIC mean better model
    
    between_ospins_AIC_df <- 
      rbind(between_ospins_AIC_df,
            curr_AIC_models)
    
    
    #Now lets choose the best model. 
    
    temp_best_AIC_model <- 
      head(curr_AIC_models %>%
             arrange(AIC), 1) %>%
      pull(model)
    
    curr_delta <- 
      curr_AIC_models %>%
      filter(model == temp_best_AIC_model) %>%
      pull(Delta_AIC)
    
    
    if(curr_delta < -2){
      
      best_AIC_model <- temp_best_AIC_model
      
    } else {
      
      best_AIC_model <- "Model simple"
    }
    
    curr_best_models_df <- 
      as.data.frame(cbind(curr_response, curr_predictor, best_AIC_model))
    colnames(curr_best_models_df) <- c("Response", "predictor", "best_model")
    
    best_models_df <- 
      rbind(best_models_df,
            curr_best_models_df)
    
    
  }
  
}





between_opsins_pgls_df <- 
  read.table(
    "PGLS_between_opsins_Astral.tsv",
    header=TRUE,
    sep="\t"
  )
between_opsins_pgls_df_genomesize <- 
  read.table(
    "PGLS_between_opsins_genome_size_Astral.tsv",
    header=TRUE,
    sep="\t"
  )
between_opsins_pgls_df_genomesize_int <- 
  read.table(
    "PGLS_between_opsins_genome_size_interaction_Astral.tsv",
    header=TRUE,
    sep="\t"
  )
between_ospins_AIC_df <- 
  read.table(
    "between_ospins_AIC_df.tsv",
    header=TRUE,
    sep="\t"
  )



### Plot a summary of pGLS results between subfamilies

between_opsins_pgls_df$Response <-
  factor(between_opsins_pgls_df$Response ,
         levels=c("rh1", "exorh", "rh2", "sws2", "sws1", "lws", "pinopsin", "va", "parapinopsin",
                  "parietopsin", "tmt1", "tmt2", "tmt3","opn3", "opn8b", "opn8a", 
                  "opn8c", "opn7b", "opn7a", "opn5", "opn9", "opn6", "rgr", "rrh", "opn4m1_3", 
                  "opn4m2", "opn4x"))
between_opsins_pgls_df$predictor <-
  factor(between_opsins_pgls_df$predictor ,
         levels=c("rh1", "exorh", "rh2", "sws2", "sws1", "lws", "pinopsin", "va", "parapinopsin",
                  "parietopsin", "tmt1", "tmt2", "tmt3","opn3", "opn8b", "opn8a", 
                  "opn8c", "opn7b", "opn7a", "opn5", "opn9", "opn6", "rgr", "rrh", "opn4m1_3", 
                  "opn4m2", "opn4x"))

between_opsins_pgls_df <- 
  between_opsins_pgls_df %>%
  mutate(numeric_pvalue = if_else(
    pvalue == "< 2.2e-16",
    "2.2e-16",
    pvalue
  )) %>%
  dplyr::select(-pvalue)

unique_response <- between_opsins_pgls_df$Response %>% unique()
between_opsins_pgls_df_corr <- as.data.frame(NULL)
for(curr_response in unique_response){
  
  curr_df <- 
    between_opsins_pgls_df %>%
    filter(Response == curr_response)
  
  
  initial_pvalues <- 
    as.numeric(
      curr_df %>% pull(numeric_pvalue)
    )
  
  corrected_pvalues_FDR <- p.adjust(initial_pvalues, method = "fdr")
  corrected_pvalues_Bonferoni <- p.adjust(initial_pvalues, method = "bonferroni")
  
  curr_df <- 
    curr_df %>%
    mutate(FDR_pvalue = corrected_pvalues_FDR,
           Bonferroni_pvalue = corrected_pvalues_Bonferoni)
  
  between_opsins_pgls_df_corr <- 
    rbind(between_opsins_pgls_df_corr, curr_df)
  
}

between_opsins_pgls_df <- between_opsins_pgls_df_corr

between_opsins_pgls_df <- 
  between_opsins_pgls_df %>%
  mutate(significant = if_else(
    as.numeric(numeric_pvalue) < 0.05,
    "Significant", 
    "N.S"
  )) %>%
  mutate(significant_FDR = if_else(
    as.numeric(FDR_pvalue) < 0.05,
    "Significant", 
    "N.S"
  )) %>%
  mutate(significant_bonferroni = if_else(
    as.numeric(Bonferroni_pvalue) < 0.05,
    "Significant", 
    "N.S"
  ))

between_opsins_pgls_df <- 
  between_opsins_pgls_df %>%
  mutate(slope_sign = case_when(
    slope <= 0 ~ "minus",
    slope > 0 ~ "plus",
    is.na(slope) ~ "cat",
  ))



slope_sign_shapes <- 
  c("minus"="\u25BC",
    "plus"="\u25B2",
    "cat"="\u25FE")


#transform R2 as well as pvalue to matrixes 


order <- 
  c("rh1", "exorh", "rh2", "sws2", "sws1", "lws", "va", "parapinopsin",
    "parietopsin", "tmt1", "tmt2", "tmt3","opn3", "opn8b", "opn8a", 
    "opn8c", "opn7b", "opn7a", "opn5", "opn9", "opn6", "rgr", "rrh", "opn4m1_3", 
    "opn4m2", "opn4x")

between_opsins_pgls_df <- 
  between_opsins_pgls_df %>%
  mutate(R = if_else(
    slope > 0,
    sqrt(R2),
    -sqrt(R2)
  ))

between_opsins_pgls_df_R2 <- 
  between_opsins_pgls_df %>%
  dplyr::select(Response, predictor, R)
between_opsins_pgls_df_R2 <- 
  as.data.frame(between_opsins_pgls_df_R2 %>%
                  pivot_wider(names_from = predictor, values_from = R, values_fill = 1))
between_opsins_pgls_df_R2 <- mutate_all(between_opsins_pgls_df_R2, ~replace_na(.,0))
rownames(between_opsins_pgls_df_R2) <- between_opsins_pgls_df_R2$Response
between_opsins_pgls_df_R2 <- between_opsins_pgls_df_R2 %>% dplyr::select(-Response)
between_opsins_pgls_df_R2 <- as.matrix(between_opsins_pgls_df_R2)
between_opsins_pgls_df_R2 <- between_opsins_pgls_df_R2[order,]
between_opsins_pgls_df_R2 <- between_opsins_pgls_df_R2[,order]




between_opsins_pgls_df_pvalue <- 
  between_opsins_pgls_df %>%
  dplyr::select(Response, predictor, Bonferroni_pvalue)
between_opsins_pgls_df_pvalue <- 
  as.data.frame(between_opsins_pgls_df_pvalue %>%
                  pivot_wider(names_from = predictor, values_from = Bonferroni_pvalue, values_fill = 1))
rownames(between_opsins_pgls_df_pvalue) <- between_opsins_pgls_df_pvalue$Response
between_opsins_pgls_df_pvalue <- between_opsins_pgls_df_pvalue %>% dplyr::select(-Response)
between_opsins_pgls_df_pvalue <- as.matrix(between_opsins_pgls_df_pvalue)
between_opsins_pgls_df_pvalue <- between_opsins_pgls_df_pvalue[order,]
between_opsins_pgls_df_pvalue <- between_opsins_pgls_df_pvalue[,order]



ggcorrplot(
  corr = between_opsins_pgls_df_R2,
  p.mat = between_opsins_pgls_df_pvalue,
  sig.level = 0.05,
  colors = c("#6D9EC1", "white", "darkred"),
  type = "lower", insig = "blank")






##### pGLS between opsins subfamilies -- only teleost   -------------------


responses_variables <- 
  all_ospins_complete_wide %>% 
  dplyr::select(-c(species, total_non_visual, total_visual, total_cone_visual))
responses_variables <- colnames(responses_variables)
responses_variables <- responses_variables[! responses_variables %in% c("pinopsin", "total_opsins")]

between_opsins_pgls_df_teleost <- as.data.frame(NULL)

for(curr_response in responses_variables){
  
  predictor_variables <- 
    responses_variables[-which(responses_variables %in% c(curr_response))]
  
  for(curr_predictor in predictor_variables){
    
    current_test_name <- paste(curr_response, curr_predictor, sep="_")
    print(current_test_name)
    
    curr_data_opsins <- 
      all_ospins_complete_wide_ecology %>%
      filter(infraclass == "Teleostei") %>%
      filter(WGD == "No") %>%
      dplyr::select(species, curr_response, curr_predictor, genome_size_mbp)
    
    
    colnames(curr_data_opsins) <- c("species", "reponse", "predictor", "genome_size")
    
    curr_astral_caper <- 
      comparative.data(phy = ASTRAL_species_tree_nodelabel, data = curr_data_opsins,
                       names.col = species, vcv = TRUE,
                       na.omit = FALSE, warn.dropped = TRUE)
    
    
    #Perform a test without the genome size into account
    #Use trycatch to prevent loop for breaking if pgls return an error
    tryCatch({
      curr_pgls_result <-
        pgls(reponse ~ predictor, 
             data = curr_astral_caper, 
             lambda = "ML",
             bounds=list(lambda=c(0,1)))
    }, error=function(e){})
    
    #Check if there was an error on the pgls computation. If yes, reduce lambda boundaries untill
    #there is no error anymore. 
    
    if(nrow(between_opsins_pgls_df_teleost) >= 1){
      sum_fit_phy <- summary(curr_pgls_result)
      current_R2 <- formatC(sum_fit_phy$adj.r.squared, digits = 2)
      current_slope <- formatC(sum_fit_phy$coefficients[2], digits = 3)
      current_lambda <- as.numeric(sum_fit_phy$param[2])
      previous_R2 <- tail(between_opsins_pgls_df_teleost, 1) %>% pull(R2)
      previous_slope <- tail(between_opsins_pgls_df_teleost, 1) %>% pull(slope)
      previous_lambda <- tail(between_opsins_pgls_df_teleost, 1) %>% pull(lambda)
      
      if(current_R2 == previous_R2 & current_slope == previous_slope & current_lambda == previous_lambda){
        tryCatch({
          curr_pgls_result <-
            pgls(reponse ~ predictor, 
                 data = curr_astral_caper, 
                 lambda = "ML",
                 bounds=list(lambda=c(0.7,1)))
        }, error=function(e){})
        
      }
      
      
      sum_fit_phy <- summary(curr_pgls_result)
      current_R2 <- formatC(sum_fit_phy$adj.r.squared, digits = 2)
      current_slope <- formatC(sum_fit_phy$coefficients[2], digits = 3)
      current_lambda <- as.numeric(sum_fit_phy$param[2])
      previous_R2 <- tail(between_opsins_pgls_df_teleost, 1) %>% pull(R2)
      previous_slope <- tail(between_opsins_pgls_df_teleost, 1) %>% pull(slope)
      previous_lambda <- tail(between_opsins_pgls_df_teleost, 1) %>% pull(lambda)
      
      if(current_R2 == previous_R2 & current_slope == previous_slope & current_lambda == previous_lambda){
        tryCatch({
          curr_pgls_result <-
            pgls(reponse ~ predictor, 
                 data = curr_astral_caper, 
                 lambda = "ML",
                 bounds=list(lambda=c(0.9,1)))
        }, error=function(e){})
        
      }
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
          curr_predictor
        )
      )
    
    colnames(curr_df) <- c("Response", "R2", "pvalue", "lambda", "slope","predictor")
    
    between_opsins_pgls_df_teleost <- 
      rbind(between_opsins_pgls_df_teleost, curr_df)
    
  }
  
}





between_opsins_pgls_df_teleost <- 
  read.table(
    "PGLS_between_opsins_Astral_onlyteleost.tsv",
    header=TRUE,
    sep="\t"
  )


between_opsins_pgls_df_teleost %>%
  group_by(Response) %>%
  summarise(n())
between_opsins_pgls_df_teleost %>%
  group_by(predictor) %>%
  summarise(n())


### Correct-pvalue

between_opsins_pgls_df_teleost$Response <-
  factor(between_opsins_pgls_df_teleost$Response ,
         levels=c("rh1", "exorh", "rh2", "sws2", "sws1", "lws", "va", "parapinopsin",
                  "parietopsin", "tmt1", "tmt2", "tmt3","opn3", "opn8b", "opn8a", 
                  "opn8c", "opn7b", "opn7a", "opn5", "opn9", "opn6", "rgr", "rrh", "opn4m1_3", 
                  "opn4m2", "opn4x"))
between_opsins_pgls_df_teleost$predictor <-
  factor(between_opsins_pgls_df_teleost$predictor ,
         levels=c("rh1", "exorh", "rh2", "sws2", "sws1", "lws", "va", "parapinopsin",
                  "parietopsin", "tmt1", "tmt2", "tmt3","opn3", "opn8b", "opn8a", 
                  "opn8c", "opn7b", "opn7a", "opn5", "opn9", "opn6", "rgr", "rrh", "opn4m1_3", 
                  "opn4m2", "opn4x"))

between_opsins_pgls_df_teleost <- 
  between_opsins_pgls_df_teleost %>%
  mutate(numeric_pvalue = if_else(
    pvalue == "< 2.2e-16",
    "2.2e-16",
    pvalue
  )) %>%
  dplyr::select(-pvalue)

unique_response <- between_opsins_pgls_df_teleost$Response %>% unique()
between_opsins_pgls_df_corr <- as.data.frame(NULL)
for(curr_response in unique_response){
  
  curr_df <- 
    between_opsins_pgls_df_teleost %>%
    filter(Response == curr_response)
  
  
  initial_pvalues <- 
    as.numeric(
      curr_df %>% pull(numeric_pvalue)
    )
  
  corrected_pvalues_FDR <- p.adjust(initial_pvalues, method = "fdr")
  corrected_pvalues_Bonferoni <- p.adjust(initial_pvalues, method = "bonferroni")
  
  curr_df <- 
    curr_df %>%
    mutate(FDR_pvalue = corrected_pvalues_FDR,
           Bonferroni_pvalue = corrected_pvalues_Bonferoni)
  
  between_opsins_pgls_df_corr <- 
    rbind(between_opsins_pgls_df_corr, curr_df)
  
}

between_opsins_pgls_df_teleost <- between_opsins_pgls_df_corr

between_opsins_pgls_df_teleost <- 
  between_opsins_pgls_df_teleost %>%
  mutate(significant = if_else(
    as.numeric(numeric_pvalue) < 0.05,
    "Significant", 
    "N.S"
  )) %>%
  mutate(significant_FDR = if_else(
    as.numeric(FDR_pvalue) < 0.05,
    "Significant", 
    "N.S"
  )) %>%
  mutate(significant_bonferroni = if_else(
    as.numeric(Bonferroni_pvalue) < 0.05,
    "Significant", 
    "N.S"
  ))

between_opsins_pgls_df_teleost <- 
  between_opsins_pgls_df_teleost %>%
  mutate(slope_sign = case_when(
    slope <= 0 ~ "minus",
    slope > 0 ~ "plus",
    is.na(slope) ~ "cat",
  ))



slope_sign_shapes <- 
  c("minus"="\u25BC",
    "plus"="\u25B2",
    "cat"="\u25FE")



#transform R2 as well as pvalue to matrixes 


order <- 
  c("rh1", "exorh", "rh2", "sws2", "sws1", "lws", "va", "parapinopsin",
    "parietopsin", "tmt1", "tmt2", "tmt3","opn3", "opn8b", "opn8a", 
    "opn8c", "opn7b", "opn7a", "opn5", "opn9", "opn6", "rgr", "rrh", "opn4m1_3", 
    "opn4m2", "opn4x")

between_opsins_pgls_df_teleost <- 
  between_opsins_pgls_df_teleost %>%
  mutate(R = if_else(
    slope > 0,
    sqrt(R2),
    -sqrt(R2)
  ))

between_opsins_pgls_df_teleost_R2 <- 
  between_opsins_pgls_df_teleost %>%
  dplyr::select(Response, predictor, R)
between_opsins_pgls_df_teleost_R2 <- 
  as.data.frame(between_opsins_pgls_df_teleost_R2 %>%
  pivot_wider(names_from = predictor, values_from = R, values_fill = 1))
between_opsins_pgls_df_teleost_R2 <- mutate_all(between_opsins_pgls_df_teleost_R2, ~replace_na(.,0))
rownames(between_opsins_pgls_df_teleost_R2) <- between_opsins_pgls_df_teleost_R2$Response
between_opsins_pgls_df_teleost_R2 <- between_opsins_pgls_df_teleost_R2 %>% dplyr::select(-Response)
between_opsins_pgls_df_teleost_R2 <- as.matrix(between_opsins_pgls_df_teleost_R2)
between_opsins_pgls_df_teleost_R2 <- between_opsins_pgls_df_teleost_R2[order,]
between_opsins_pgls_df_teleost_R2 <- between_opsins_pgls_df_teleost_R2[,order]




between_opsins_pgls_df_teleost_pvalue <- 
  between_opsins_pgls_df_teleost %>%
  dplyr::select(Response, predictor, Bonferroni_pvalue)
between_opsins_pgls_df_teleost_pvalue <- 
  as.data.frame(between_opsins_pgls_df_teleost_pvalue %>%
                  pivot_wider(names_from = predictor, values_from = Bonferroni_pvalue, values_fill = 1))
rownames(between_opsins_pgls_df_teleost_pvalue) <- between_opsins_pgls_df_teleost_pvalue$Response
between_opsins_pgls_df_teleost_pvalue <- between_opsins_pgls_df_teleost_pvalue %>% dplyr::select(-Response)
between_opsins_pgls_df_teleost_pvalue <- as.matrix(between_opsins_pgls_df_teleost_pvalue)
between_opsins_pgls_df_teleost_pvalue <- between_opsins_pgls_df_teleost_pvalue[order,]
between_opsins_pgls_df_teleost_pvalue <- between_opsins_pgls_df_teleost_pvalue[,order]




ggcorrplot(
  corr = between_opsins_pgls_df_teleost_R2,
  p.mat = between_opsins_pgls_df_teleost_pvalue,
  sig.level = 0.05,
  colors = c("#6D9EC1", "white", "darkred"),
  type = "lower", insig = "blank")





##### Complete / Incomplete / Pseudo corr ? -------------------


#Prepare data and caper object
all_opsins_table_count_summary <- 
  all_opsins_table_count %>%
  group_by(species, gene_state) %>%
  summarise(count_total = sum(count))
all_opsins_table_count_summary_wide <- 
  all_opsins_table_count_summary %>%
  pivot_wider(names_from = gene_state, values_from = count_total)
all_opsins_table_count_summary_wide[is.na(all_opsins_table_count_summary_wide)] <- 0
all_opsins_table_count_summary_wide <- 
  left_join(all_opsins_table_count_summary_wide, species_table, by="species")
all_opsins_table_count_summary_wide <- as.data.frame(all_opsins_table_count_summary_wide)

all_opsins_table_count_summary_wide_teleost <- 
  all_opsins_table_count_summary_wide %>%
  filter(infraclass == "Teleostei") %>%
  filter(WGD == "No")

caper_data_opsins_ecology_ASTRAL_pseudo <- 
  comparative.data(phy = ASTRAL_species_tree_nodelabel, data = all_opsins_table_count_summary_wide,
                   names.col = species, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)


caper_data_opsins_ecology_ASTRAL_pseudo_onlyteleost <- 
  comparative.data(phy = ASTRAL_species_tree_nodelabel, data = all_opsins_table_count_summary_wide_teleost,
                   names.col = species, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)

#Run the pGLS and graph 

pgls_pseudo_complete <-
  pgls(Pseudogene ~ Complete, 
       data = caper_data_opsins_ecology_ASTRAL_pseudo, 
       lambda = "ML",
       bounds=list(lambda=c(0,1)))
pgls_pseudo_complete_teleost <-
  pgls(Pseudogene ~ Complete, 
       data = caper_data_opsins_ecology_ASTRAL_pseudo_onlyteleost, 
       lambda = "ML",
       bounds=list(lambda=c(0,1)))

summary(pgls_pseudo_complete)
summary(pgls_pseudo_complete_teleost)


sum_fit_phy <- summary(pgls_pseudo_complete)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC(sum_fit_phy$adj.r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
lambda <- as.numeric(sum_fit_phy$param[2])
slope <- formatC(sum_fit_phy$coefficients[2], digits = 3)


curr_lm <- 
  lm(data = all_opsins_table_count_summary_wide,
     formula = Pseudogene ~ Complete)
#curr_function <- GLS_function(curr_lm)
GLS_cc <- coef(summary(pgls_pseudo_complete))
curr_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

curr_lm_teleost <- 
  lm(data = all_opsins_table_count_summary_wide %>% filter(infraclass == "Teleostei") %>% filter(WGD == "No"),
     formula = Pseudogene ~ Complete)
curr_function_teleost <- GLS_function(curr_lm_teleost)
#GLS_cc <- coef(summary(pgls_pseudo_complete_teleost))
#curr_function_teleost <- function(x) GLS_cc[1] + GLS_cc[2]*x


all_opsins_table_count_summary_wide <- 
  all_opsins_table_count_summary_wide %>%
  mutate(teleost_WGD = if_else(
    (infraclass == "Teleostei") & (WGD == "Yes"),
    "Teleost_WGD",
    "Other"
  ))



all_opsins_table_count_summary_wide %>%
  ggplot(., aes(x=Complete, y=Pseudogene, color=teleost_WGD)) +
  geom_point() +
  theme_classic() +
  stat_function(fun = curr_function, color="gray") +
  stat_function(fun = curr_function_teleost, color="black") +
  xlab("# of genes") +
  ylab("# of pseudogenes") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none")  +
  scale_color_manual(values = c("Teleost_WGD" = "gray", "Other" = "black")) +
  ylim(0, 30)











all_opsins_table_count_summary_wide %>%
  ggplot(., aes(x=Pseudogene, y=Incomplete)) + 
  geom_point()




tail(all_opsins_table_count_summary_wide %>%
  arrange(Pseudogene), 50) %>%
  dplyr::select(species, Complete, Pseudogene, Incomplete, WGD)


head(all_opsins_table_count_summary_wide %>%
       arrange(Pseudogene), 50) %>%
  dplyr::select(species, Complete, Pseudogene, Incomplete, WGD)

all_ospins_complete_wide_ecology %>%
  filter(species == "Dissostichus_mawsoni") %>%
  dplyr::select(species, depth_mean)



all_ospins_complete_wide_ecology %>%
  filter(! is.na(depth_mean)) %>%
  arrange(depth_mean) %>%
  dplyr::select(species, depth_mean)



all_ospins_complete_wide_ecology %>%
  filter(DemersPelag == "bathydemersal") %>%
  pull(species)


all_ospins_complete_wide_ecology %>%
  filter(species == "Dissostichus_mawsoni") %>%
  pull(DemersPelag)



all_opsins_table_count_summary_wide %>%
  filter(infraclass == "Teleostei") %>%
  ggplot(., aes(x=WGD, y=Pseudogene)) +
  geom_violin() +
  geom_jitter(shape=16, color="gray", alpha=0.8, position=position_jitter(0.2)) +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") +
    ylab("# of pseudogenes")
  





all_opsins_table_count_summary_wide_teleost_all <- 
  all_opsins_table_count_summary_wide %>%
  filter(infraclass == "Teleostei") 

caper_data_opsins_ecology_ASTRAL_pseudo_onlyteleost_all <- 
  comparative.data(phy = ASTRAL_species_tree_nodelabel, data = all_opsins_table_count_summary_wide_teleost_all,
                   names.col = species, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)

pgls_pseudo_WGD <-
  pgls(Pseudogene ~ WGD, 
       data = caper_data_opsins_ecology_ASTRAL_pseudo_onlyteleost_all, 
       lambda = "ML",
       bounds=list(lambda=c(0,1)))



pgls_pseudo_WGD <-
  pgls(Complete ~ WGD, 
       data = caper_data_opsins_ecology_ASTRAL_pseudo_onlyteleost_all, 
       lambda = "ML",
       bounds=list(lambda=c(0,1)))




