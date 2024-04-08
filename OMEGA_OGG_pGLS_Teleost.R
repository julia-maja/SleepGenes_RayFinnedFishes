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
library(RColorBrewer)

split_tibble <- function(tibble, column = 'col') {
  tibble %>% split(., .[,column]) %>% lapply(., function(x) x[,setdiff(names(x),column)])
}



#### My functions + Colors palettes  ---------------------------------

args = commandArgs(trailingOnly=TRUE)


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


#Barbus_barbus => recent specific whole genome duplication 
#ref1 = https://serval.unil.ch/resource/serval:BIB_8E5536BAC604.P001/REF.pdf
#ref2 = https://academic.oup.com/gbe/article/13/7/evab131/6320064


#Oxygymnocypris_stewartii => recent specific whole genome duplicatiion
#ref1 = https://www.nature.com/articles/sdata20199/tables/9



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

##### Data load - Species trees   ---------------------------------

ASTRAL_species_tree_nodelabel <- 
  read.tree("AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel")
ASTRAL_species_tree_nodelabel$edge.length <- ASTRAL_species_tree_nodelabel$edge.length * 1000


ASTRAL_species_tree_nodelabel_teleost <- 
  
  
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


all_ospins_complete_wide_ecology %>% 
  filter(! is.na(depth_mean)) %>%
  group_by(depth_mean_source) %>% summarise(n())


#Prepare Caper data for the two species tree

caper_data_opsins_ecology_ASTRAL <- 
  comparative.data(phy = ASTRAL_species_tree_nodelabel, data = all_ospins_complete_wide_ecology,
                   names.col = species, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)



##### pGLS between opsins and OGG omega  -------------------


mean_omega_per_sp_per_ogg_df_wide <- 
  read.table("mean_omega_per_sp_per_ogg_df_wide.tsv",
             header=TRUE,
             sep=",")

all_ospins_complete_wide_ecology <- 
  left_join(mean_omega_per_sp_per_ogg_df_wide, all_ospins_complete_wide_ecology, by="species")

mean_omega_per_sp_per_ogg_df




list_OGG_sup20 <- 
  scan("list_OGG_sup20.txt",
       what="character")


#Lets launch pGLS



list_OGG_sup20 <- c(list_OGG_sup20[-1], list_OGG_sup20[1])

responses_variables <- c("total_non_visual", "total_visual", "total_opsins")


opsins_OGG_pgls_df <- as.data.frame(NULL)
for(curr_response in responses_variables){
  
  for(curr_predictor in list_OGG_sup20){
    
    print(curr_predictor)
    
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
    
    if(nrow(opsins_OGG_pgls_df) >= 1){
      sum_fit_phy <- summary(curr_pgls_result)
      current_R2 <- formatC(sum_fit_phy$adj.r.squared, digits = 2)
      current_slope <- formatC(sum_fit_phy$coefficients[2], digits = 3)
      current_lambda <- as.numeric(sum_fit_phy$param[2])
      previous_R2 <- tail(opsins_OGG_pgls_df, 1) %>% pull(R2)
      previous_slope <- tail(opsins_OGG_pgls_df, 1) %>% pull(slope)
      previous_lambda <- tail(opsins_OGG_pgls_df, 1) %>% pull(lambda)
      
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
      previous_R2 <- tail(opsins_OGG_pgls_df, 1) %>% pull(R2)
      previous_slope <- tail(opsins_OGG_pgls_df, 1) %>% pull(slope)
      previous_lambda <- tail(opsins_OGG_pgls_df, 1) %>% pull(lambda)
      
      if(current_R2 == previous_R2 & current_slope == previous_slope & current_lambda == previous_lambda){
        tryCatch({
          curr_pgls_result <-
            pgls(reponse ~ predictor, 
                 data = curr_astral_caper, 
                 lambda = "ML",
                 bounds=list(lambda=c(0.9,1)))
        }, error=function(e){})
        
      }
      
      sum_fit_phy <- summary(curr_pgls_result)
      current_R2 <- formatC(sum_fit_phy$adj.r.squared, digits = 2)
      current_slope <- formatC(sum_fit_phy$coefficients[2], digits = 3)
      current_lambda <- as.numeric(sum_fit_phy$param[2])
      previous_R2 <- tail(opsins_OGG_pgls_df, 1) %>% pull(R2)
      previous_slope <- tail(opsins_OGG_pgls_df, 1) %>% pull(slope)
      previous_lambda <- tail(opsins_OGG_pgls_df, 1) %>% pull(lambda)
      
      
      if(current_R2 == previous_R2 & current_slope == previous_slope & current_lambda == previous_lambda){
        tryCatch({
          curr_pgls_result <-
            pgls(reponse ~ predictor, 
                 data = curr_astral_caper, 
                 lambda = 1)
        }, error=function(e){})
        
      }
      
    }
    
    if(exists("curr_pgls_result")){
      
      sum_fit_phy <- summary(curr_pgls_result)
      PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
      PGLS_r2 =      formatC(sum_fit_phy$adj.r.squared, digits = 2)
      PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
      if (PGLS_pvalue == "   0"){ PGLS_pvalue = 2.2e-16}
      lambda <- as.numeric(sum_fit_phy$param[2])
      slope <- formatC(sum_fit_phy$coefficients[2], digits = 3)
      
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
      
      opsins_OGG_pgls_df <- 
        rbind(opsins_OGG_pgls_df, curr_df)
      
      rm(sum_fit_phy)
      
      
    }
  }
  
}







