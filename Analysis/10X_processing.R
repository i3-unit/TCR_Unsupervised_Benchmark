# Importation and processing of 10X Genomics data from Donor1

source("../Functions/general_functions.R")
source("../Functions/subset_functions.R")

#library
library(readr)
library(Seurat)
library(dplyr)
library(reshape2)
library(tidyr)
library(data.table)

#Import Donor 1 data
donor1 <- read_delim(file="../Data/10X/Donor1/Contig_vdj_TRB_filtered.csv", delim = ";")
matrix <- Read10X_h5("../Data/10X/Donor1/vdj_v1_hs_aggregated_donor1_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)

#Processing
donor1_filtered <- filter(donor1, is_cell==TRUE & high_confidence ==TRUE & productive =="True" & full_length ==TRUE) 
donor1_filtered <- donor1_filtered %>% select(-c("is_cell", "high_confidence", "productive", "full_length"))
donor1_filtered <- filter(donor1_filtered, chain != "Multi")
donor1_filtered <- donor1_filtered[unsplit(table(donor1_filtered$barcode), donor1_filtered$barcode) < 4, ]#take cells with <4 CDR3
donor1_filtered <- group_by(donor1_filtered, barcode, chain) %>% top_n(1, reads)#if 2 TRA with 1 TRB in one cell, take TRA with highest nb of reads
donor1_filtered <- subset(donor1_filtered, nchar(as.character(donor1_filtered$cdr3))<=22 & nchar(as.character(donor1_filtered$cdr3))>=6) #take CDR3 with size between 6 & 22 AA

tmp <- t(data.frame((matrix$`Antibody Capture`))) %>% data.frame()#transpose of Antibody capture matrix
tmp_2 <- rownames_as_first_column(tmp, "cells")
tmp_2$cells <- gsub(pattern = "\\.",replacement = "-", x = tmp_2$cells)#replace "." by "-" in cell names

combined_df <- merge(donor1_filtered, tmp_2, by.x="barcode", by.y="cells") #combine dfs
#combined_df_TRB <- filter(combined_df, chain=="TRB") 

##bind epitopes
epitope <- relocate(combined_df, cdr3, .before = VTEHDTLLY_IE.1_CMV)
epitope <- epitope[,28:78] #df with cdr3 and score binding

#mean on binding scores of a same CDR3
plyr_base_fn <- function(){
  group_by(epitope, cdr3) %>% dplyr::summarise_all(
    function(x)mean(x)
  )
}
agg <- plyr_base_fn ()
#agg <- read_delim(file="./Dataset/10X_REF/Donor1/agg.csv", delim = ",")


liste_cdr3_all <- data.frame(cdr3 = donor1_filtered$cdr3, chain = donor1_filtered$chain)
liste_cdr3_all <- liste_cdr3_all %>% filter(chain =="TRB") %>% unique() %>% select(-chain) #26 561 CDR3B

binding_donor1 <- tmp_2

#matrix of score
labelled_cdr3b <- merge(liste_cdr3_all, agg, by="cdr3", all.x =TRUE) %>% remove_first_column()
labelled_cdr3b[is.na(labelled_cdr3b)] <- 0 #replace NA value by 0, for CDR3 non characterized by dextramer in the single cell
labelled_cdr3b <- (labelled_cdr3b > 10) +0 #consider binding when score >10 
labelled_cdr3b <- labelled_cdr3b %>% rownames_as_first_column(., "cdr3")

#write.csv(labelled_cdr3b, paste(directory, "Data/10X/Donor1/labelled_cdr3b.csv", sep=""), quote =F)

liste_cdr3_annotated <- labelled_cdr3b %>% filter(., rowSums(labelled_cdr3b[,2:ncol(labelled_cdr3b)])!=0)# 2876 annotated CDR3B
liste_cdr3_annotated_epitope <- melt(liste_cdr3_annotated) %>% filter(value != 0)
liste_cdr3_annotated <- data.frame("cdr3" = unique(liste_cdr3_annotated$cdr3))

liste_cdr3_annotated_species <- separate(liste_cdr3_annotated_epitope, col = variable, into = c("epitope","protein","species"), sep = "_")
liste_cdr3_annotated_species$species[liste_cdr3_annotated_species$species == "Y"| liste_cdr3_annotated_species$protein == "WT.1"] <- "Cancer"
liste_cdr3_annotated_species$species[liste_cdr3_annotated_species$protein == "PSA146.154"] <- "Human"#self entigen
liste_cdr3_annotated_species$species[liste_cdr3_annotated_species$protein == "Ca2.indepen.Plip.A2"] <- "Human"
liste_cdr3_annotated_species$species[liste_cdr3_annotated_species$protein == "HTLV.1"] <- "HTLV.1"
liste_cdr3_annotated_species$species[liste_cdr3_annotated_species$protein == "Kanamycin.B.dioxygenase"] <- "Streptomyces kanamyceticus"

liste_cdr3_annotated_species$species <- replace_na(liste_cdr3_annotated_species$species, "NC")


#creation of subsets
liste_subset <- create_subset(liste_cdr3_all, liste_cdr3_annotated)

size_subset <- c()
for (j in(1:length(liste_subset))){
  size_subset <- append(size_subset, dim(liste_subset[[j]])[1])
  #myfile <- file.path("../Data/10X/Donor1/subsets/", paste0("subset", "_", j, ".txt"))
  #write.table(liste_subset[[j]], file = myfile, sep = "", row.names = FALSE, col.names = TRUE, quote = FALSE, append = FALSE)
}
subset_cdr3_df <- data.frame(subset_1=size_subset[1], subset_2=size_subset[2], subset_3=size_subset[3], subset_4=size_subset[4], subset_5=size_subset[5], subset_6=size_subset[6]) %>% 
  melt(., value.name="nb_cdr3", variable.name="subset")

#annotations des CDR3
score_binding <- merge(liste_subset[[1]], labelled_cdr3b) %>% remove_first_column()
score_binding <- score_binding[colSums(score_binding)!=0] %>% rownames_as_first_column(., "cdr3")
annotation_cdr3 <- merge(liste_subset[[1]], labelled_cdr3b) %>% select(., cdr3) %>% mutate(, annotation =1)
#write.table(annotation_cdr3,"./Dataset/10X_REF/Donor1/subset/annotated_cdr3.txt", row.names = FALSE, quote = FALSE, sep="\t")
#write.table(score_binding,"./Dataset/10X_REF/Donor1/subset/score_binding_cdr3.txt", row.names = FALSE, quote = FALSE, sep="\t")


##creation of input files
#for HD, LD, TCRMatch, clusTCR, GIANA and iSMART: simple list of CDR3 (without TRBV)

list_input_gliph <- list()
for (i in 1:length(liste_subset)){
  input_gliph <- merge(liste_subset[[i]], donor1_filtered, by ="cdr3", all.x=T) %>% select(., c(cdr3, v_gene, j_gene))
  input_gliph <- input_gliph %>% mutate(., CDR3a = NA) %>% mutate (., Subject_Condition= "Donor1:healthy") %>%
    plyr::count(., vars = c("cdr3", "v_gene", "j_gene","CDR3a", "Subject_Condition"))
  colnames(input_gliph) <- c("CDR3b", "TRBV", "TRBJ", "CDR3a",  "Subject_Condition", "count")
  list_input_gliph[[i]] <- input_gliph
  #myfile <- file.path("../Data/10X/Donor1/input/Gliph2/", paste0("subset", "_", i, ".csv"))
  #write.table(input_gliph, file = myfile, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, append = FALSE)
}

list_input_deeptcr <- list()
for (i in 1:length(liste_subset)){
  input_deeptcr <- merge(liste_subset[[i]], donor1_filtered, by ="cdr3", all.x=T) %>% select(., c(cdr3, v_gene, j_gene)) %>%
    plyr::count(., vars = c("cdr3","v_gene","j_gene"))
  colnames(input_deeptcr) <- c("CDR3_beta", "V_beta", "J_beta", "count")
  list_input_deeptcr[[i]] <- input_deeptcr
  #myfile <- file.path("../Data/10X/Donor1/input/DeepTCR/", paste0("subset", "_", i, ".tsv"))
  #write.table(input_deeptcr, file = myfile, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, append = FALSE)
}

list_input_tcrdist <- list()
for (i in 1:length(liste_subset)){
  input_tcrdist <- merge(liste_subset[[i]], donor1_filtered, by ="cdr3", all.x=TRUE) %>% select(., c(cdr3, v_gene, j_gene)) %>%
    plyr::count(., vars = c("cdr3","v_gene","j_gene"))
  colnames(input_tcrdist) <- c("CDR3b", "TRBV", "TRBJ", "count")
  input_tcrdist$TRBV <- paste(input_tcrdist$TRBV, "*01", sep = "")
  list_input_tcrdist[[i]] <- input_tcrdist
  #myfile <- file.path("../Data/10X/Donor1/input/TCRdist3/", paste0("subset", "_", i, ".tsv"))
  #write.table(input_tcrdist, file = myfile, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, append = FALSE)
}

#train data for DeepTCR
directory_database <- "../Data"
all_database_human <- fread(paste0(directory_database, "/Database/database_pooled_human_2023_03_15.txt"))

data_human <- all_database_human %>% filter(Verified_score ==2 & Identification_score >4.3 & !is.na(Epitope))
data_train <- data_human %>% filter(PubMed_ID != "https://pages.10xgenomics.com/rs/446-PBO-704/images/10x_AN047_IP_A_New_Way_of_Exploring_Immunity_Digital.pdf" & Cell_subset =="CD8")
data_train <- data_train %>% filter(!is.na(CDR3_alpha) & !is.na(CDR3_beta) & !is.na(V_alpha) & !is.na(V_beta) & !is.na(J_alpha) & !is.na(J_beta))

train_deeptcr <- data.frame(CDR3_beta = data_train$CDR3_beta,
                            V_beta = data_train$V_beta,
                            J_beta = data_train$J_beta,
                            Epitope = data_train$Epitope) %>%
  group_by(CDR3_beta, V_beta, J_beta) %>%
  mutate(count = n()) %>% unique()
 
input_deeptcr$clonotype <- paste(input_deeptcr$V_beta,input_deeptcr$CDR3_beta, input_deeptcr$J_beta, sep = "_")
train_deeptcr$clonotype <- paste(train_deeptcr$V_beta,train_deeptcr$CDR3_beta, train_deeptcr$J_beta, sep = "_")

#remove common sequences between train and test sets from train set
common_seq <- intersect(input_deeptcr$clonotype, train_deeptcr$clonotype)#82 clonotypes in common
train_deeptcr <- subset(train_deeptcr,!(clonotype %in% common_seq))

#write.table(train_deeptcr, file = "../Data/10X/TrainData_DeepTCR.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, append = FALSE)


#nb of TCR per epitope
labelled_cdr3b <- read.csv(paste(directory, "Data/10X/Donor1/labelled_cdr3b.csv", sep=""))
labelled_cdr3b <- labelled_cdr3b[-1]

labelled_cdr3b_list <- labelled_cdr3b %>% melt() %>% filter(value!=0)#2876 annotated CDR3b

count_CDR3bAg_human <- labelled_cdr3b_list %>% 
  group_by(variable) %>% 
  summarise(nb_tcr = n_distinct(cdr3))


count_CDR3bAg_human <- count_CDR3bAg_human %>% group_by(nb_tcr) %>% mutate(TCR_epitope_nb_group = n())
count_CDR3bAg_human <- count_CDR3bAg_human %>% select(nb_tcr, TCR_epitope_nb_group) %>%unique() %>%  mutate(nb_epitope = 45)
count_CDR3bAg_human <- count_CDR3bAg_human %>% mutate(TCR_nb_pc = TCR_epitope_nb_group/nb_epitope *100)

nb_epitope_TCR_human <-ggplot(count_CDR3bAg_human , aes(x= nb_tcr, y = TCR_nb_pc))+
  geom_bar(stat="identity", color="black", fill="lightblue", bins=10)+
  xlim(0,11)+
  theme_light()+
  labs(y= "Percentage of Epitope (%)", x = "Number of TCRs") +
  theme(plot.title = element_text(hjust = 0.00001, size = 25, face = "bold"), 
        plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
        axis.text=element_text(size=texte_size),
        text=element_text(size=texte_size))
