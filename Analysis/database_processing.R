#Pooled database processing
#create input for methods
library(data.table)
library(dplyr)
library(tidyr)
source("../Functions/general_functions.R")

directory_database <- "../Data"
all_database_human <- fread(paste0(directory_database, "/Database/database_pooled_human_2023_03_15.txt"))

#Database filtering
data_human <- all_database_human %>% 
  filter(Verified_score ==2 & Identification_score >4.3 & !is.na(Epitope)) %>% 
  filter(PubMed_ID != "https://pages.10xgenomics.com/rs/446-PBO-704/images/10x_AN047_IP_A_New_Way_of_Exploring_Immunity_Digital.pdf") %>% 
  filter(Cell_subset == "CD8") %>% 
  drop_na(V_alpha) %>% drop_na(J_alpha) %>% drop_na(V_beta) %>% drop_na(J_beta)

data_human <- subset(data_human, nchar(as.character(data_human$CDR3_alpha))<=23 & nchar(as.character(data_human$CDR3_alpha))>=6 & nchar(as.character(data_human$CDR3_beta))<=23 & nchar(as.character(data_human$CDR3_beta))>=6)

#check if at least 2 TCR bind each epitope
data_human <- data_human %>%
  mutate(pair = paste0(CDR3_beta, '_',CDR3_alpha)) %>%
  group_by(Epitope) %>%
  mutate(nb_tcr = n_distinct(pair)) %>%
  filter(nb_tcr != 1) %>% unique() %>% mutate(origin = "data_human")

data_human <- data_human %>% filter(pair != "CASSDSRGTEAFF_CASSDSRGTEAFF") %>% filter(pair != "CASSPITGTGAYGYTF_CASSSVNEQYF") %>% unique()#remove doublons
data_human <- data_human %>% mutate(CDR3b_CDR3a= paste0(CDR3_beta, '',CDR3_alpha))

#separate and concatenate pairing alpha/beta
dataframe_beta <- data_human %>% select(CDR3_beta, V_beta, J_beta, Epitope, Antigen, Antigen_organism, Cell_subset, PubMed_ID, Sequencing_method, Antigen_identification, Antigen_type, Verified_score, Database, Identification_score) %>% as.data.frame() %>% mutate(chain = "beta")
dataframe_beta <- dataframe_beta %>% dplyr::rename("CDR3" ="CDR3_beta",
                                                   "V"= "V_beta",
                                                   "J" = "J_beta")

dataframe_alpha <- data_human %>% select(CDR3_alpha, V_alpha, J_alpha, Epitope, Antigen, Antigen_organism, Cell_subset, PubMed_ID, Sequencing_method, Antigen_identification, Antigen_type, Verified_score, Database, Identification_score) %>% mutate(chain = "alpha")
dataframe_alpha <- dataframe_alpha %>% dplyr::rename( "CDR3" = "CDR3_alpha",
                                                      "V" = "V_alpha",
                                                      "J" = "J_alpha")
dataframe_alpha_beta <- rbind(dataframe_alpha, dataframe_beta)#concatenate alpha beta
(rm(dataframe_beta, dataframe_alpha))

#unique sequences
liste_sequence_TRA_TRB <- data.frame(CDR3=unique(dataframe_alpha_beta$CDR3))
liste_sequence_TRA_TRB_chain <- data.frame(CDR3 = dataframe_alpha_beta$CDR3, chain = dataframe_alpha_beta$chain) %>% unique()
liste_sequence_TRA_TRB_chain <- liste_sequence_TRA_TRB_chain %>% #add column with color for the chain
  mutate(color_type = case_when(
    chain == "alpha" ~ "chartreuse3",
    chain == "beta" ~ "darkorange1"))


liste_sequence_TRA <- liste_sequence_TRA_TRB_chain %>% filter(chain == "alpha")
liste_sequence_TRB <- liste_sequence_TRA_TRB_chain %>% filter(chain == "beta")

#Create input for methods
# input_gliph <- data.frame(CDR3_beta = data_human$CDR3_beta, 
#                           V_beta = data_human$V_beta, 
#                           J_beta = data_human$J_beta, 
#                           CDR3_alpha = data_human$CDR3_alpha) %>% 
#   mutate(Subject_Condition = "database:CD8+") %>%
#   group_by(CDR3_beta, V_beta, J_beta, CDR3_alpha) %>% 
#   mutate(freq = n()) %>% unique()
# 
# 
# input_deeptcr <- data.frame(CDR3_alpha = data_human$CDR3_alpha,
#                             CDR3_beta = data_human$CDR3_beta, 
#                             V_alpha = data_human$V_alpha,
#                             V_beta = data_human$V_beta, 
#                             J_alpha = data_human$J_alpha, 
#                             J_beta = data_human$J_beta) %>% 
#   group_by(CDR3_beta, V_beta, J_beta, CDR3_alpha, V_alpha, J_alpha) %>% 
#   mutate(count = n()) %>% unique()
# 
# 
# input_clustcr <- data.frame(CDR3b = data_human$CDR3_beta, 
#                             CDR3a = data_human$CDR3_alpha) %>% unique()
# 
# 
# input_tcrdist <- data.frame(CDR3a = data_human$CDR3_alpha,
#                             CDR3b = data_human$CDR3_beta, 
#                             TRAV = data_human$V_alpha,
#                             TRBV = data_human$V_beta, 
#                             TRAJ = data_human$J_alpha, 
#                             TRBJ = data_human$J_beta) %>% 
#   group_by(CDR3b, TRBV, TRBJ, CDR3a, TRAV, TRAJ) %>% 
#   mutate(count = n()) %>% unique()
# 
# input_tcrdist$TRAV <- paste(input_tcrdist$TRAV, "*01", sep = "")
# input_tcrdist$TRBV <- paste(input_tcrdist$TRBV, "*01", sep = "")


#DeepTCR - train data
# data_human <- all_database_human %>% filter(Verified_score ==2 & Identification_score >4.3 & !is.na(Epitope))
# data_train <- data_human %>% filter(PubMed_ID == "https://pages.10xgenomics.com/rs/446-PBO-704/images/10x_AN047_IP_A_New_Way_of_Exploring_Immunity_Digital.pdf" & Cell_subset =="CD8")
# data_train <- data_train %>% filter(!is.na(CDR3_alpha) & !is.na(CDR3_beta) & !is.na(V_alpha) & !is.na(V_beta) & !is.na(J_alpha) & !is.na(J_beta))
# 
# train_deeptcr <- data.frame(CDR3_alpha = data_train$CDR3_alpha,
#                             CDR3_beta = data_train$CDR3_beta, 
#                             V_alpha = data_train$V_alpha,
#                             V_beta = data_train$V_beta, 
#                             J_alpha = data_train$J_alpha, 
#                             J_beta = data_train$J_beta, 
#                             Epitope = data_train$Epitope) %>% 
#   group_by(CDR3_beta, V_beta, J_beta, CDR3_alpha, V_alpha, J_alpha) %>% 
#   mutate(count = n()) %>% unique()
# 
# input_deeptcr$clonotype <- paste(input_deeptcr$V_alpha, input_deeptcr$CDR3_alpha, input_deeptcr$J_alpha , input_deeptcr$V_beta,input_deeptcr$CDR3_beta, input_deeptcr$J_beta, sep = "_")
# train_deeptcr$clonotype <- paste(train_deeptcr$V_alpha, train_deeptcr$CDR3_alpha, train_deeptcr$J_alpha , train_deeptcr$V_beta,train_deeptcr$CDR3_beta, train_deeptcr$J_beta, sep = "_")
# 
# #remove common sequences between train and test sets from train set
# common_seq <- intersect(input_deeptcr$clonotype, train_deeptcr$clonotype)
# train_deeptcr <- subset(train_deeptcr,!(clonotype %in% common_seq))


#unique pairs
unique_pairing <- data.frame(CDR3_beta = data_human$CDR3_beta, CDR3_alpha = data_human$CDR3_alpha) %>% unique()
unique_pairing$CDR3b_CDR3a <- paste(unique_pairing$CDR3_beta, unique_pairing$CDR3_alpha, sep = "")
unique_pairing <- unique_pairing %>% select(-c(CDR3_beta, CDR3_alpha))

unique_clone <- data_human %>% select(CDR3_beta, CDR3_alpha, V_alpha, J_alpha, V_beta, J_beta) %>%  unique()
unique_clone$V_alpha <- paste(unique_clone$V_alpha, "*01", sep ="")
unique_clone$V_beta <- paste(unique_clone$V_beta, "*01", sep ="")
unique_clone$clone <- paste(unique_clone$V_alpha, unique_clone$CDR3_alpha, unique_clone$J_alpha, unique_clone$V_beta, unique_clone$CDR3_beta, unique_clone$J_beta, sep = "_")
unique_clone <- data.frame(clone=unique_clone$clone) %>% unique()


#create a specificity matrix
list_epitope <- c("GILGFVFTL","NLVPMVATV","GLCTLVAML","FLCMKALLL", "NYNYLYRLF", "LLFGYPVYV")

matrix_alpha_beta <- dataframe_alpha_beta %>% select(CDR3, Epitope) %>% unique() %>% mutate(bind = 1) %>% filter(Epitope %in% list_epitope) %>% reshape2::dcast(CDR3~Epitope)
matrix_alpha_beta[is.na(matrix_alpha_beta)] <- 0
matrix_alpha_beta <- matrix_alpha_beta[,append("CDR3", list_epitope)]

matrix_alpha_beta_other <- dataframe_alpha_beta %>% select(CDR3, Epitope) %>% unique() %>% mutate(bind = 1) %>%filter(Epitope %!in% list_epitope) %>% reshape2::dcast(CDR3~Epitope)
matrix_alpha_beta_other[is.na(matrix_alpha_beta_other)] <- 0
matrix_alpha_beta_other <- matrix_alpha_beta_other %>% mutate(other_specificities = rowSums(matrix_alpha_beta_other[-1]))

matrix_alpha_beta_final <- merge(matrix_alpha_beta, matrix_alpha_beta_other[,c("CDR3", "other_specificities")], by = "CDR3", all = T)
matrix_alpha_beta_final[is.na(matrix_alpha_beta_final)] <- 0
matrix_alpha_beta_final <- matrix_alpha_beta_final %>% mutate(polyspe = ifelse(rowSums(matrix_alpha_beta_final[-1])>1, "poly", "mono"))


matrix_pair <- data_human %>% select(pair, Epitope) %>% unique() %>% mutate(bind = 1) 
matrix_pair$pair <- gsub(x=matrix_pair$pair, pattern = "_", replacement ="")
matrix_pair <- matrix_pair %>% filter(Epitope %in% list_epitope) %>% reshape2::dcast(pair~Epitope)
matrix_pair[is.na(matrix_pair)] <- 0
matrix_pair <- matrix_pair[,append("pair", list_epitope)]

matrix_pair_other <- data_human %>% select(pair, Epitope) %>% unique() %>% mutate(bind = 1) 
matrix_pair_other$pair <- gsub(x=matrix_pair_other$pair, pattern = "_", replacement ="")
matrix_pair_other <- matrix_pair_other %>% filter(Epitope %!in% list_epitope) %>% reshape2::dcast(pair~Epitope)
matrix_pair_other[is.na(matrix_pair_other)] <- 0
matrix_pair_other <- matrix_pair_other %>% mutate(other_specificities = rowSums(matrix_pair_other[-1]))

matrix_pair_final <- merge(matrix_pair, matrix_pair_other[,c("pair", "other_specificities")], by = "pair", all = T)
matrix_pair_final[is.na(matrix_pair_final)] <- 0
matrix_pair_final <- matrix_pair_final %>% mutate(polyspe = ifelse(rowSums(matrix_pair_final[-1])>1, "poly", "mono"))

rm(matrix_alpha_beta, matrix_alpha_beta_other, matrix_pair, matrix_pair_other)

epitope_organism <- dataframe_alpha_beta[,c("Epitope", "Antigen_organism")] %>% unique() %>% drop_na()

