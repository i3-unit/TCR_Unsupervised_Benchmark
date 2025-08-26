#Description of human/mouse pooled databases

library(ggplot2)
library(dplyr)
library(gridExtra)
library(cowplot)
library(ggforce)
library(readr)
library(tidyr)
library(ggvenn)

directory <- "../Data"
all_database_human <- read_delim(paste0(directory, "/Database/database_pooled_human_2023_03_15.txt"), delim = "\t")
all_database_human$Cell_subset[all_database_human$Cell_subset=="CD4,CD8"] <- "CD4, CD8"

texte_size <- 18

#Figure 1B: VennDiagram of overlap of studies between the three public databases
db_iedb <- all_database_human %>% filter(Database == "IEDB") %>% drop_na(PubMed_ID)
db_vdjdb <- all_database_human %>% filter(Database == "VDJdb") %>% drop_na(PubMed_ID)
db_mcpas <- all_database_human %>% filter(Database == "Mc-PAS") %>% drop_na(PubMed_ID)

list_db_study <-list(IEDB= db_iedb$PubMed_ID,
                     VDJdb = db_vdjdb$PubMed_ID,
                     McPAS = db_mcpas$PubMed_ID)

venn_diagram <- ggvenn(
  list_db_study, 
  fill_color = c("#0073C2FF", "#868686FF","#EFC000FF"),
  stroke_size = 0, set_name_size = 6, stroke_color = "white", stroke_alpha=0,
  show_elements=F, show_percentage=F, text_size =5) +ggtitle("B") + 
  theme(plot.title = element_text(hjust= 0.01, size = 25, face = "bold"),
        plot.margin=unit(c(0,1,0,1),"cm"))


#Figure 1C: participation of each database in the pooled database
database_perc_human <- all_database_human %>% 
  group_by(Database) %>% 
  summarise(count_db = n(), .groups = 'drop') %>% 
  arrange(desc(Database)) %>%
  mutate(perc_db = count_db/sum(count_db)*100) %>% 
  mutate(db_ypos = cumsum(perc_db) - 0.5*perc_db) %>% 
  mutate(species = "Human")

mycols <- c("#0073C2FF", "#EFC000FF", "#868686FF")


pieChart_db_human <- ggplot(database_perc_human, aes(x = "", y = perc_db, fill = Database)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = db_ypos, label = paste0(round(perc_db, 0), "%")), color = "black", size = 8)+
  scale_fill_manual(values = mycols) +
  theme_void()+
  ggtitle("C")+
  theme(legend.position = "none",
        plot.title = element_text(hjust= 0.01, size = 25, face = "bold"),
        plot.margin=unit(c(0,1,0,1),"cm"))


all_database_human <- all_database_human %>% 
  mutate(chain_type = ifelse(!is.na(CDR3_alpha) & !is.na(CDR3_beta), "pairing",#split by chain type
                             ifelse(!is.na(CDR3_alpha) & is.na(CDR3_beta), "alpha", 
                                    ifelse(is.na(CDR3_alpha) & !is.na(CDR3_beta), "beta", NA))))

typechain_perc_human <- all_database_human %>% 
  group_by(chain_type, Database) %>% 
  summarise(count = n(), .groups = 'drop') %>% 
  mutate(perc = count/sum(count)*100) %>% 
  mutate(species = "Human")


barplot_database_human <- ggplot(data=typechain_perc_human, aes(x=chain_type, y=perc, fill=Database)) +
  geom_bar(stat="identity", color="white")+
  scale_fill_manual(values = mycols)+
  theme_bw()+
  ggtitle("")+
  ylab("Percentage (%)") +
  xlab("")+
  theme(legend.position="right",
        plot.margin=unit(c(1,0,1,0),"cm"),
        legend.text=element_text(size=texte_size),
        legend.title=element_text(size=texte_size),
        text=element_text(size=texte_size),
        axis.text=element_text(size=texte_size),
        axis.text.x = element_text(angle = 45, hjust = 1))


#Figure 1D: cell type
typecell_perc_human <- all_database_human %>% 
  group_by(Cell_subset) %>% 
  summarise(count = n(), .groups = 'drop') %>% 
  arrange(desc(Cell_subset)) %>%
  mutate(perc = count/sum(count)*100) %>% 
  mutate(ypos = cumsum(perc) - 0.5*perc) %>% 
  mutate(species = "Human")

piechart_cellsubset_human <- ggplot(typecell_perc_human, aes(x = "", y = perc, fill = Cell_subset)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = ypos, label = ifelse(perc > 1, paste0(round(perc, 0), "%"), NA)), color = "black", size = 8)+
  scale_fill_brewer(palette="Dark2", na.value="grey")+
  theme_void()+
  ggtitle("D")+
  theme(legend.position = "none",
        plot.title = element_text(hjust= 0.01, size = 25, face ="bold"),
        strip.text = element_text(size = 15),
        plot.margin=unit(c(0,1,0,1),"cm"))

typecellchain_perc_human <- all_database_human %>% 
  group_by(chain_type, Cell_subset) %>% 
  summarise(count = n(), .groups = 'drop') %>% 
  mutate(perc = count/sum(count)*100) %>% 
  mutate(species = "Human")

barplot_cellsubset_human <- ggplot(data=typecellchain_perc_human, aes(x=chain_type, y=perc, fill=Cell_subset)) +
  geom_bar(stat="identity", color="white")+
  scale_fill_brewer(palette="Dark2", na.value="grey")+
  theme_bw()+
 ggtitle(" ")+
  ylab("Percentage (%)") +
  xlab("")+
  theme(legend.position="right",plot.margin=unit(c(1,0,1,0),"cm"),
        legend.text=element_text(size=texte_size),
        legend.title=element_text(size=texte_size),
        text=element_text(size=texte_size),
        axis.text=element_text(size=texte_size),
        axis.text.x = element_text(angle = 45, hjust = 1))

#Figure 1H: Verified score
# count_verified_score_h <- all_database_human %>% 
#   group_by(Verified_score) %>% 
#   summarise(count = n(), .groups = 'drop') %>% 
#   mutate(perc = count/sum(count)*100) %>% 
#   mutate(species = "Human")

# count_verified_score_h$Verified_score <- as.character(count_verified_score_h$Verified_score)

# barplot_verifiedscore_all_human <- ggplot(data=count_verified_score_h, aes(x=species, y=perc, fill=Verified_score)) +
#   geom_bar(stat="identity", color="white")+
#   scale_fill_brewer(palette="Pastel2", na.value="grey")+
#   theme_bw()+
#   ggtitle("E")+
#   ylab("Percentage (%)")+
#   xlab("Pooled database")+
#   theme(plot.title = element_text(hjust = 0.00001, size=25, face = "bold"), 
#         legend.position = "none",plot.margin=unit(c(1,1,1,1),"cm"))

verifiedscore_perc_human <- all_database_human %>% 
  group_by(Verified_score) %>% 
  summarise(count = n(), .groups = 'drop') %>% 
  arrange(desc(Verified_score)) %>%
  mutate(perc = count/sum(count)*100) %>% 
  mutate(ypos = cumsum(perc) - 0.5*perc) %>% 
  mutate(species = "Human")

verifiedscore_perc_human$Verified_score <- as.character(verifiedscore_perc_human$Verified_score)

piechart_verifiedscore_all_human <- ggplot(data=verifiedscore_perc_human, aes(x = "", y = perc, fill = Verified_score)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = ypos, label = ifelse(perc > 1, paste0(round(perc, 0), "%"), NA)), color = "black", size = 8)+
  scale_fill_brewer(palette="Pastel2", na.value="grey")+
  theme_void()+
  ggtitle("E")+
  theme(legend.position = "none",
        plot.title = element_text(hjust= 0.01, size = 25, face = "bold"),
        plot.margin=unit(c(0,1,0,1),"cm"))


verified_score_perc_human <- all_database_human %>% 
  group_by(chain_type, Verified_score) %>% 
  summarise(count = n(), .groups = 'drop') %>% 
  mutate(perc = count/sum(count)*100) %>% 
  mutate(species = "Human")

verified_score_perc_human$Verified_score <- as.character(verified_score_perc_human$Verified_score)

barplot_verifiedscore_split_human <- ggplot(data=verified_score_perc_human, aes(x=chain_type, y=perc, fill=Verified_score)) +
  geom_bar(stat="identity", color="white")+
  scale_fill_brewer(palette="Pastel2", na.value="grey")+
  theme_bw()+
  ggtitle(" ")+
  ylab("Percentage (%)")+
  xlab("")+
  theme(legend.position="right",
        plot.margin=unit(c(1,0,1,0),"cm"),
        legend.text=element_text(size=texte_size),
        legend.title=element_text(size=texte_size),
        text=element_text(size=texte_size),
        axis.text=element_text(size=texte_size),
        axis.text.x = element_text(angle = 45, hjust = 1))


#Figure 1I: Ag identification score
all_database_human <- all_database_human %>% mutate(interval_score = ifelse(Identification_score ==3.1 | Identification_score ==4.1, "3.1,4.1",
                                                                            ifelse(Identification_score ==3.2 | Identification_score ==4.2, "3.2,4.2",
                                                                                   ifelse(Identification_score ==3.3 | Identification_score ==4.3, "3.3,4.3",
                                                                                          ifelse(Identification_score ==3.4 | Identification_score ==4.4, "3.4,4.4",
                                                                                                 ifelse(Identification_score ==3.5 | Identification_score ==4.5, "3.5,4.5",Identification_score))))))


# count_ag_score_h <- all_database_human %>% 
#   group_by(interval_score) %>% 
#   summarise(count = n(), .groups = 'drop') %>% 
#   mutate(perc = count/sum(count)*100) %>% 
#   mutate(species = "Human")
# 
# count_ag_score_h$interval_score <- as.character(count_ag_score_h$interval_score)

# barplot_agscore_all_human <- ggplot(data=count_ag_score_h, aes(x=species, y=perc, fill=interval_score)) +
#   geom_bar(stat="identity", color="white",alpha =0.5)+
#   scale_fill_brewer(palette = "BrBG") +
#   theme_bw()+
#   ggtitle("F")+
#   ylab("Percentage (%)")+
#   xlab("Pooled database")+
#   theme(plot.title = element_text(hjust = 0.001, size =25, face="bold"), 
#         legend.position = "none",
#         plot.margin=unit(c(1,1,1,1),"cm"))

agscore_perc_human <- all_database_human %>% 
  group_by(interval_score) %>% 
  summarise(count = n(), .groups = 'drop') %>% 
  arrange(desc(interval_score)) %>%
  mutate(perc = count/sum(count)*100) %>% 
  mutate(ypos = cumsum(perc) - 0.5*perc) %>% 
  mutate(species = "Human")


piechart_agscore_all_human <- ggplot(data=agscore_perc_human, aes(x = "", y = perc, fill = interval_score)) +
  geom_bar(width = 1, stat = "identity", color = "white", alpha =1) +
  coord_polar("y", start = 0)+
  geom_text(aes(y = ypos, label = ifelse(perc > 2, paste0(round(perc,0), "%"), NA)), color = "black", size =8)+
  scale_fill_brewer(palette = "Paired")+
  theme_void()+
  ggtitle("F")+
  theme(legend.position = "none",
        plot.title = element_text(hjust= 0.01, size = 25, face = "bold"),
        plot.margin=unit(c(0,1,0,1),"cm"))

count_ag_score_h <- all_database_human %>% 
  group_by(chain_type, interval_score) %>% 
  summarise(count = n(), .groups = 'drop') %>% 
  mutate(perc = count/sum(count)*100) %>% 
  mutate(species = "Human")

count_ag_score_h$interval_score <- as.character(count_ag_score_h$interval_score)

barplot_agscore_split_human <- ggplot(data=count_ag_score_h, aes(x=chain_type, y=perc, fill=interval_score)) +
  geom_bar(stat="identity", color="white", alpha = 1)+
  scale_fill_brewer(palette="Paired", na.value="grey")+
  theme_bw()+
  theme(plot.title = element_text(hjust= 0.5, size = 16))+
  ggtitle(" ")+
  ylab("Percentage (%)")+
  xlab("")+
  theme(legend.position="right",
        plot.margin=unit(c(1,0,1,0),"cm"),
        legend.text=element_text(size=texte_size),
        legend.title=element_text(size=texte_size),
        text=element_text(size=texte_size),
        axis.text=element_text(size=texte_size),
        axis.text.x = element_text(angle = 45, hjust = 1))


#Fgure 1E: antigen organism
count_organism_h <- all_database_human %>% 
  group_by(Antigen_organism) %>% 
  summarise(count = n(), .groups = 'drop') %>% 
  mutate(perc = count/sum(count)*100) %>% 
  mutate(species = "Human") 

count_organism_h_biggest <- count_organism_h %>% filter(perc >=1) %>% filter(!is.na(Antigen_organism))

count_organism_h_biggest$Antigen_organism[count_organism_h_biggest$Antigen_organism=="SARS coronavirus Tor2 (Severe acute respiratory syndrome-related coronavirus Tor2)"] <- "SARS coronavirus Tor2"

list_organism_color <- list("SARS-CoV2"="#FFFF99",
                            "Cytomegalovirus (CMV)"= "#7FC97F",
                            "Epstein Barr virus (EBV)"= "#BEAED4",
                            "Influenza A virus" = "#FDC086",
                            "Mycobacterium tuberculosis" = "#F0027F",
                            "Yellow fever virus (YFV)" = "#386CB0",
                            "Homo sapiens (human)"= "#666666",
                            "SARS-CoV1"= "darkcyan",
                            "Hepatitis B virus (HBV)"= "#BF5B17",
                            "SARS coronavirus Tor2" = "aquamarine2")

count_organism_human <- ggplot(data=count_organism_h_biggest, aes(x= reorder(Antigen_organism, perc), y=perc, fill = Antigen_organism)) +
  geom_bar(stat="identity", alpha = 1)+
  geom_text(aes(label=round(perc,1)), hjust=-0.1, size=4)+
  theme_classic()+
  scale_fill_manual(values = list_organism_color)+
  theme(plot.title = element_text(hjust= 0.001, size = 25, face = "bold"),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), 
        legend.position = "none",
        axis.text=element_text(size=texte_size),
        text=element_text(size=texte_size))+
  ggtitle("G")+
  ylab("Percentage of TCRs (%)")+
  xlab("")+
  coord_flip()


#Figure 1F: antigen epitope
count_epitope_h <- all_database_human %>% 
  group_by(Epitope, Antigen_organism) %>% 
  summarise(count = n(), .groups = 'drop') %>% 
  mutate(perc = count/sum(count)*100) %>% 
  mutate(species = "Human") 


count_epitope_h_biggest <- count_epitope_h %>% filter(perc >=1) %>% filter(!is.na(Epitope))

count_epitope_human <- ggplot(data=count_epitope_h_biggest, aes(x= reorder(Epitope, perc), y=perc, fill = Antigen_organism)) +
  geom_bar(stat="identity", alpha = 1)+
  geom_text(aes(label=round(perc,1)), hjust=-0.1, size=4)+
  theme_classic()+
  scale_fill_brewer(palette="Accent", na.value="grey")+
  theme(plot.title = element_text(hjust= 0.001, size = 25, face = "bold"),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
        legend.position = "bottom",
        legend.title=element_text(size=texte_size), 
        legend.text=element_text(size=texte_size),
        axis.text=element_text(size=texte_size),
        text=element_text(size=texte_size))+
  ggtitle("H")+
  ylab("Percentage of TCRs (%)")+
  xlab("")+
  coord_flip()
legend <- get_legend(count_epitope_human)
count_epitope_human <- count_epitope_human + theme(legend.position="none")


#Figure 1G: Number of TCRs (alpha/beta) which recognise one Ag
count_TCRAg_human <- all_database_human %>% 
  mutate(pair = paste0(CDR3_alpha, '_',CDR3_beta)) %>% 
  group_by(Epitope) %>% 
  summarise(nb_tcr = n_distinct(pair))


count_TCRAg_human <- count_TCRAg_human %>% group_by(nb_tcr) %>% mutate(TCR_epitope_nb_group = n())
count_TCRAg_human <- count_TCRAg_human %>% select(nb_tcr, TCR_epitope_nb_group) %>%unique() %>%  mutate(nb_epitope = 2314)
count_TCRAg_human <- count_TCRAg_human %>% mutate(TCR_nb_pc = TCR_epitope_nb_group/nb_epitope *100)

nb_epitope_TCR_human <-ggplot(count_TCRAg_human, aes(x= nb_tcr, y = TCR_nb_pc))+
  geom_bar(stat="identity", color="black", fill="lightblue", bins=10)+
  xlim(0,11)+
  theme_light()+
  labs(y= "Percentage of Epitope (%)", x = "Number of TCRs") +
  ggtitle("G")+
  theme(plot.title = element_text(hjust = 0.00001, size = 25, face = "bold"), 
        plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
        axis.text=element_text(size=texte_size),
        text=element_text(size=texte_size))


###############################################################################$######
#Figure 1 C-F
grid.arrange(pieChart_db_human, barplot_database_human, piechart_cellsubset_human, barplot_cellsubset_human,
             piechart_verifiedscore_all_human, barplot_verifiedscore_split_human, piechart_agscore_all_human, barplot_agscore_split_human,
             widths = c(1, 1, 1,1),heights=c(1,1),
             layout_matrix = rbind(c(1, 2, 3, 4),
                                   c(5, 6, 7, 8)))

#Figure 1 G-I
grid.arrange(count_organism_human, count_epitope_human, nb_epitope_TCR_human, legend,
             widths = c(1,1,1),heights=c(1, 0.1),
             layout_matrix = rbind(c(1, 2, 3),
                                   c(4,4, NA)))

grid.arrange(count_organism_human, count_epitope_human, legend,
             widths = c(1,1),heights=c(1, 0.1),
             layout_matrix = rbind(c(1, 2),
                                   c(4,4)))
