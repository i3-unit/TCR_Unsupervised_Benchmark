# Processing of clustering outputs
library(dplyr)
library(readr)
library(stringdist)
library(reshape2)
library(igraph)
library(stringr)


directory <- "D:/Labo_i3/Paper_Benchmark_v2/Rcode_Github/TCR_Unsupervised_Benchmark_modif/"
source(paste(directory, "/Functions/general_functions.R", sep=""))

liste_subset <- list()
for (i in 1:6) {
  print(i)
  myfile <- file.path(paste0(directory,"Data/10X/Donor1/subsets"), paste0("subset", "_", i, ".txt"))
  liste_subset[[i]] <- read_delim(file=myfile, delim = ",") %>% data.frame()
  i <- i+1
}


## Analysis of outputs

#1. Hamming output
version_list <- c("a","b","c")
output_list_HD <- list()
clustered_cdr3_HD <- data.frame()
for (i in 1:6) {
  print(i)
  for(k in version_list){
    myfile <- file.path(paste0(directory,"Data/10X_Output_methods/output_HD/"), paste0("subset", "_", i,k, ".csv"))
    output_list_HD[[paste("subset_",as.character(i),k,sep="")]] <- read_delim(file=myfile, delim = "\t") %>% data.frame()
    
    df <- data.frame(version = k, subset = paste("subset",i,sep="_"), nb_cdr3 = n_distinct(output_list_HD[[paste("subset_",as.character(i),k,sep="")]]$CDR3))
    clustered_cdr3_HD <- rbind(clustered_cdr3_HD, df)
  }
}



#2. Levenshtein output
output_list_LD <- list()
clustered_cdr3_LD <- data.frame()
for (i in 1:6) {
  print(i)
  for(k in version_list){
    myfile <- file.path(paste0(directory,"Data/10X_Output_methods/output_LD/"), paste0("subset", "_", i,k, ".csv"))
    output_list_LD[[paste("subset_",as.character(i),k,sep="")]] <- read_delim(file=myfile, delim = "\t") %>% data.frame()
    
    df <- data.frame(version = k, subset = paste("subset",i,sep="_"), nb_cdr3 = n_distinct(output_list_LD[[paste("subset_",as.character(i),k,sep="")]]$CDR3))
    clustered_cdr3_LD <- rbind(clustered_cdr3_LD, df)
  }
}


#3. TCRMatch
output_list_tcrm <- list()
clustered_cdr3_tcrm <- data.frame()
for (i in 1:6) {
  print(i)
  for(k in version_list){
    myfile <- file.path(paste0(directory,"Data/10X_Output_methods/output_TCRMatch"), paste0("output_TCRMatch_sample", i,k, ".txt"))
    edge_list_tcrm <- data.table::fread(myfile, header = F) %>% as.data.frame()
    colnames(edge_list_tcrm) <- c("input_sequence", "match_sequence", "score")
    edge_list_tcrm <- edge_list_tcrm %>% filter(score != 1)
    
    graph_tcrm <- graph_from_data_frame(edge_list_tcrm, directed = F)
    component <- igraph::components(graph_tcrm)
    output_tcrm <- data.frame(clusters = component$membership) %>% rownames_as_first_column(., "CDR3") %>% 
      dplyr::group_by(clusters) %>% dplyr::mutate(size = n())
    output_list_tcrm[[paste("subset_",as.character(i),k,sep="")]] <- output_tcrm
    
    df <- data.frame(version = k, subset = paste("subset",i,sep="_"), nb_cdr3 = n_distinct(output_list_tcrm[[paste("subset_",as.character(i),k,sep="")]]$CDR3))
    clustered_cdr3_tcrm <- rbind(clustered_cdr3_tcrm, df)
  }
}


#4. iSMART
output_list_ismart <- list()
clustered_cdr3_ismart <- data.frame()
for (i in 1:6) {
  print(i)
  for(k in version_list){
    myfile <- file.path(paste0(directory,"Data/10X_Output_methods/output_iSMART"), paste0("subset_", i,k,"_results.txt_ClusteredCDR3s_7.5.txt"))
    output_ismart <- read.table(myfile, header=T,sep='\t',stringsAsFactors = F)
    colnames(output_ismart) <- c('CDR3', 'clusters')
    output_ismart <- output_ismart %>% group_by(clusters) %>% dplyr::mutate(size = n())
    output_list_ismart[[paste("subset_",as.character(i),k,sep="")]] <- output_ismart
    
    df <- data.frame(version = k, subset = paste("subset",i,sep="_"), nb_cdr3 = n_distinct(output_list_ismart[[paste("subset_",as.character(i),k,sep="")]]$CDR3))
    clustered_cdr3_ismart <- rbind(clustered_cdr3_ismart, df)
  }
}


#5. GIANA
output_list_giana <- list()
clustered_cdr3_giana <- data.frame()
for (i in 1:6) {
  print(i)
  for(k in version_list){
    myfile <- file.path(paste0(directory,"Data/10X_Output_methods/output_GIANA"), paste0("subset_", i,k,"_results.txt"))
    output_giana <- read.table(myfile, header=F,sep='\t',stringsAsFactors = F)
    colnames(output_giana) <- c('CDR3', 'clusters', 'TRBV')
    output_giana <- output_giana %>% group_by(clusters) %>% dplyr::mutate(size = n()) %>% select(-TRBV)
    output_list_giana[[paste("subset_",as.character(i),k,sep="")]] <- output_giana
    
    df <- data.frame(version = k, subset = paste("subset",i,sep="_"), nb_cdr3 = n_distinct(output_list_giana[[paste("subset_",as.character(i),k,sep="")]]$CDR3))
    clustered_cdr3_giana <- rbind(clustered_cdr3_giana, df)
  }
}


#6. ClusTCR
output_list_clustcr <- list()
clustered_cdr3_clustcr <- data.frame()
for (i in 1:6) {
  print(i)
  for(k in version_list){
    directory_clustcr <- paste0(directory,"Data/10X_Output_methods/output_clusTCR/subset", i, k, "/")
    output_clustcr <- read_delim(file=paste0(directory_clustcr,"subset_",i, "_output.csv"), delim = ",")
    colnames(output_clustcr) <- c("CDR3","clusters")
    output_clustcr_summary <- read_delim(file=paste0(directory_clustcr, "subset_",i, "_output_summary.csv"), delim = ",")
    colnames(output_clustcr_summary)[1] <- "clusters"
  
    output_clustcr_merge <- merge(output_clustcr, output_clustcr_summary[, c("clusters", "size")], by="clusters")
    output_clustcr_merge <- output_clustcr_merge %>% relocate(CDR3, .before=clusters)
    
    output_list_clustcr[[paste("subset_",as.character(i),k,sep="")]] <- output_clustcr_merge
    
    df <- data.frame(version = k, subset = paste("subset",i,sep="_"), nb_cdr3 = n_distinct(output_list_clustcr[[paste("subset_",as.character(i),k,sep="")]]$CDR3))
    clustered_cdr3_clustcr <- rbind(clustered_cdr3_clustcr, df)
  }
}


#7. GLIPH2
output_list_gliph2 <- list()
clustered_cdr3_gliph2 <- data.frame()
for (i in 1:6) {
  print(i)
  for(k in version_list){
    directory_gliph2 <- paste0(directory,"Data/10X_Output_methods/output_gliph2/output_subset", i, k, "/")
    output_gliph2 <- read_delim(file=paste0(directory_gliph2,"output_cluster.csv"), delim = ",")

    output_gliph2 <- output_gliph2[1:18]
    output_gliph2 <- output_gliph2 %>% select(-TcRa)
    output_gliph2 <- na.omit(output_gliph2)
    output_gliph2 <- filter(output_gliph2, pattern != "single") 
    output_gliph2 <- output_gliph2 %>% group_by(index) %>% dplyr::mutate(size = n())
    output_gliph2 <- dplyr::rename(output_gliph2, "CDR3" = "TcRb",
                                  "clusters" = "index")
    output_gliph2 <- output_gliph2 %>% select(CDR3, clusters, size) %>% unique()
    output_list_gliph2[[paste("subset_",as.character(i),k,sep="")]] <- output_gliph2
    
    df <- data.frame(version = k, subset = paste("subset",i,sep="_"), nb_cdr3 = n_distinct(output_list_gliph2[[paste("subset_",as.character(i),k,sep="")]]$CDR3))
    clustered_cdr3_gliph2 <- rbind(clustered_cdr3_gliph2, df)
    
  }
}


#8. DeepTCR
output_list_deeptcr <- list()
clustered_cdr3_deeptcr <- data.frame()
for (i in 1:6) {
  print(i)
  for(k in version_list){
    myfile <- file.path(paste0(directory,"Data/10X_Output_methods/output_DeepTCR"), paste0("results_clusteringTRB_subset", i,k,".csv"))
    output_deeptcr <- read_delim(file=myfile, delim = ",")
    output_deeptcr <- output_deeptcr %>% select(cdr3, cluster,size)
    colnames(output_deeptcr) <- c("CDR3", "clusters", "size")
    
    output_list_deeptcr[[paste("subset_",as.character(i),k,sep="")]] <- output_deeptcr
    
    df <- data.frame(version = k, subset = paste("subset",i,sep="_"), nb_cdr3 = n_distinct(output_list_deeptcr[[paste("subset_",as.character(i),k,sep="")]]$CDR3))
    clustered_cdr3_deeptcr <- rbind(clustered_cdr3_deeptcr, df)
  }
}



#9. TCRdist3
output_list_tcrdist <- list()
clustered_cdr3_tcrdist <- data.frame()
for (i in 1:6) {
  print(i)
  for(k in version_list){
    directory_tcrdist <- paste0(directory,"Data/10X_Output_methods/output_TCRdist3/subset", i, k, "/")
    output_tcrdist <- read_delim(file=paste0(directory_tcrdist,"output_hclust_beta_cdr3.csv"), delim = ";")
    output_tcrdist <- output_tcrdist %>% group_by(clusters) %>% dplyr::mutate(size= n()) %>% filter(size != 1)
    output_tcrdist <- output_tcrdist %>% select(CDR3b, clusters,size)
    colnames(output_tcrdist) <- c("CDR3", "clusters", "size")
    
    
    output_list_tcrdist[[paste("subset_",as.character(i),k,sep="")]] <- output_tcrdist
    
    df <- data.frame(version = k, subset = paste("subset",i,sep="_"), nb_cdr3 = n_distinct(output_list_tcrdist[[paste("subset_",as.character(i),k,sep="")]]$CDR3))
    clustered_cdr3_tcrdist <- rbind(clustered_cdr3_tcrdist, df)
  }
}



## Performance metrics

#1. Retention
subset_cdr3_df <- data.frame(subset_1=nrow(liste_subset[[1]]), subset_2=nrow(liste_subset[[2]]), 
                             subset_3=nrow(liste_subset[[3]]), subset_4=nrow(liste_subset[[4]]), 
                             subset_5=nrow(liste_subset[[5]]), subset_6=nrow(liste_subset[[6]])) %>%
  melt(., value.name="nb_cdr3", variable.name="subset")


retention_HD <- merge(clustered_cdr3_HD, subset_cdr3_df, by = "subset") %>% mutate(., retention = nb_cdr3.x/nb_cdr3.y) %>% mutate(methods ="Hamming Distance")
retention_HD_2 <- data_summary(retention_HD, varname="retention", groupnames=c("subset", "methods"))

retention_LD <- merge(clustered_cdr3_LD, subset_cdr3_df, by = "subset") %>% mutate(., retention = nb_cdr3.x/nb_cdr3.y) %>% mutate(methods ="Levenshtein Distance")
retention_LD_2 <- data_summary(retention_LD, varname="retention", groupnames=c("subset", "methods"))

retention_tcrm <- merge(clustered_cdr3_tcrm, subset_cdr3_df, by = "subset") %>% mutate(., retention = nb_cdr3.x/nb_cdr3.y) %>% mutate(methods ="TCRMatch")
retention_tcrm_2 <- data_summary(retention_tcrm, varname="retention", groupnames=c("subset", "methods"))

retention_ismart <- merge(clustered_cdr3_ismart, subset_cdr3_df, by = "subset") %>% mutate(., retention = nb_cdr3.x/nb_cdr3.y) %>% mutate(methods ="iSMART")
retention_ismart_2 <- data_summary(retention_ismart, varname="retention", groupnames=c("subset", "methods"))

retention_giana <- merge(clustered_cdr3_giana, subset_cdr3_df, by = "subset") %>% mutate(., retention = nb_cdr3.x/nb_cdr3.y) %>% mutate(methods ="GIANA")
retention_giana_2 <- data_summary(retention_giana, varname="retention", groupnames=c("subset", "methods"))

retention_clustcr <- merge(clustered_cdr3_clustcr, subset_cdr3_df, by = "subset") %>% mutate(., retention = nb_cdr3.x/nb_cdr3.y) %>% mutate(methods ="clusTCR")
retention_clustcr_2 <- data_summary(retention_clustcr, varname="retention", groupnames=c("subset", "methods"))

retention_gliph2 <- merge(clustered_cdr3_gliph2, subset_cdr3_df, by = "subset") %>% mutate(., retention = nb_cdr3.x/nb_cdr3.y) %>% mutate(methods ="GLIPH2")
retention_gliph2_2 <- data_summary(retention_gliph2, varname="retention", groupnames=c("subset", "methods"))

retention_deeptcr <- merge(clustered_cdr3_deeptcr, subset_cdr3_df, by = "subset") %>% mutate(., retention = nb_cdr3.x/nb_cdr3.y) %>% mutate(methods ="DeepTCR")
retention_deeptcr_2 <- data_summary(retention_deeptcr, varname="retention", groupnames=c("subset", "methods"))

retention_tcrdist <- merge(clustered_cdr3_tcrdist, subset_cdr3_df, by = "subset") %>% mutate(., retention = nb_cdr3.x/nb_cdr3.y) %>% mutate(methods ="TCRdist3")
retention_tcrdist_2 <- data_summary(retention_tcrdist, varname="retention", groupnames=c("subset", "methods"))




#2. Cluster Purity
labelled_cdr3b <- read.csv(paste(directory, "Data/10X/Donor1/labelled_cdr3b.csv", sep=""))
labelled_cdr3b <- labelled_cdr3b[-1]

labelled_cdr3b_list <- labelled_cdr3b %>% melt() %>% filter(value!=0)#2876 annotated CDR3b
colnames(labelled_cdr3b_list)[1] <- "CDR3"
colnames(labelled_cdr3b_list)[2] <- "Epitope"


purity_function <- function(output, th){
  output <- output %>% filter(size >1)
  output_spe <- merge(output, labelled_cdr3b_list, by = "CDR3", all.x=T) %>% unique()
  output_spe <- output_spe %>% group_by(clusters, Epitope) %>% dplyr::mutate(nb_epi_cluster = sum(!is.na(Epitope)))
  output_spe <- output_spe %>% group_by(clusters) %>% dplyr::mutate(nb_majoritaire = max(nb_epi_cluster))
  output_spe_sum <- output_spe %>% select(c(clusters, size, nb_majoritaire)) %>% unique()
  
  return(sum(output_spe_sum$nb_majoritaire)/sum(output_spe_sum$size))
}


purity_HD <- data.frame()
for(subset in names(output_list_HD)){
  score <- purity_function(output_list_HD[[subset]])
  df <- data.frame(subset = str_sub(subset,-9,-2), version = str_sub(subset,-1,-1), purity = score, methods = "Hamming Distance")
  purity_HD <- rbind(purity_HD , df)
}

purity_HD_2 <- data_summary(purity_HD, varname="purity", groupnames=c("subset", "methods"))


purity_LD <- data.frame()
for(subset in names(output_list_LD)){
  score <- purity_function(output_list_LD[[subset]])
  df <- data.frame(subset = str_sub(subset,-9,-2), version = str_sub(subset,-1,-1), purity = score, methods = "Levenshtein Distance")
  purity_LD <- rbind(purity_LD , df)
}

purity_LD_2 <- data_summary(purity_LD, varname="purity", groupnames=c("subset", "methods"))


purity_tcrm <- data.frame()
for(subset in names(output_list_tcrm)){
  score <- purity_function(output_list_tcrm[[subset]])
  df <- data.frame(subset = str_sub(subset,-9,-2), version = str_sub(subset,-1,-1), purity = score, methods = "TCRMatch")
  purity_tcrm <- rbind(purity_tcrm , df)
}

purity_tcrm_2 <- data_summary(purity_tcrm, varname="purity", groupnames=c("subset", "methods"))


purity_ismart <- data.frame()
for(subset in names(output_list_ismart)){
  score <- purity_function(output_list_ismart[[subset]])
  df <- data.frame(subset = str_sub(subset,-9,-2), version = str_sub(subset,-1,-1), purity = score, methods = "iSMART")
  purity_ismart <- rbind(purity_ismart , df)
}

purity_ismart_2 <- data_summary(purity_ismart, varname="purity", groupnames=c("subset", "methods"))


purity_giana <- data.frame()
for(subset in names(output_list_giana)){
  score <- purity_function(output_list_giana[[subset]])
  df <- data.frame(subset = str_sub(subset,-9,-2), version = str_sub(subset,-1,-1), purity = score, methods = "GIANA")
  purity_giana <- rbind(purity_giana , df)
}

purity_giana_2 <- data_summary(purity_giana, varname="purity", groupnames=c("subset", "methods"))


purity_clustcr <- data.frame()
for(subset in names(output_list_clustcr)){
  score <- purity_function(output_list_clustcr[[subset]])
  df <- data.frame(subset = str_sub(subset,-9,-2), version = str_sub(subset,-1,-1), purity = score, methods = "clusTCR")
  purity_clustcr <- rbind(purity_clustcr , df)
}

purity_clustcr_2 <- data_summary(purity_clustcr, varname="purity", groupnames=c("subset", "methods"))


purity_gliph2 <- data.frame()
for(subset in names(output_list_gliph2)){
  score <- purity_function(output_list_gliph2[[subset]])
  df <- data.frame(subset = str_sub(subset,-9,-2), version = str_sub(subset,-1,-1), purity = score, methods = "GLIPH2")
  purity_gliph2 <- rbind(purity_gliph2 , df)
}

purity_gliph2_2 <- data_summary(purity_gliph2, varname="purity", groupnames=c("subset", "methods"))


purity_deeptcr <- data.frame()
for(subset in names(output_list_deeptcr)){
  score <- purity_function(output_list_deeptcr[[subset]])
  df <- data.frame(subset = str_sub(subset,-9,-2), version = str_sub(subset,-1,-1), purity = score, methods = "DeepTCR")
  purity_deeptcr <- rbind(purity_deeptcr , df)
}

purity_deeptcr_2 <- data_summary(purity_deeptcr, varname="purity", groupnames=c("subset", "methods"))


purity_tcrdist <- data.frame()
for(subset in names(output_list_tcrdist)){
  score <- purity_function(output_list_tcrdist[[subset]])
  df <- data.frame(subset = str_sub(subset,-9,-2), version = str_sub(subset,-1,-1), purity = score, methods = "TCRdist3")
  purity_tcrdist <- rbind(purity_tcrdist , df)
}

purity_tcrdist_2 <- data_summary(purity_tcrdist, varname="purity", groupnames=c("subset", "methods"))

#test
# clustering_purity_bis <- function(list_output, matrix_epitope_sequence){
#   global_purity <- c()
#   for (i in 1:length(list_output)) {
#     print(i)
#     metadata_epitope_sequence <- merge(list_output[[i]], matrix_epitope_sequence, by = "cdr3") %>% select(., - cdr3) %>% sum_by_cluster(., "clusters", "cluster_size") %>% filter(., cluster_size != 1)
#     purity_epitope_sequence <- data.frame(cluster = metadata_epitope_sequence$clusters, size = metadata_epitope_sequence$cluster_size, epitope_maj = do.call(pmax, metadata_epitope_sequence[,3:ncol(metadata_epitope_sequence)])) %>% filter(., epitope_maj!=0)
#     global_purity <- append(global_purity, sum(purity_epitope_sequence$epitope_maj)/sum(purity_epitope_sequence$size))
#   } 
#   return(global_purity)
# }

#dans le second calcul de la pureté, les clusters qui ne possède que des CDR non nooté on été retiré. C'est pour cela que la pureté est supérieure par rapport au premier calcul

dataframe_retention_all <- rbind(retention_LD_2, retention_deeptcr_2, retention_HD_2, retention_gliph2_2, retention_tcrm_2, retention_clustcr_2, retention_ismart_2, retention_giana_2, retention_tcrdist_2)
dataframe_purity_all <- rbind(purity_LD_2, purity_deeptcr_2, purity_HD_2, purity_gliph2_2, purity_tcrm_2, purity_clustcr_2, purity_ismart_2, purity_giana_2, purity_tcrdist_2)

dataframe_retention_all$methods <- factor(dataframe_retention_all$methods, level = c("DeepTCR", "GLIPH2", "Hamming Distance", "Levenshtein Distance","TCRMatch", "clusTCR","GIANA", "iSMART", "TCRdist3"))
dataframe_purity_all$methods <- factor(dataframe_purity_all$methods, level = c("DeepTCR", "GLIPH2", "Hamming Distance", "Levenshtein Distance","TCRMatch", "clusTCR","GIANA", "iSMART", "TCRdist3"))


library(ggplot2)
plot_retention_all <- ggplot(dataframe_retention_all, aes(x= subset, y= retention, group = methods, color = methods)) + 
  geom_point(aes(color = methods)) +
  geom_line(aes(color = methods)) +
  scale_color_manual(values=c("chartreuse3", "pink", "darkorange", "#21908CFF", "#440154FF","#FFCC00", "brown1", "grey", "black")) +
  theme_bw() +
  labs(title= "Retention")+
  theme(legend.position = "right")+
  geom_errorbar(aes(ymin=retention-sd, ymax=retention+sd), width=.2,
                position=position_dodge(0.05))

plot_purity_all <- ggplot(dataframe_purity_all, aes(x= subset, y= purity, group = methods, color = methods)) + 
  geom_point(aes(color = methods)) +
  geom_line(aes(color = methods)) +
  scale_color_manual(values=c("chartreuse3", "pink", "darkorange", "#21908CFF", "#440154FF","#FFCC00", "brown1", "grey", "black")) +
  theme_bw() +
  labs(title= "Purity")+
  theme(legend.position = "right")+
  geom_errorbar(aes(ymin=purity-sd, ymax=purity+sd), width=.2,
                position=position_dodge(0.05))+
  ylim(0,1)


#3. Proportions of non annotated CDR3 per clusters
output_all <- list("Hamming Distance"=output_list_HD,
                   "Levenshtein Distance"=output_list_LD,
                   "TCRMatch"=output_list_tcrm,
                   "iSMART"=output_list_ismart,
                   "GIANA"=output_list_giana,
                   "ClusTCR"=output_list_clustcr,
                   "GLIPH2"=output_list_gliph2,
                   "TCRdist3"=output_list_tcrdist,
                   "DeepTCR"=output_list_deeptcr)

list_method <- c("Hamming Distance", "Levenshtein Distance", "TCRMatch", "iSMART","GIANA", "ClusTCR","GLIPH2", "DeepTCR", "TCRdist3")

output_all_df <- data.frame()
for(tool in names(output_all)){
  print(tool)
  for(i in 2:6){
    for(k in version_list){
      output <- output_all[[tool]][[paste("subset_",as.character(i),k,sep="")]]
      output_spe <- merge(output, labelled_cdr3b_list, by = "CDR3", all.x=T) %>% unique()
      output_spe <- output_spe %>% mutate(interval_size = ifelse(size<4, "[2-3]", 
                                                                 ifelse(size >3 & size<11, "[4-10]",
                                                                        ifelse(size>10 & size<51, "[11-50]",
                                                                               ifelse(size>50 & size <101, "[50-100]",
                                                                                      ifelse(size>100, ">100", NA))))))
      
      
      output_spe <- output_spe %>% group_by(interval_size) %>% dplyr::mutate(nb_non_annotated = sum(is.na(Epitope)))
      output_spe <- output_spe %>% mutate(percentage = nb_non_annotated/sum(is.na(output_spe$Epitope))*100) %>% 
        mutate(methods = tool, subset = paste0("subset_",i,sep=""), version =k)
      
      output_all_df <-rbind(output_all_df, output_spe)
    }
  }
}

test <- output_all_df %>% select(interval_size, percentage, methods, subset, version) %>% unique()

v <- "c"

output_all_df_v <- test %>% filter(version==v)
output_all_df_v$interval_size <- factor(output_all_df_v$interval_size, levels= rev(c("[2-3]", "[4-10]", "[11-50]", "[50-100]",">100")))
output_all_df_v$methods <- factor(output_all_df_v$methods, level = list_method)

ggplot(data = output_all_df_v, aes(x=methods, y = percentage, fill = interval_size))+#y=count
  geom_bar(stat="identity")+
  #scale_fill_grey(start=0.8, end=0.2) +
  scale_fill_brewer(palette="Dark2", direction=-1)+
  theme_bw()+
  ylab("Percentage of non annotated CDR3b (%)")+
  xlab("")+
  theme(legend.position ="bottom",
        #axis.text.x = element_blank())+
        axis.text.x = element_text(hjust=1, angle = 45, size =10),
        plot.title = element_text(hjust=0.5, face = "bold", size = 20),
        text=element_text(size=15),
        axis.text=element_text(size=13),
        axis.title.y = element_text(size = 12))+
  facet_wrap(~subset, ncol=5)+ggtitle(paste0("Version ", v, sep=""))



#4. Proportions of annotated CDR3 per clusters
output_all_df_annotated <- data.frame()
perc_annotated_cdr3_df <- data.frame()
for(tool in names(output_all)){
  print(tool)
  for(i in 2:6){
    for(k in version_list){
      output <- output_all[[tool]][[paste("subset_",as.character(i),k,sep="")]]
      output_spe <- merge(output, labelled_cdr3b_list, by = "CDR3") %>% unique()
      output_spe <- output_spe %>% mutate(interval_size = ifelse(size<4, "[2-3]", 
                                                                 ifelse(size >3 & size<11, "[4-10]",
                                                                        ifelse(size>10 & size<51, "[11-50]",
                                                                               ifelse(size>50 & size <101, "[50-100]",
                                                                                      ifelse(size>100, ">100", NA))))))
      
      
      output_spe <- output_spe %>% group_by(interval_size) %>% dplyr::mutate(nb_annotated = sum(!is.na(Epitope)))
      output_spe <- output_spe %>% mutate(percentage = nb_annotated/sum(!is.na(output_spe$Epitope))*100) %>% 
        mutate(methods = tool, subset = paste0("subset_",i,sep=""), version =k)
      
      output_all_df_annotated <-rbind(output_all_df_annotated, output_spe)
      
      nb_annotated <- n_distinct(output_spe$CDR3)
      df_2 <- data.frame(methods = tool, version = k, subset = paste0("subset_",i,sep=""), nb_annotated_cdr3 =nb_annotated, perc_annotated_cdr3 = nb_annotated/2876)
      perc_annotated_cdr3_df <- rbind(perc_annotated_cdr3_df , df_2)
    }
  }
}

test <- output_all_df_annotated %>% select(interval_size, percentage, methods, subset, version) %>% unique()

v <- "c"

output_all_df_v <- test %>% filter(version==v)
output_all_df_v$interval_size <- factor(output_all_df_v$interval_size, levels= rev(c("[2-3]", "[4-10]", "[11-50]", "[50-100]",">100")))
output_all_df_v$methods <- factor(output_all_df_v$methods, level = list_method)

ggplot(data = output_all_df_v, aes(x=methods, y = percentage, fill = interval_size))+#y=count
  geom_bar(stat="identity")+
  #scale_fill_grey(start=0.8, end=0.2) +
  scale_fill_brewer(palette="Dark2", direction=-1)+
  theme_bw()+
  ylab("Percentage of non annotated CDR3b (%)")+
  xlab("")+
  theme(legend.position ="bottom",
        #axis.text.x = element_blank())+
        axis.text.x = element_text(hjust=1, angle = 45, size =10),
        plot.title = element_text(hjust=0.5, face = "bold", size = 20),
        text=element_text(size=15),
        axis.text=element_text(size=13),
        axis.title.y = element_text(size = 12))+
  facet_wrap(~subset, ncol=5)+ggtitle(paste0("Version ", v, sep=""))


perc_annotated_cdr3_df_2 <- data_summary(perc_annotated_cdr3_df , varname="perc_annotated_cdr3", groupnames=c("subset", "methods"))
perc_annotated_cdr3_df_2$methods <- factor(perc_annotated_cdr3_df_2$methods, level = c("DeepTCR", "GLIPH2", "Hamming Distance", "Levenshtein Distance","TCRMatch", "ClusTCR","GIANA", "iSMART", "TCRdist3"))

ggplot(perc_annotated_cdr3_df_2, aes(x= subset, y= perc_annotated_cdr3, group = methods, color = methods)) + 
  geom_point(aes(color = methods)) +
  geom_line(aes(color = methods)) +
  scale_color_manual(values=c("chartreuse3", "pink", "darkorange", "#21908CFF", "#440154FF","#FFCC00", "brown1", "grey", "black")) +
  theme_bw() +
  labs(title= "Percentage_annotated CDR3")+
  theme(legend.position = "right")+
  geom_errorbar(aes(ymin=perc_annotated_cdr3-sd, ymax=perc_annotated_cdr3+sd), width=.2,
                position=position_dodge(0.05))

#tableau
perc_annotated_cdr3_df
retention_all <- rbind(retention_LD, retention_deeptcr, retention_HD, retention_gliph2, retention_tcrm, retention_clustcr, retention_ismart, retention_giana, retention_tcrdist)
purity_all <- rbind(purity_LD, purity_deeptcr, purity_HD, purity_gliph2, purity_tcrm, purity_clustcr, purity_ismart, purity_giana, purity_tcrdist)

big_table <- merge(retention_all, purity_all, by=c("methods", "subset","version")) %>% merge(., perc_annotated_cdr3_df, by=c("methods", "subset","version"))
big_table <- big_table %>% select(-c(nb_cdr3.x, nb_cdr3.y))

write.csv(big_table, paste(directory, "Plot/metric_table.csv", sep=""), quote = F)

#5. Retention of annotated CDR3
test <- merge(output_list_gliph2[["subset_1a"]], output_list_gliph2[["subset_2a"]], by="CDR3")
n_distinct(test$CDR3)/n_distinct(output_list_gliph2[["subset_1a"]]$CDR3)

retention_annotated_all <- data.frame()
for(tool in names(output_all)){
  print(tool)
  for(i in 2:6){
    print(i)
    for(k in version_list){
      print(k)
      output_annotated <- output_all[[tool]][[paste("subset_1",k,sep="")]]
      output <- output_all[[tool]][[paste("subset_",as.character(i),k,sep="")]]
      output_merge <- merge(output_annotated, output, by="CDR3")
      score <- n_distinct(output_merge$CDR3)/n_distinct(output_annotated$CDR3)*100
      df <- data.frame(methods = tool, version = k, subset = paste("subset_",as.character(i),sep=""), retention_annotated = score)
      retention_annotated_all <- rbind(retention_annotated_all, df)
    }
  }
}
retention_annotated_all_2 <- data_summary(retention_annotated_all, varname="retention_annotated", groupnames=c("subset", "methods"))

retention_annotated_all_2$methods <- factor(retention_annotated_all_2$methods, level = c("DeepTCR", "GLIPH2", "Hamming Distance", "Levenshtein Distance","TCRMatch", "ClusTCR","GIANA", "iSMART", "TCRdist3"))

ggplot(retention_annotated_all_2, aes(x= subset, y= retention_annotated, group = methods, color = methods)) + 
  geom_point(aes(color = methods)) +
  geom_line(aes(color = methods)) +
  scale_color_manual(values=c("chartreuse3", "pink", "darkorange", "#21908CFF", "#440154FF","#FFCC00", "brown1", "grey", "black")) +
  theme_bw() +
  labs(title= "Percentage of annotated clustered CDR3s relative to subset 1")+
  theme(legend.position = "right")+ylab("Percentage")+
  geom_errorbar(aes(ymin=retention_annotated-sd, ymax=retention_annotated+sd), width=.2,
                position=position_dodge(0))


#6. Percentage of non annotated CDR3 in clusters
perc_non_annotated_cluster <- data.frame()
for(tool in names(output_all)){
  print(tool)
  for(i in 2:6){
    for(k in version_list){
      output <- output_all[[tool]][[paste("subset_",as.character(i),k,sep="")]]
      output_spe <- merge(output, labelled_cdr3b_list, by = "CDR3", all.x=T) %>% unique()
      output_spe <- output_spe %>% group_by(clusters) %>% dplyr::mutate(nb_non_annotated_per_cluster = sum(is.na(Epitope)))
      output_spe <- output_spe %>% dplyr::mutate(perc_non_annotated_per_cluster = nb_non_annotated_per_cluster/size*100)
      
      output_spe <- output_spe %>% mutate(interval_per = ifelse(perc_non_annotated_per_cluster<10, "<10%", 
                                                                 ifelse(perc_non_annotated_per_cluster >10 & perc_non_annotated_per_cluster<31, "10-30%",
                                                                        ifelse(perc_non_annotated_per_cluster>30 & perc_non_annotated_per_cluster<61, "30-60%",
                                                                               ifelse(perc_non_annotated_per_cluster>60 & perc_non_annotated_per_cluster <91, "60-90%",
                                                                                      ifelse(perc_non_annotated_per_cluster>90 & perc_non_annotated_per_cluster<100, "90-99%", 
                                                                                             ifelse(perc_non_annotated_per_cluster==100, "100%", NA)))))))
      
      
      output_spe <- output_spe %>% group_by(interval_per) %>% dplyr::mutate(nb_cluster= n_distinct(clusters))
      output_spe <- output_spe %>% dplyr::mutate(nb_total = n_distinct(output_spe$clusters)) %>% 
        mutate(methods = tool, subset = paste0("subset_",i,sep=""), version =k)
      output_spe <- output_spe %>% dplyr::mutate(percentage = nb_cluster/nb_total*100)
      
      perc_non_annotated_cluster  <-rbind(perc_non_annotated_cluster, output_spe)
    }
  }
}


test <- perc_non_annotated_cluster%>% select(interval_per, percentage, methods, subset, version) %>% unique()

v <- "c"

output_all_df_v <- test %>% filter(version==v)
output_all_df_v$interval_per <- factor(output_all_df_v$interval_per, levels= rev(c("<10%", "10-30%", "30-60%", "60-90%","90-99%", "100%")))
output_all_df_v$methods <- factor(output_all_df_v$methods, level = list_method)

ggplot(data = output_all_df_v, aes(x=methods, y = percentage, fill = interval_per))+#y=count
  geom_bar(stat="identity")+
  #scale_fill_grey(start=0.8, end=0.2) +
  scale_fill_brewer(palette="Dark2", direction=-1)+
  theme_bw()+
  ylab("Percentage of cluster with non annotated CDR3b (%)")+
  xlab("")+
  theme(legend.position ="bottom",
        #axis.text.x = element_blank())+
        axis.text.x = element_text(hjust=1, angle = 45, size =10),
        plot.title = element_text(hjust=0.5, face = "bold", size = 20),
        text=element_text(size=15),
        axis.text=element_text(size=13),
        axis.title.y = element_text(size = 12))+
  facet_wrap(~subset, ncol=5)+ggtitle(paste0("Version ", v, sep=""))


                                    