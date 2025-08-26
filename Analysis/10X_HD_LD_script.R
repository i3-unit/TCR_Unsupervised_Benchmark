# Processing of clustering outputs
library(dplyr)
library(readr)
library(stringdist)
library(reshape2)

directory <- "E:/Labo_i3/Paper_Benchmark_v2/Rcode_Github/TCR_Unsupervised_Benchmark_modif/"


liste_subset <- list()
for (i in 1:6) {
  print(i)
  myfile <- file.path(paste0(directory,"Data/10X/Donor1/subsets"), paste0("subset", "_", i, ".txt"))
  liste_subset[[i]] <- read_delim(file=myfile, delim = ",") %>% data.frame()
  i <- i+1
}

# 1. HD output

file <- paste0(directory, "Data/10X_Output_methods/output_HD")
version_list <- c("a","b","c")

for(version in version_list){
  print(version)
    for(k in 1:length(liste_subset)){
    print(k)
    HD_matrix <- as.matrix(stringdistmatrix(liste_subset[[k]]$cdr3, method = "hamming",useNames = "strings", useBytes=TRUE ))
    output_HD <- melt(HD_matrix)
    colnames(output_HD) <- c("input_sequence", "match_sequence", "score")
    output_HD <- filter(output_HD, score == 1)
    graph_HD <- graph_from_data_frame(output_HD, directed = F)
    component <- components(graph_HD)
    output_HD <- data.frame(clusters = component$membership) %>% rownames_as_first_column(., "CDR3") %>% 
      group_by(clusters) %>% dplyr:: mutate(size = n())
    myfile <- file.path(file, paste0("subset", "_", k,version, ".csv"))
    write.table(output_HD, file = myfile, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, append = FALSE)
  }
}


# 2. LD output
file <- paste0(directory, "Data/10X_Output_methods/output_LD")

for(version in version_list){
  print(version)
  for(k in 1:length(liste_subset)){
    print(k)
    LD_matrix <- as.matrix(stringdistmatrix(liste_subset[[k]]$cdr3, method = "lv",useNames = "strings", useBytes=TRUE ))
    output_LD <- melt(LD_matrix)
    colnames(output_LD) <- c("input_sequence", "match_sequence", "score")
    output_LD <- filter(output_LD, score == 1)
    graph_LD <- graph_from_data_frame(output_LD, directed = F)
    component <- components(graph_LD)
    output_LD <- data.frame(clusters = component$membership) %>% rownames_as_first_column(., "CDR3") %>% 
      group_by(clusters) %>% dplyr:: mutate(size = n())
    myfile <- file.path(file, paste0("subset", "_", k,version, ".csv"))
    write.table(output_LD, file = myfile, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, append = FALSE)
  }
}


