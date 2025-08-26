#function to create 6 subsets of CDR3
# - liste_cdr3_all : list of all CDR3
# - liste_cdr3_filtered : list of annotated CDR3 
create_subset <- function(liste_cdr3_all, liste_cdr3_filtered){
  subset <- liste_cdr3_filtered
  data <- anti_join(liste_cdr3_all, liste_cdr3_filtered)
  interval <- (nrow(liste_cdr3_all) - nrow(liste_cdr3_filtered))/5
  liste_subset <- list()
  liste_subset[[1]] <- subset
  for (i in 2:6){
    sample <- sample_n(data, interval)
    data <- anti_join(data, sample)
    liste_subset[[i]] <- rbind(subset, sample)
    subset <- liste_subset[[i]]
    print(nrow(subset))
  }
  return(liste_subset)
}