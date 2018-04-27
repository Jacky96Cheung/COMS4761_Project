# A script that loads all the sleuth objects and outputs a volcano plot
# of the Wald beta-parameters and a csv file of the gene expression results

library("sleuth")
library("biomaRt")

# First, gather all the filepaths/names of the sleuth objects to open
setwd("/Users/elizabethboylesobolik/Desktop/sleuth_objects/")
p <- "/Users/elizabethboylesobolik/Desktop"
r <- "sleuth_objects"
sleuth_objs <- (dir(file.path(p,r)))

# Now open, plot each one!
for (s in sleuth_objs){
  # load the sleuth object
  so <- sleuth_load(s)
  
  # get a summary of the models computed
  models <- models(so)
  
  # get the pointer to the full model
  full_model <- models[1]
  # then get first information in the full model (beta values)
  fm <- full_model[[1]]
  
  # then get the value of the first beta value
  k <- fm[[1]]
  list_info <- k[[1]]
  condition_name <- list_info$b1
  test <- names(condition_name)
  print(test)
  volcano <- plot_volcano(so, test, test_type = "wt", which_model = "full",
               sig_level = 0.1, point_alpha = 0.2, sig_color = "red",
               highlight = NULL)
  s_results <- sleuth_results(so, test, test_type = "wt", which_model = "full",
                 rename_cols = TRUE, show_all = TRUE)
  file_path <- "/Users/elizabethboylesobolik/Desktop/sr/"
  volc_file_name <- paste(file_path,"volcano", s, ".pdf", sep = "")
  table_file_name <- paste(file_path, s,".csv", sep = "")
  ggsave(volc_file_name)
  write.csv(s_results, table_file_name)
}