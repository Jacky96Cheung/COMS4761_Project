# A script that loads all the sleuth objects and outputs a volcano plot
# of the Wald beta-parameters and a csv file of the gene expression results

library("sleuth")
library("biomaRt")

# File path specifications:
# where are the sleuth objects?
p <- "/Users/elizabethboylesobolik/Desktop"
r <- "sleuth_objects"

# where do you want to save the results?
file_path <- "/Users/elizabethboylesobolik/Desktop/COMS4761_Project/sr/"

# First, gather all the filepaths/names of the sleuth objects to open
setwd(p)
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

  # make a volcano plot from the beta values of a Wald Test
  volcano <- plot_volcano(so, test, test_type = "wt", which_model = "full",
               sig_level = 0.1, point_alpha = 0.2, sig_color = "red",
               highlight = NULL)
  s_results <- sleuth_results(so, test, test_type = "wt", which_model = "full",
                 rename_cols = TRUE, show_all = TRUE)
  # save it!
  volc_file_name <- paste(file_path,"volcano", s, ".pdf", sep = "")
  table_file_name <- paste(file_path, s,".csv", sep = "")
  
  # as an image 
  ggsave(volc_file_name)
  
  # and a table of results to match the data points on the plot to genes
  write.csv(s_results, table_file_name)
}