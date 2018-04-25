library("sleuth")
library("biomaRt")
setwd("/Users/elizabethboylesobolik/Desktop/sleuth_objects/")
p <- "/Users/elizabethboylesobolik/Desktop"
r <- "sleuth_objects"
sleuth_objs <- (dir(file.path(p,r)))
print(sleuth_objs)
for (s in sleuth_objs){
  so <- sleuth_load(s)
  models <- models(so)
  full_model <- models[1]
  fm <- full_model[[1]]
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