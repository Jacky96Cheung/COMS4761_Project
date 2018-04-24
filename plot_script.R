library("sleuth")
library("biomaRt")
p <- "/Users/elizabethboylesobolik/Desktop/COMS4761_Project"
r <- "sleuth_objects"
sleuth_objs <- (dir(file.path(p,r)))
sos <- sleuth_objs[6:10]
for (s in sos){
  so <- sleuth_load(s)
  print(models(so))
  test <- "conditionIP_Sucrose"
  volcano <- plot_volcano(so, test, test_type = "wt", which_model = "full",
               sig_level = 0.1, point_alpha = 0.2, sig_color = "red",
               highlight = NULL)
  file_name <- paste("volcano", s, ".pdf", sep = "")
  ggsave(file_name)
}