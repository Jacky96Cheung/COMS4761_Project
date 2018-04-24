get_sleuthobj <- function(p,r, cond1, cond2, t2g){
  sample_ids <- dir(file.path(p,r))
  print(sample_ids)
  cond1_id <- sample_ids[grepl(cond1, sample_ids)]
  print(cond1_id)
  cond2_id <- sample_ids[grepl(cond2, sample_ids)]
  sample_id <- c(cond1_id, cond2_id)
  conds1 <- rep(cond1, length(cond1_id))
  conds2 <- rep(cond2, length(cond2_id))
  conds <- c(conds1, conds2)
  kal_dirs <- file.path(p,r, sample_id)
  print(kal_dirs)
  s2c <- data_frame(sample = sample_id, condition = conds, path = kal_dirs)
  print(s2c)
  so <- sleuth_prep(s2c, target_mapping = t2g, extra_bootstrap_summary = TRUE)
  so <- sleuth_fit(so, ~condition, 'full')
  so <- sleuth_fit(so, ~1, 'reduced')
  beta <- paste('condition', cond1, sep = "")
  so <- sleuth_wt(so, beta, which_model = "full")
  f <- paste("/Users/elizabethboylesobolik/Desktop/COMS4761_Project/sleuth_objects/",cond1,cond2, sep = "")
  print(f)
  sleuth_save(so, f)
}

library("sleuth")
library("biomaRt")

print(listMarts())
mart <- useMart("ensembl",dataset="mmusculus_gene_ensembl")
t2g <- getBM(attributes = c("ensembl_transcript_id", 
                            "ensembl_gene_id", "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, 
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
p <- "/Users/elizabethboylesobolik/Desktop/COMS4761_Project"
r <- "results"
cr <- "CR"
ip <- "IP"
stims <- c("Intralipid", "LiCl", "RACEQTI", "GW7647", "IV")
for (s in stims ){
  cond1 <- paste(ip,"Sucrose", sep = "_")
  cond2 <- paste(ip, s, sep = "_")
  get_sleuthobj(p,r,cond1,cond2,t2g)
}


