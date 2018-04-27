# a function that takes a path (p), the directory where the results are stored (r)
# two conditions to compare (cond1, cond2) and an Ensembl gene id conversion table (t2g)

get_sleuthobj <- function(p,r, cond1, cond2, t2g){
  # collect all the conditions/replicate names (references to samples)
  sample_ids <- dir(file.path(p,r))
  
  # because of the way the files are named, we have the conditions in the filenames!
  # (collect the relevant data by searching the conditions as substrings) 
  cond1_id <- sample_ids[grepl(cond1, sample_ids)]
  cond2_id <- sample_ids[grepl(cond2, sample_ids)]
  sample_id <- c(cond1_id, cond2_id)
  
  # make a table that matches the samples to their conditions
  conds1 <- rep(cond1, length(cond1_id))
  conds2 <- rep(cond2, length(cond2_id))
  conds <- c(conds1, conds2)
  
  # sleuth needs the filepaths to the relevant sample data (kallisto directories)
  kal_dirs <- file.path(p,r, sample_id)
  
  # make the data frame that will be passed into the sleuth object constructer
  s2c <- data_frame(sample = sample_id, condition = conds, path = kal_dirs)

  # construct a sleuth object
  # the aggregation column is so that we are assessing differential expression at the
  # gene level and not the transcript level -- this jives with our later linear regression 
  # against the Allen Brain Altas ISH data, which is at the gene level ... 
  so <- sleuth_prep(s2c, target_mapping = t2g, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE)
  
  # make a general linear model (full means the gene identities are known)
  so <- sleuth_fit(so, ~condition, 'full')
  
  # make a general linear model (reduced is irrespective of gene identity)
  so <- sleuth_fit(so, ~1, 'reduced')
  
  # if not otherwise specified, sleuth selects the dependent condition based on alphabetical order
  # which works for us because CR comes before IP -- this is a way of identifying the dependent cond
  # for IP/IP comparisons
  
  if(cond1 < cond2){
    cond <- cond2
  }
  else{
    cond <- cond1
  }
  
  # make the string needed to do the wald test
  beta <- paste('condition', cond, sep = "")
  so <- sleuth_wt(so, beta, which_model = "full")
  
  # designate a file path, name for the sleuth object
  f <- paste("/Users/elizabethboylesobolik/Desktop/sleuth_objects/",cond1,cond2, sep = "")
  
  # save it!
  sleuth_save(so, f)
  print("Saved sleuth object!")
}

# load relevant libraries
library("sleuth")
library("biomaRt")

print(listMarts())

# Get the mouse Ensembl mart and the relevant labels into a transcript to genome dataframe
mart <- useMart("ensembl",dataset="mmusculus_gene_ensembl")
t2g <- getBM(attributes = c("ensembl_transcript_id", 
                            "ensembl_gene_id", "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, 
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

p <- "/Users/elizabethboylesobolik/Desktop/COMS4761_Project"
r <- "results"

# read in text file of the comparisons we want
s_meta <- read.table("sleuth_metadata.txt", header = TRUE)
for (row in 1:nrow(s_meta) ){
  cond1 <- as.character(s_meta[row, "condition1"])
  cond2 <- as.character(s_meta[row, "condition2"])
  get_sleuthobj(p,r,cond1,cond2,t2g)
}


