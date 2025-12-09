#get argument
args <- commandArgs(trailingOnly = TRUE)

#check
if (length(args) == 0) {
  stop("Please enter a valid argument!")
}

#as numeric
input_value <- as.numeric(args[1])

#input matrix
norm_counts <- read.csv("RNMF_InputMatrix.csv", row.names = 1, header = T)

library(NMF)
nmf_result <- nmf(norm_counts, rank = input_value, method = "brunet", nrun = 100)

saveRDS(nmf_result, paste0("nmf_result_k", input_value, ".rds"))
