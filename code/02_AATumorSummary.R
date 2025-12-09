#-------------------------------------------------------------------------------
# 02_AATumorSummary.R
#-------------------------------------------------------------------------------
library(dplyr)
library(data.table)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(ggsignif)
library(scales)
theme_cb <- function(color="black"){theme_classic() + theme(axis.text = element_text(color=color), axis.ticks = element_line(color=color))}

SampleInfo <- read.csv("output/Tables/00_SampleListSummary_Curated.csv", row.names = 1)

#----- AA output summary -----#
AA_result <- mclapply(SampleInfo$SampleID_T,
                      function(x) fread(paste0("../data/RE-01_AAresults/", x, "/", x, "_classification/", x,  "_result_table.tsv")),
                      mc.cores = 20) 
AA_result <- do.call(rbind, AA_result)
write.csv(AA_result, "output/Tables/02_AA_Summary_Tumor.csv")

AA_result <- read.csv("output/Tables/02_AA_Summary_Tumor.csv", row.names = 1)
AA_result_T <- AA_result[!is.na(AA_result$Classification),]
length(unique(AA_result_T$Sample.name)) # 581 samples, 2088 amplicons (prefilter)

#----- Compare amplicons between normal and tumor -----#
AA_result_N <- read.csv("output/Tables/01_AA_Summary_Normal.csv", row.names = 1)
AA_result_N <- AA_result_N[!is.na(AA_result_N$Classification),]

#number of genes in each amplicons
AA_result_T_geneNo <- lapply(c(1:nrow(AA_result_T)), function(i){
  input_string <- AA_result_T[i,]$All.genes
  clean_string <- gsub("\\[", "", input_string)
  clean_string <- gsub("\\]", "", clean_string)
  clean_string <- gsub("'", "", clean_string)
  clean_string <- gsub(" ", "", clean_string)
  result <- unlist(strsplit(clean_string, ","))
  result <- length(result)
  return(result)
}) %>% unlist

AA_result_N_geneNo <- lapply(c(1:nrow(AA_result_N)), function(i){
  input_string <- AA_result_N[i,]$All.genes
  clean_string <- gsub("\\[", "", input_string)
  clean_string <- gsub("\\]", "", clean_string)
  clean_string <- gsub("'", "", clean_string)
  clean_string <- gsub(" ", "", clean_string)
  result <- unlist(strsplit(clean_string, ","))
  result <- length(result)
  return(result)
}) %>% unlist

#compare normal vs tumor
df1 <- data.frame(amplicon = AA_result_T$Feature.ID,
                  class = AA_result_T$Classification,
                  CN = AA_result_T$Feature.maximum.copy.number,
                  length = AA_result_T$Captured.interval.length,
                  NumGenes = AA_result_T_geneNo,
                  sample = "Tumor")
df2 <- data.frame(amplicon = AA_result_N$Feature.ID,
                  class = AA_result_N$Classification,
                  CN = AA_result_N$Feature.maximum.copy.number,
                  length = AA_result_N$Captured.interval.length,
                  NumGenes = AA_result_N_geneNo,
                  sample = "Normal")
df3 <- rbind(df1, df2)

p1 <- ggplot(df3, aes(x = sample, y = CN)) + geom_jitter(height = 0, width = 0.25, shape = 3, alpha = 0.5, aes(color = sample)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() + 
  scale_color_manual(values = c("Tumor" = "#D51F26", "Normal" = "#272E6A")) +
  labs(x = "Sample Type", y = "WGS-derived CN (Maximum segment)") + 
  scale_y_continuous(trans = "log2", breaks = c(0,2,4,8,16,32,64,128,256)) + 
  geom_signif(comparisons = list(c("Tumor", "Normal")), test = "wilcox.test")
p2 <- ggplot(df3, aes(x = sample, y = length)) + geom_jitter(height = 0, width = 0.25, shape = 3, alpha = 0.5, aes(color = sample)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() + 
  scale_y_log10(breaks=10^(4:8), labels=trans_format("log10",math_format(10^.x))) +
  scale_color_manual(values = c("Tumor" = "#D51F26", "Normal" = "#272E6A")) +
  labs(x = "Sample Type", y = "Amplicon length (bp)") +
  geom_signif(comparisons = list(c("Tumor", "Normal")), test = "wilcox.test")
p3 <- ggplot(df3, aes(x = sample, y = NumGenes)) + geom_jitter(height = 0, width = 0.25, shape = 3, alpha = 0.5, aes(color = sample)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() + 
  scale_color_manual(values = c("Tumor" = "#D51F26", "Normal" = "#272E6A")) +
  labs(x = "Sample Type", y = "# of genes on amplicons") +
  geom_signif(comparisons = list(c("Tumor", "Normal")), test = "wilcox.test")

pdf("output/Plots/02_Boxplot_CompareNormalTumor.pdf", width = 3.5, height = 5)
p1
p2
p3
dev.off()

p4 <- ggplot(df3, aes(x = NumGenes, fill = sample)) + geom_bar() + theme_cb() + facet_grid(~sample) + 
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(breaks = c(seq(0,350,50))) +
  scale_fill_manual(values = c("Tumor" = "#D51F26", "Normal" = "#272E6A")) +
  labs(x = "# of genes on amplicons", y = "# of amplicons")
p5 <- ggplot(df3, aes(x = NumGenes, fill = sample)) + geom_bar() + theme_cb() + facet_grid(~sample) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(limits = c(-1,50), breaks = c(seq(0,50,10))) +
  scale_fill_manual(values = c("Tumor" = "#D51F26", "Normal" = "#272E6A")) +
  labs(x = "# of genes on amplicons", y = "# of amplicons")

pdf("output/Plots/02_Barplot_NumGenes_CompareNormalTumor.pdf", width = 6, height = 4)
p4
p5
dev.off()

## Overlap amplicons between normal and cancer (tech artifacts) ##
gr_N <- readRDS("rds/01_gr_amplicon_normal.rds")

AA_result_T <- AA_result_T[which(AA_result_T$Location != "[]"),]
gr_T <- mclapply(c(1:nrow(AA_result_T)), function(i){
  print(i)
  bedPath <- gsub("/home/output/", "", AA_result_T$Feature.BED.file[i])
  bedPath <- paste0("../data/RE-01_AAresults/", AA_result_T$Sample.name[i], "/", bedPath)
  gr <- rtracklayer::import.bed(bedPath)
  gr$FeatureID <- AA_result_T$Feature.ID[i]
  gr$Classification <- AA_result_T$Classification[i]
  gr$CN <- AA_result_T$Feature.maximum.copy.number[i]
  gr$length <- AA_result_T$Captured.interval.length[i]
  return(gr)
}, mc.cores = 20)
gr_T <- GRangesList(gr_T)

fo1 <- findOverlaps(gr_T, gr_N)
idx1 <- unique(queryHits(fo1)) #397 amplicons are overlapped
gr_T_filtered <- gr_T[-idx1] #1690 amplicons

idx2 <- lapply(gr_T_filtered, function(x) unique(x$FeatureID))
idx2 <- unlist(idx2)

AA_result_T_filt <- AA_result_T[which(AA_result_T$Feature.ID %in% idx2),]
write.csv(AA_result_T_filt, "output/Tables/02_AA_Summary_Tumor_filtered.csv")
saveRDS(gr_T_filtered, "rds/02_gr_T_filtered.rds")

#----- After filtering: Compare amplicons between normal and tumor -----#
#number of genes in each amplicons
AA_result_T_geneNo_Filt <- lapply(c(1:nrow(AA_result_T_filt)), function(i){
  input_string <- AA_result_T_filt[i,]$All.genes
  clean_string <- gsub("\\[", "", input_string)
  clean_string <- gsub("\\]", "", clean_string)
  clean_string <- gsub("'", "", clean_string)
  clean_string <- gsub(" ", "", clean_string)
  result <- unlist(strsplit(clean_string, ","))
  result <- length(result)
  return(result)
}) %>% unlist

df1 <- data.frame(amplicon = AA_result_T_filt$Feature.ID,
                  class = AA_result_T_filt$Classification,
                  CN = AA_result_T_filt$Feature.maximum.copy.number,
                  length = AA_result_T_filt$Captured.interval.length,
                  NumGenes = AA_result_T_geneNo_Filt,
                  sample = "Tumor")
df2 <- data.frame(amplicon = AA_result_N$Feature.ID,
                  class = AA_result_N$Classification,
                  CN = AA_result_N$Feature.maximum.copy.number,
                  length = AA_result_N$Captured.interval.length,
                  NumGenes = AA_result_N_geneNo,
                  sample = "Normal")
df3 <- rbind(df1, df2)

p6 <- ggplot(df3, aes(x = sample, y = CN)) + geom_jitter(height = 0, width = 0.25, shape = 3, alpha = 0.5, aes(color = sample)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() + 
  scale_color_manual(values = c("Tumor" = "#D51F26", "Normal" = "#272E6A")) +
  labs(x = "Sample Type", y = "WGS-derived CN (Maximum segment)") + 
  scale_y_continuous(trans = "log2", breaks = c(0,2,4,8,16,32,64,128,256)) + 
  geom_signif(comparisons = list(c("Tumor", "Normal")), test = "wilcox.test")
p7 <- ggplot(df3, aes(x = sample, y = length)) + geom_jitter(height = 0, width = 0.25, shape = 3, alpha = 0.5, aes(color = sample)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() + 
  scale_y_log10(breaks=10^(4:8), labels=trans_format("log10",math_format(10^.x))) +
  scale_color_manual(values = c("Tumor" = "#D51F26", "Normal" = "#272E6A")) +
  labs(x = "Sample Type", y = "Amplicon length (bp)") +
  geom_signif(comparisons = list(c("Tumor", "Normal")), test = "wilcox.test")
p8 <- ggplot(df3, aes(x = sample, y = NumGenes)) + geom_jitter(height = 0, width = 0.25, shape = 3, alpha = 0.5, aes(color = sample)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() + 
  scale_color_manual(values = c("Tumor" = "#D51F26", "Normal" = "#272E6A")) +
  labs(x = "Sample Type", y = "# of genes on amplicons") +
  geom_signif(comparisons = list(c("Tumor", "Normal")), test = "wilcox.test")

pdf("output/Plots/02_Boxplot_CompareNormalTumor_AfterFilt.pdf", width = 3.5, height = 5)
p6
p7
p8
dev.off()

p9 <- ggplot(df3, aes(x = NumGenes, fill = sample)) + geom_bar() + theme_cb() + facet_grid(~sample) + 
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(breaks = c(seq(0,350,50))) +
  scale_fill_manual(values = c("Tumor" = "#D51F26", "Normal" = "#272E6A")) +
  labs(x = "# of genes on amplicons", y = "# of amplicons")
p10 <- ggplot(df3, aes(x = NumGenes, fill = sample)) + geom_bar() + theme_cb() + facet_grid(~sample) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(limits = c(-1,50), breaks = c(seq(0,50,10))) +
  scale_fill_manual(values = c("Tumor" = "#D51F26", "Normal" = "#272E6A")) +
  labs(x = "# of genes on amplicons", y = "# of amplicons")

pdf("output/Plots/02_Barplot_NumGenes_CompareNormalTumor_AfterFilt.pdf", width = 6, height = 4)
p9
p10
dev.off()
