#-------------------------------------------------------------------------------
# 01_AANormalSummary.R
#-------------------------------------------------------------------------------
library(dplyr)
library(data.table)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

SampleInfo <- read.csv("output/Tables/00_SampleListSummary_Curated.csv", row.names = 1)

#----- AA summary for normal -----#
AA_result <- mclapply(SampleInfo$SampleID_N,
                      function(x) fread(paste0("../data/RE-01_AAresults_N/", x, "/", x, "_classification/", x,  "_result_table.tsv")),
                      mc.cores = 20) 
AA_result <- do.call(rbind, AA_result)
write.csv(AA_result, "output/Tables/01_AA_Summary_Normal.csv")

AA_result2 <- AA_result[!is.na(AA_result$Classification),]
length(unique(AA_result2$`Sample name`)) # 134 normal samples, 141 amplicons

#----- GrangesList -----#
gr_amplicons <- lapply(c(1:nrow(AA_result2)), function(i){
  bedPath <- gsub("/home/output/", "", AA_result2$`Feature BED file`[i])
  bedPath <- paste0("../data/RE-01_AAresults_N/", AA_result2$`Sample name`[i], "/", bedPath)
  gr <- rtracklayer::import.bed(bedPath)
  gr$CaseID <- SampleInfo$CaseID[which(SampleInfo$SampleID_N == AA_result2$`Sample name`[i])]
  gr$FeatureID <- AA_result2$`Feature ID`[i]
  gr$Classification <- AA_result2$Classification[i]
  return(gr)
})
gr_amplicons <- unlist(GRangesList(gr_amplicons))

#----- Overlap intervals -----#
overlapMtx <- matrix(data = 0, nrow = length(gr_amplicons), ncol = length(gr_amplicons))
for(i in seq_along(gr_amplicons)){
  idy <- findOverlaps(gr_amplicons[i], gr_amplicons) %>% subjectHits()
  overlapMtx[i, idy] <- 1
}
overlapMtxPfct <- matrix(data = 0, nrow = length(gr_amplicons), ncol = length(gr_amplicons))
for(i in seq_along(gr_amplicons)){
  idy <- findOverlaps(gr_amplicons[i], gr_amplicons, type = "equal") %>% subjectHits()
  overlapMtxPfct[i, idy] <- 1
}

#heatmap
fh <- function(x) hclust(dist(x), method = "ward.D2")
col_fun1 <- colorRamp2(c(0,1), c("lightgray", "#89288F"))
chr_colors <- c(ArchR::ArchRPalettes$stallion2, ArchR::ArchRPalettes$calm)[c(1:24)]
names(chr_colors) <- paste0("chr", c(1:22,"X","Y"))

c("ecDNA" = "#E7298A", "BFB"= "#7570B3", "Complex-non-cyclic" = "#D95F02", "Linear" = "#1B9E77")

ha1 <- HeatmapAnnotation(class = gr_amplicons$Classification, chr = factor(seqnames(gr_amplicons), levels = paste0("chr", c(1:22,"X","Y"))),
                         col = list(class = c("ecDNA" = "#E7298A", "BFB"= "#7570B3", "Complex-non-cyclic" = "#D95F02", "Linear" = "#1B9E77"), chr = chr_colors))
ra1 <- rowAnnotation(class = gr_amplicons$Classification, chr = factor(seqnames(gr_amplicons), levels = paste0("chr", c(1:22,"X","Y"))),
                     col = list(class = c("ecDNA" = "#E7298A", "BFB"= "#7570B3", "Complex-non-cyclic" = "#D95F02", "Linear" = "#1B9E77"), chr = chr_colors))

ht1 <- Heatmap(overlapMtx, name = "overlap", cluster_columns = fh, cluster_rows = fh, col = col_fun1,
               top_annotation = ha1, left_annotation = ra1, column_split = gr_amplicons$Classification, row_split = gr_amplicons$Classification,
               row_title = "352 amplified intervals")
p1 <- draw(ht1)

ht2 <- Heatmap(overlapMtx, name = "overlap", cluster_columns = fh, cluster_rows = fh, col = col_fun1,
               top_annotation = ha1, left_annotation = ra1, row_title = "352 amplified intervals")
p2 <- draw(ht2)

ht3 <- Heatmap(overlapMtxPfct, name = "overlap, perfect-match", cluster_columns = fh, cluster_rows = fh, col = col_fun1,
               top_annotation = ha1, left_annotation = ra1, column_split = gr_amplicons$Classification, row_split = gr_amplicons$Classification,
               row_title = "352 amplified intervals")
p3 <- draw(ht3)

pdf("output/Plots/01_Heatmap_AmplifiedIntervals_Overlap_Normal.pdf", width = 9, height = 8)
p1
p2
p3
dev.off()

#karyoplot
library(karyoploteR)
pdf("output/Plots/01_KaryotypePlot_AmplifiedIntervals_Normal.pdf", width = 6, height = 8)
kp <- plotKaryotype(genome = "hg38", chromosomes = paste0("chr", c(1:22, "X", "Y")))
kpPlotRegions(kp, data=gr_amplicons, col = "black")
dev.off()

saveRDS(gr_amplicons, "rds/01_gr_amplicon_normal.rds")
