#-------------------------------------------------------------------------------
# 04_prepareRNAdata.R
#-------------------------------------------------------------------------------
library(data.table)
library(parallel)
library(dplyr)
library(GenomicRanges)
library(GenomicFeatures)
library(SummarizedExperiment)
theme_cb <- function(color="black"){theme_classic() + theme(axis.text = element_text(color=color), axis.ticks = element_line(color=color))}
SampleInfo <- read.csv("output/Tables/00_SampleListSummary_Curated.csv", row.names = 1)

#----- generate count matrix (STAR counts) -----#
rna_info <- SampleInfo$UUID_R
names(rna_info) <- SampleInfo$StudyID

counts_ls <- mclapply(as.character(rna_info), function(i){
  tmp <- paste0("../data/rna_count/", i, ".ReadsPerGene.out.tab")
  out <- fread(tmp)
  out <- data.frame(gene = out$V1, count = out$V2)[-c(1:4),]
  colnames(out) <- c("gene", i)
  return(out)
}, mc.cores = 26)
counts <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "gene", all = TRUE), counts_ls)

#checking data
counts[c(1:5), c(1:5)]
which(is.na(counts))

rownames(counts) <- counts$gene
counts <- counts[,-1]

tail(rownames(counts)) # c(59583:59586)
# [1] "ENSG00000288696.1" "ENSG00000288697.1" "IGH-.g@-ext"       "IGH.g@-ext"        "IGL-.g@-ext"       "IGL.g@-ext"       rownames(counts)[c(59583:59586)]
counts <- counts[-c(59583:59586),]
saveRDS(counts, "rds/04_RNA_rawcounts.rds")

#----- colData & rowData -----#
#colData
SampleInfo2 <- SampleInfo[names(rna_info),]

#rowData
gff_anno <- rtracklayer::import.gff("ref/gencode.v37.annotation.gtf")
# gff_anno$gene_type %>% table
# gff_anno$gene_name[which(seqnames(gff_anno) == "chrM")]

geneData <- DataFrame(row.names = gff_anno$gene_id, 
                      gene_type = gff_anno$gene_type, 
                      gene_name = gff_anno$gene_name)
geneData <- geneData[!duplicated(rownames(geneData)),]
geneData <- geneData[rownames(counts),]

#----- calculate TPM -----#
#gene length: https://www.biostars.org/p/83901/
txdb <- makeTxDbFromGFF("ref/gencode.v37.annotation.gtf",format="gtf")
exons.list.per.gene <- exonsBy(txdb,by="gene") # then collect the exons per gene id
exonic.gene.sizes <- sum(width(reduce(exons.list.per.gene))) # then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then

intersect(names(exonic.gene.sizes), rownames(counts)) %>% length #59582
len <- exonic.gene.sizes[rownames(counts)] %>% as.numeric()

countToTpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}
tpms <- countToTpm(counts, len)

#----- finishing RNA-seq SummarizedExperiment -----#
geneData2 <- genes(txdb)[rownames(geneData)]
mcols(geneData2) <- geneData

se <- SummarizedExperiment(assays = list(counts = counts, log2tpm = log2(tpms+1)),
                           colData = DataFrame(SampleInfo2),
                           rowData = geneData2)
saveRDS(se, "rds/04_se.rds")
