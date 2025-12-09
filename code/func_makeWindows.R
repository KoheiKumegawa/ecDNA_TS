makeWindows <- function(genome, blacklist, windowSize = 10e6, slidingSize = 2e6){
  chromSizes <- GRanges(names(seqlengths(genome)), IRanges(1, seqlengths(genome)))
  chromSizes <- GenomeInfoDb::keepStandardChromosomes(chromSizes, pruning.mode = "coarse")
  windows <- slidingWindows(x = chromSizes, width = windowSize, step = slidingSize) %>% unlist %>% .[which(width(.)==windowSize),]
  mcols(windows)$wSeq <- as.character(seqnames(windows))
  mcols(windows)$wStart <- start(windows)
  mcols(windows)$wEnd <- end(windows)
  message("Subtracting Blacklist...")
  windowsBL <- lapply(seq_along(windows), function(x){
    if(x %% 100 == 0){
      message(sprintf("%s of %s", x, length(windows)))
    }
    gr <- GenomicRanges::setdiff(windows[x,], blacklist)
    mcols(gr) <- mcols(windows[x,])
    return(gr)
  })
  names(windowsBL) <- paste0("w",seq_along(windowsBL))
  windowsBL <- unlist(GenomicRangesList(windowsBL), use.names = TRUE)
  mcols(windowsBL)$name <- names(windowsBL)
  message("Adding Nucleotide Information...")
  windowSplit <- split(windowsBL, as.character(seqnames(windowsBL)))
  windowNuc <- lapply(seq_along(windowSplit), function(x){
    message(sprintf("%s of %s", x, length(windowSplit)))
    chrSeq <- Biostrings::getSeq(genome,chromSizes[which(levels(seqnames(chromSizes))==names(windowSplit)[x])])
    grx <- windowSplit[[x]]
    aFreq <- alphabetFrequency(Biostrings::Views(chrSeq[[1]], ranges(grx)))
    mcols(grx)$GC <- rowSums(aFreq[, c("G","C")]) / rowSums(aFreq)
    mcols(grx)$AT <- rowSums(aFreq[, c("A","T")]) / rowSums(aFreq)
    return(grx)
  }) %>% GenomicRangesList %>% unlist %>% sortSeqlevels %>% sort
  windowNuc$N <- 1 - (windowNuc$GC + windowNuc$AT)
  windowNuc
}
