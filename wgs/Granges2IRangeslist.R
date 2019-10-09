library(GenomicRanges)

gr.obj <- GRanges(seqnames = c('chr1','chr4','chr3'),
                  IRanges(start = c(100, 200, 300),
                          width = 2000))

Granges2IRangeslist <- function(gr.obj){
  irangeslist <- IRangesList()
  name <- vector('character',length(gr.obj))
  for (i in seq_len(length(gr.obj))) {
    name[i] <- as.character(gr.obj@seqnames[i])
    irangeslist[i] <-  IRanges(start(gr.obj@ranges[i]),end(gr.obj@ranges[i]))
  }
  names(irangeslist) <- name
  irangeslist
}

