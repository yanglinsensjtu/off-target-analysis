changeIRanges <- function(granges = granges, upstream = integer, width = integer){
  GRanges(seqnames = granges@seqnames,
          IRanges(start = start(granges@ranges) - upstream, 
                  width = width),
          strand = granges@strand,
          seqinfo = granges@seqinfo,
          mcols(granges))
}