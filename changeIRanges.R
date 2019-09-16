
# change granges obj ranges -----------------------------------------------

changeIRanges <- function(granges.obj = granges.obj, upstream = integer, width = integer){
  GRanges(seqnames = granges.obj@seqnames,
          IRanges(start = start(granges.obj@ranges) - upstream, 
                  width = width),
          strand = granges.obj@strand,
          seqinfo = granges.obj@seqinfo,
          mcols(granges.obj))
}