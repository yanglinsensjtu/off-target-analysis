library(BSgenome.Hsapiens.UCSC.hg19)
BS.hg19 <- BSgenome.Hsapiens.UCSC.hg19
getGranges <- function(chr = chr,
                       sequence = p1){
  p1 <- DNAString(sequence)
  chrseq <-  BS.hg19[[chr]]
  m1 <- matchPattern(p1, chrseq, max.mismatch = 0)
  if (!setequal(nchar(m1), integer(0))) {
    m1.gr <- GRanges(seqnames = chr,
                     ranges = m1@ranges,
                     strand = '+')
    return(m1.gr)
  }else{
    p1.R <- reverseComplement(p1)
    m.R <- matchPattern(p1.R, chrseq, max.mismatch = 0)
    m.R.gr <- GRanges(seqnames = chr,
                      ranges = m.R@ranges,
                      strand = '-')
    return(m.R.gr)
  }
}


