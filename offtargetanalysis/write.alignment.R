library(Biostrings)
library(readr)
library(BSgenome.Hsapiens.UCSC.hg19)
BS.hg19 <- BSgenome.Hsapiens.UCSC.hg19
write.alignment2txt <- function(granges.ann = granges.ann,
                                sgRNA = sgRNA){
  off.target.seq <- getSeq(BS.hg19, granges.ann)
  write_lines('****The sgRNA alignment with the off target sites***', 
              path = '../align sgRNA target.txt')
  for (i in seq_len(length(granges.ann))) {
    cat(paste(i, ' '))
    paln <- pairwiseAlignment(sgRNA, 
                              off.target.seq[i],
                              gapOpening = 0,
                              gapExtension = 1)
    write_lines('------------------------', 
                path = '../align sgRNA target.txt',
                append = T)
    write_lines(names(sgRNA), 
                path = '../align sgRNA target.txt',
                append = T)
    write_lines(paste(granges.ann$GENEID[i], 
                      as.character(granges(granges.ann)[i])), 
                path = '../align sgRNA target.txt',
                append = T)
    write_lines('.....', 
                path = '../align sgRNA target.txt',
                append = T)
    write_lines(alignedPattern(paln), 
                path = '../align sgRNA target.txt',
                append = T)
    write_lines(alignedSubject(paln), 
                path = '../align sgRNA target.txt',
                append = T)
    
  }
}
