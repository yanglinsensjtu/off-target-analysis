library(BSgenome.Hsapiens.UCSC.hg19)
BS.hg19 <- BSgenome.Hsapiens.UCSC.hg19
# sgRNA aliganment --------------------------------------------------------
targetgene <- DNAString(read_lines('../sgRNAseq.txt'))
sgRNA.alignment <- function(grange = Grange.obj,
                            sgRNA = targetgene){
  off.target.seq <- getSeq(BS.hg19, unique(grange))
  for (i in seq_len(length(unique(grange)))) {
    print(i)
    print(pairwiseAlignment(sgRNA, 
                            off.target.seq[i],
                            gapOpening = 0,
                            gapExtension = 1))
  }
}