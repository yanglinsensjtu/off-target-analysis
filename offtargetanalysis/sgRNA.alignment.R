library(BSgenome.Hsapiens.UCSC.hg19)
BS.hg19 <- BSgenome.Hsapiens.UCSC.hg19
# sgRNA aliganment --------------------------------------------------------
sgRNA.alignment <- function(grange = Grange.obj,
                            sgRNA = targetgene){
  off.target.seq <- getSeq(BS.hg19, unique(grange))
  if (is.null(grange$GENEID)) {
    print('There is no geneid')
    for (i in seq_len(length(unique(grange)))) {
      print(i)
      print(pairwiseAlignment(sgRNA, 
                              off.target.seq[i],
                              gapOpening = 0,
                              gapExtension = 1))
    }
  }else{
    for (i in seq_len(length(unique(grange)))) {
      print(i)
      
        if (is.na(grange$GENEID[i])) {
          print('The geneid is NA')
        }else{
          print(grange$GENEID[i])
        }
      
      print(pairwiseAlignment(sgRNA, 
                              off.target.seq[i],
                              gapOpening = 0,
                            gapExtension = 1))
    }
  }
}

