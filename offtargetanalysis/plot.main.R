coding.obj <- cco.ann.co
for (i in seq_along(coding.obj)) {
  chr <- as.character(coding.obj@seqnames[i])
  strand <- as.character(coding.obj@strand[i])
  start <- IRanges::start(coding.obj)[i]
  end <- IRanges::end(coding.obj)[i]
  filetitle <- str_c(i,chr,'-',start,'-',end,'.jpg')
  if (file.exists(paste0('../off target map/',filetitle))) {
    print(paste(filetitle,'exists'))
    next()
  }else{
    plotGviz(coding.obj = coding.obj,
             i=i,
             bounds = 10000,
             plottitle = str_c(coding.obj$Score[i],
                               '-',
                               coding.obj$Query.type[i],
                               '-',
                               coding.obj$Mismatch[i]),
             folder = 'off target map')
  }
  
}
