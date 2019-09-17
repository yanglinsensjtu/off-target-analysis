
gr.obj <- predict.off.target.coding.site.ann.screen
irangeslist <- vector('list',length(gr.obj))
for (i in seq_len(length(gr.obj))) {
  chr <- deparse(substitute(as.character(gr.obj@seqnames[i])))
  temp <- IRangesList(chr = IRanges(start(gr.obj@ranges[i]),end(gr.obj@ranges[i])))
  temp
  irangeslist[i]  <-  temp
}
irangeslist
