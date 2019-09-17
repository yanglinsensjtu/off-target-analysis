library(Rsamtools)
predict.off.target.coding.site.ann.screen <- changeIRanges(granges.obj = predict.off.target.coding.site.ann.screen,
                                                           upstream = 1000,
                                                           width = 1000)


which <- IRangesList(predict.off.target.coding.site.ann.screen)
which <- IRangesList(chr1  = IRanges(10000,20000))
what <- c("rname", "strand", "pos", "qwidth", "seq")
param <- ScanBamParam(which = which, what = what)
file = 'F:/WGS/T0501 T0502 WGS/s1382x10001Somatic/Alignment/BAM/R19026469LD01-T0502-CAR_sorted_dedup_realign.bam'
bam <- scanBam(file = file,
               param = param)

bam 
