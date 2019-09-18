library(Rsamtools)


which <- Granges2IRangeslist(wgs.and.pre.sites.ex)
which1 <- IRangesList(chr1 = IRanges(111000, 111023))
what <- c("rname", "strand", "pos", "qwidth", "seq")
param <- ScanBamParam(which = which1, what = what)
file <-  'F:/WGS/T0501 T0502 WGS/s1382x10001JointCall/Alignment/BAM/R19026468LD01-T0501-MOCK_sorted_dedup_realign.bam'
index <- 'F:/WGS/T0501 T0502 WGS/s1382x10001JointCall/Alignment/BAM/R19026468LD01-T0501-MOCK_sorted_dedup_realign.bam.bai'
bam <- scanBam(file = file,
               index = index,
               param = param)

bam
library(BSgenome.Hsapiens.UCSC.hg19)
strack <- SequenceTrack(Hsapiens)
library(Gviz)
altrack <- AlignmentsTrack(file, 
                           start = 111000, 
                           end = 111043, 
                           chromosome = 'chr1', 
                           genome = 'hg19',
                           isPaired = T
                             )
plotTracks(c(altrack,strack), chromosome = 'chr1', from = 111000, to = 111043)
