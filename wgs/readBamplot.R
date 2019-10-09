library(Rsamtools)
library(stringr)
 

num <- '10'
chr <- 'chr'
chrom <- str_c(chr,num)
#which <- Granges2IRangeslist(wgs.and.pre.sites.ex)
which1 <- IRangesList(chr1 = IRanges(20000, 22000))
what <- c("rname", "strand", "pos", "qwidth", "seq")
param <- ScanBamParam(which = which1, what = what)
file <-  str_c('../BAM/',chr,num,'con.bam')
a <- str_c('../BAM/',chr,num,'.bam.bai')
bam <- scanBam(file = file,
               index = a,
               param = param)

bam
library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg19)
strack <- SequenceTrack(Hsapiens)

altrack <- AlignmentsTrack(file, 
                           start = 111000, 
                           end = 111067, 
                           chromosome = chrom, 
                           genome = 'hg19',
                           isPaired = T
                             )
plotTracks(c(altrack,strack), chromosome = chrom, from = 111000, to = 111067)
