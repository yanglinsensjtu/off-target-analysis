
library(Gviz)
library(biomaRt)

plotGviz <- function(){
  
# annotation track creation -----------------------------------------------

anntrack  <- AnnotationTrack(start = stt,
                             end = ed,
                             strand = st, 
                             chromosome = chr, 
                             genome ="hg19", 
                             name ="CRISPR")

# biomartGeneRegiontrack creation -----------------------------------------


f <- as.integer(stt - 60000)
t <- as.integer(stt + 60000)
bm <- useMart(host="grch37.ensembl.org", 
              biomart="ENSEMBL_MART_ENSEMBL", 
              dataset="hsapiens_gene_ensembl")
biomTrack <- BiomartGeneRegionTrack(genome="hg19", 
                                    chromosome=chr, 
                                    start=f, 
                                    end=t,
                                    name="ENSEMBL", 
                                    biomart=bm)

# Ideogramtrack creation --------------------------------------------------

itrack <- IdeogramTrack(genome = 'hg19',chromosome = chr)

# genome axis track creation ----------------------------------------------

gatrack <- GenomeAxisTrack(distFromAxis = 15,labelPos="below")

# save the image ----------------------------------------------------------
oldpath <- getwd()
setwd(dir = '../')
jpeg(str_c(i,"output.jpg"), width = 5000, height =3090,res = 720)
plotTracks(trackList = ls,
           transcriptAnnotation="symbol",
           sizes = sz,
           from = f,
           to = t,
           main = mn,
           cex.main = 0.8)
dev.off()
setwd(dir = oldpath)

}