library(Gviz)
library(biomaRt)
library(stringr)

plotGviz <- function(coding.obj, i, bounds, plottitle = '', folder){
  tempfolder <- file.path('..', folder)
  if (!file.exists(tempfolder)) {
    print(paste('The upper directory does not exist', folder))
    print('Creating folder!')
    dir.create(tempfolder)
  }
  chr <- as.character(coding.obj@seqnames[i])
  strand <- as.character(coding.obj@strand[i])
  start <- IRanges::start(coding.obj)[i]
  end <- IRanges::end(coding.obj)[i]
  filetitle <- str_c(i,chr,'-',start,'-',end,'.jpg')
  print(filetitle)
# annotation track creation -----------------------------------------------

  anntrack  <- AnnotationTrack(start = start,
                               end = end,
                               strand = strand, 
                               chromosome = chr, 
                               genome ="hg19", 
                               name ="CRISPR")

# biomartGeneRegiontrack creation -----------------------------------------

  f <- as.integer(start - bounds)
  t <- as.integer(start + bounds)
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
  setwd(dir = tempfolder)
  sz <- c(1,1,4,1)
  ls <- list(itrack,gatrack,biomTrack,anntrack)
  mn <- str_c(i,'-',chr,'...',plottitle)
  jpeg(filetitle, width = 5000, height =3090,res = 720)
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

