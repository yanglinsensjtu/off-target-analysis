# read the raw data into the R --------------------------------------------
library(readr)
library(tibble)
library(readxl)
library(dplyr)
off.spotter<- read.delim('../Off-Spotter.hg19.txt', 
                         header = F, 
                         stringsAsFactors = F) %>% 
  as_tibble()

write.csv(off.spotter, 
          file = '../Off-Spotter.csv')
cas.offfinder <- read.delim('../Cas-OFFinder.hg19.txt',
                            header = T,
                            stringsAsFactors = F) %>% 
  as_tibble() 

cosmid <- read_xlsx('../cosmid_hg19_qz5depo5i9_output.xlsx')

cosmid.t <- read.delim('../comsid hg19.txt',
                       header = F, 
                       stringsAsFactors = F) %>% 
  as_tibble()


# transform the data into granges obj -------------------------------------

library(GenomicRanges)
library(stringr)
library(tidyr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqinfo <- Seqinfo(genome = 'hg19')
off.spotter$V1 <- str_c('chr',off.spotter$V1)
names(off.spotter)[5:6] <- c('query','hit')
off.spotter.gr <- GRanges(seqnames = off.spotter$V1,
                          strand = off.spotter$V2,
                          ranges = IRanges(start = off.spotter$V3,
                                           end = off.spotter$V4),
                          seqinfo = seqinfo,
                          dplyr::select(off.spotter, 
                                        query, 
                                        hit))

cas.offfinder.gr <- GRanges(seqnames = cas.offfinder$Chromosome,
                            strand = cas.offfinder$Direction,
                            ranges = IRanges(start = cas.offfinder$Position,
                                             width = 23),
                            seqinfo = seqinfo,
                            dplyr::select(cas.offfinder,
                                          X.Bulge.type,
                                          crRNA,
                                          DNA,
                                          Mismatches,
                                          Bulge.Size))

cosmid <- separate(cosmid, `Chr Position`, sep = ':', into = c('chromosome', 'position'))
cosmid$chromosome <- stringr::str_replace_all(cosmid$chromosome, 'Chr', 'chr')
cosmid <- separate(cosmid, position, sep = '-', into = c('position.start', 'position.end'))
names(cosmid) <- str_replace_all(names(cosmid), ' ', '.')
cosmid.gr <- GRanges(seqnames = cosmid$chromosome,
                     strand = str_sub(cosmid$Strand,1,1),
                     ranges = IRanges(start = as.integer(cosmid$position.start),
                                      end = as.integer(cosmid$position.end)),
                     seqinfo = seqinfo,
                     dplyr::select(cosmid,
                                   Score,
                                   Query.tag,
                                   Search.result,
                                   Query.type,
                                   Mismatch,
                                   Ends.with.RG,
                                   Cut.site))


# gene annotation ---------------------------------------------------------

library(VariantAnnotation)
cosmid.gr.ann <- locateVariants(cosmid.gr, txdb, AllVariants())
table(cosmid.gr.ann$LOCATION)
cas.offfinder.gr.ann  <- locateVariants(cas.offfinder.gr, txdb, AllVariants())
table(cas.offfinder.gr.ann$LOCATION)
off.spotter.gr.ann <- locateVariants(off.spotter.gr, txdb, AllVariants())
table(off.spotter.gr.ann$LOCATION)
cas.offfinder.gr.ann.c <- cas.offfinder.gr.ann[cas.offfinder.gr.ann$LOCATION == 'coding']

# pairwiseAlignment -------------------------------------------------------
i <- 3
library(BSgenome.Hsapiens.UCSC.hg19)
BS.hg19 <- BSgenome.Hsapiens.UCSC.hg19
off.target.seq <- getSeq(BS.hg19, cas.offfinder.gr.ann.c)

TRAC <- DNAString(read_lines('../sgRNAseq.txt'))
a <- pairwiseAlignment(TRAC, off.target.seq[i])
print(a)

# gene location visualization ---------------------------------------------

chr <- as.character(cas.offfinder.gr.ann.c@seqnames[i])
st <- as.character(cas.offfinder.gr.ann.c@strand[i])
stt <-  IRanges::start(cas.offfinder.gr.ann.c)[i]
ed <- IRanges::end(cas.offfinder.gr.ann.c)[i]

# annotation track creation -----------------------------------------------

library(Gviz)
anntrack  <- AnnotationTrack(start = stt,
                             end = ed,
                             strand = st, 
                             chromosome = chr, 
                             genome ="hg19", 
                             name ="CRISPR")




# biomartGeneRegiontrack creation -----------------------------------------


library(biomaRt)
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



# plot the jpg ------------------------------------------------------------


sz <- c(1,1,4,1)
ls <- list(itrack,gatrack,biomTrack,anntrack)
mn <- str_c(i,a)

plotTracks(trackList = ls,
           transcriptAnnotation="symbol",
           sizes = sz,
           from = f,
           to = t,
           main = mn,
           cex.main = 0.8)

# save the image ----------------------------------------------------------

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




