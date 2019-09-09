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
cas.offfinder.gr.ann  <- locateVariants(cas.offfinder.gr, txdb, AllVariants())
off.spotter.gr.ann <- locateVariants(off.spotter.gr, txdb, AllVariants())
cas.offfinder.gr.ann.c <- cas.offfinder.gr.ann[cas.offfinder.gr.ann$LOCATION == 'coding']

# gene location visualization ---------------------------------------------

  
  
  
anntrack  <- AnnotationTrack(start = cas.offfinder.gr.ann.c,
                             width = wd,
                             strand=st, 
                             chromosome=chr, 
                             genome="hg19", 
                             name="Crispr")

library(biomaRt)
f <- as.integer(stt - 30000)
t <- as.integer(stt + 30000)
bm <- useMart(host="grch37.ensembl.org", 
              biomart="ENSEMBL_MART_ENSEMBL", 
              dataset="hsapiens_gene_ensembl")
biomTrack <- BiomartGeneRegionTrack(genome="hg19", 
                                    chromosome=chr, 
                                    start=f, 
                                    end=t,
                                    name="ENSEMBL", 
                                    biomart=bm)

itrack <- IdeogramTrack(genome = 'hg19',chromosome = chr)

gatrack <- GenomeAxisTrack(distFromAxis = 15,labelPos="below")

jpeg(str_c(i,"output.jpg"), width = 5000, height =3090,res = 720)
plotTracks(list(itrack,gatrack,biomTrack,anntrack),
           transcriptAnnotation="symbol",
           sizes = c(1,1,5,1),
           from = f,
           to = t,
           main = str_c(i,'-',refseq,'  ->  ',altseq),
           cex.main = 0.8)
dev.off()






