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
off.spotter$V1 <- str_c('chr',off.spotter$V1)
names(off.spotter)[5:6] <- c('query','hit')
off.spotter.gr <- GRanges(seqnames = off.spotter$V1,
                          strand = off.spotter$V2,
                          ranges = IRanges(start = off.spotter$V3,
                                           end = off.spotter$V4),
                          dplyr::select(off.spotter, 
                                        query, 
                                        hit))

cas.offfinder.gr <- GRanges(seqnames = cas.offfinder$Chromosome,
                            strand = cas.offfinder$Direction,
                            ranges = IRanges(start = cas.offfinder$Position,
                                             width = 22),
                            dplyr::select(cas.offfinder,
                                          X.Bulge.type,
                                          crRNA,
                                          DNA,
                                          Mismatches,
                                          Bulge.Size))

cosmid <- separate(cosmid, `Chr Position`, sep = ':', into = c('chromosome', 'position'))
cosmid <- separate(cosmid, position, sep = '-', into = c('position.start', 'position.end'))
names(cosmid) <- str_replace_all(names(cosmid), ' ', '.')
cosmid.gr <- GRanges(seqnames = cosmid$chromosome,
                     strand = str_sub(cosmid$Strand,1,1),
                     ranges = IRanges(start = as.integer(cosmid$position.start),
                                      end = as.integer(cosmid$position.end)),
                     dplyr::select(cosmid,
                                   Score,
                                   Query.tag,
                                   Search.result,
                                   Query.type,
                                   Mismatch,
                                   Ends.with.RG,
                                   Cut.site))












