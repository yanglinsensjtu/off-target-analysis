# gene annotation ---------------------------------------------------------

library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

cosmid.gr.ann <- locateVariants(cosmid.gr, txdb, AllVariants())
table(unique(cosmid.gr.ann)$LOCATION)
cosmid.gr.ann.c <- cosmid.gr.ann[cosmid.gr.ann$LOCATION == 'coding']
cosmid.gr.ann.c <- unique(cosmid.gr.ann.c)
cosmid.gr.ann.c.check.back<- cosmid.gr[findOverlaps(cosmid.gr,cosmid.gr.ann.c)@from,]

off.spotter.gr.ann <- locateVariants(off.spotter.gr, txdb, AllVariants())
table(unique(off.spotter.gr.ann)$LOCATION)
off.spotter.gr.ann.c <- off.spotter.gr.ann[off.spotter.gr.ann$LOCATION == 'coding']

cas.offfinder.gr.ann  <- locateVariants(cas.offfinder.gr, txdb, AllVariants())
table(unique(cas.offfinder.gr.ann)$LOCATION)
cas.offfinder.gr.ann.c <- cas.offfinder.gr.ann[cas.offfinder.gr.ann$LOCATION == 'coding']
cas.offfinder.gr.ann.c <- unique(cas.offfinder.gr.ann.c)
cas.offfinder.gr.ann.c.check.back <- cas.offfinder.gr[findOverlaps(cas.offfinder.gr,cas.offfinder.gr.ann.c)@from,]
cas.offfinder.gr.ann.c.check.back <- unique(cas.offfinder.gr.ann.c.check.back)

cc <- union(cas.offfinder.gr.ann.c, cosmid.gr.ann.c)
cco <- union(cc, off.spotter.gr.ann.c) #The union of the three
cco.ann <- unique(locateVariants(cco, txdb, AllVariants()))


