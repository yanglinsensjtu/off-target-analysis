source('data-tidy.R')
source('plot.function.R')
source('wgs.R')
source('sgRNA alignment.R')
source('changeIRanges.R')
source('Granges2IRangeslist.R')


wgs <- c(wgs.less100.gr,wgs.more100.gr)
findOverlaps(cosmid.gr, wgs)
findOverlaps(cas.offfinder.gr, wgs)

findOverlaps(cosmid.gr.ann, wgs)
findOverlaps(cas.offfinder.gr.ann, wgs)

a <- findOverlaps(cosmid.gr.ann.c.check.back, cas.offfinder.gr.ann.c.check.back)
offtarget.gr <- c(cosmid.gr, cas.offfinder.gr, off.spotter.gr)
wgs.offtarget <- findOverlaps(wgs.jc, offtarget.gr)
offtarget.gr.wgs <- unique(offtarget.gr[wgs.offtarget@to])
sgRNA.alignment(grange = offtarget.gr.wgs,
                sgRNA = targetgene)#offtarget predict sgRNA align with the target sgRNA


sgRNA.alignment(grange = wgs.offtarget.joint.gr,
                sgRNA = targetgene)
b <- findOverlaps(wgs.old.jc, wgs.offtarget.joint.gr)

wgs.and.pre.sites <- wgs.offtarget.joint.gr[b@to]
wgs.and.pre.sites.ex <- changeIRanges(granges.obj = wgs.and.pre.sites,
                            upstream = 300,
                            width = 600)
wgs.and.pre.sites.ex.seq <- getSeq(BS.hg19, wgs.and.pre.sites.ex)
wgs.and.pre.sites.ex.seq.str <- toString(wgs.and.pre.sites.ex.seq)
oldpath <- getwd()
setwd(dir = '../')
write(wgs.and.pre.sites.ex.seq.str, 'wgs.and.pre.sites.ex.seq.str.txt')
setwd(dir = oldpath)

# off target predict site -------------------------------------------------

predict.off.target.coding.site <- union(cas.offfinder.gr.ann.c, cosmid.gr.ann.c)

sgRNA.alignment(grange = predict.off.target.coding.site.ann,
                sgRNA = targetgene)



sequences <- changeIRanges(granges.obj = predict.off.target.coding.site.ann.screen,
                           upstream = 300, 
                           width = 600)
                     
sequences.DNA <- toString(getSeq(BS.hg19, sequences))
oldpath <- getwd()
setwd(dir = '../')
write(sequences.DNA, 'predict off target genes sequences.txt')
setwd(dir = oldpath)

source('readBamplot.R')

