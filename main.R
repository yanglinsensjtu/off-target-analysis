source('data-tidy.R')
source('plot.function.R')
source('wgs.R')
source('sgRNA alignment.R')
source('changeIRanges.R')

coding.obj <- cosmid.gr.ann.c.check.back


for (i in seq_along(coding.obj)) {
  
  plotGviz(coding.obj = coding.obj,
           i=i,
           bounds = 10000,
           title = str_c(coding.obj$Score[i],
                         '-',
                         coding.obj$Query.type[i],
                         '-',
                         coding.obj$Mismatch[i]),
           path = '../cosmid/')
}

coding.obj <- cas.offfinder.gr.ann.c.check.back


for (i in seq_along(coding.obj)) {
  plotGviz(coding.obj = coding.obj,
           i=i,
           bounds = 20000,
           title = cas.offfinder.gr.ann.c.check.back$DNA[i],
           path = '../casofffinder/')
}

wgs <- c(wgs.less100.gr,wgs.more100.gr)
findOverlaps(cosmid.gr, wgs)
findOverlaps(cas.offfinder.gr, wgs)

findOverlaps(cosmid.gr.ann, wgs)
findOverlaps(cas.offfinder.gr.ann, wgs)

a <- findOverlaps(cosmid.gr.ann.c.check.back, cas.offfinder.gr.ann.c.check.back)
offtarget.gr <- c(cosmid.gr, cas.offfinder.gr, off.spotter.gr)
wgs.offtarget <- findOverlaps(wgs.jc, offtarget.gr)
wgs.offtarget.joint.gr <- unique(wgs.jc[wgs.offtarget@from])

for (i in seq_along(wgs.offtarget.joint.gr)) {
  plotGviz(coding.obj = wgs.offtarget.joint.gr,
           i=i,
           bounds = 20000,
           title = '',
           path = '../wgs-offtarget-common/')
}
sgRNA.alignment(grange = wgs.offtarget.joint.gr,
                sgRNA = targetgene)
b <- findOverlaps(wgs.old.jc, wgs.offtarget.joint.gr)

wgs.and.pre.sites <- wgs.offtarget.joint.gr[b@to]
wgs.and.pre.sites.ex <- changeIRanges(granges.obj = wgs.and.pre.sites,
                            upstream = 200,
                            width = 400)
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

for (i in seq_along(predict.off.target.coding.site.ann)) {
  plotGviz(coding.obj = predict.off.target.coding.site.ann,
           i=i,
           bounds = 20000,
           title = '',
           path = '../offtarget predict site common/')
}
for (i in seq_along(predict.off.target.coding.site)) {
  plotGviz(coding.obj = predict.off.target.coding.site,
           i=i,
           bounds = 20000,
           title = '',
           path = '../offtarget predict site common/2/')
}
for (i in seq_along(predict.off.target.coding.site.ann.screen)) {
  plotGviz(coding.obj = predict.off.target.coding.site.ann.screen,
           i=i,
           bounds = 20000,
           title = '',
           path = '../offtarget predict site common/3/')
}

sequences <- changeIRanges(granges.obj = predict.off.target.coding.site.ann.screen,
                           upstream = 300, 
                           width = 600)
                     
sequences.DNA <- toString(getSeq(BS.hg19, sequences))
oldpath <- getwd()
setwd(dir = '../')
write(sequences.DNA, 'predict off target genes sequences.txt')
setwd(dir = oldpath)


