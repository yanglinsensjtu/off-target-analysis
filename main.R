source('data-tidy.R')
source('plot.function.R')
source('wgs.R')

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
