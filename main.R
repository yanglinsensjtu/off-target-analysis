source('data-tidy.R')
source('plot.function.R')

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

findOverlaps(cosmid.gr, wgs.less100.gr)
findOverlaps(cosmid.gr,wgs.more100.gr)
findOverlaps(cas.offfinder.gr, wgs.less100.gr)
findOverlaps(cas.offfinder.gr, wgs.more100.gr)


