source('data-tidy.R')
source('plot.function.R')

coding.obj <- cosmid.gr.ann.c


for (i in seq_along(coding.obj)) {
  plotGviz(coding.obj = coding.obj,
           i=i,
           bounds = 10000,
           title = '',
           path = '../cosmid/')
}

coding.obj <- cas.offfinder.gr.ann.c


for (i in seq_along(coding.obj)) {
  plotGviz(coding.obj = coding.obj,
           i=i,
           bounds = 10000,
           title = '',
           path = '../casofffinder/')
}