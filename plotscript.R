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
wgs.offtarget.joint.gr <- unique(wgs.jc[wgs.offtarget@from])
for (i in seq_along(wgs.offtarget.joint.gr)) {
  plotGviz(coding.obj = wgs.offtarget.joint.gr,
           i=i,
           bounds = 20000,
           title = '',
           path = '../wgs-offtarget-common/')
}

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