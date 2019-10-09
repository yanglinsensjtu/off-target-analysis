source('offtargetanalysis/sgRNA.alignment.R')
source('offtargetanalysis/data.tidy.R')
source('offtargetanalysis/target.anotation.R')
source('offtargetanalysis/sgRNA.alignment.R')
source('offtargetanalysis/write.alignment.R')
sgRNA <- readDNAStringSet('../sgRNA sequence.fasta')
sgRNA.alignment(cco.ann,sgRNA)
write.alignment2txt(cco.ann,sgRNA)




