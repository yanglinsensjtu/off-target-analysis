source('offtargetanalysis/sgRNA.alignment.R')
source('offtargetanalysis/data.tidy.R')
source('offtargetanalysis/target.anotation.R')
source('offtargetanalysis/sgRNA.alignment.R')
source('offtargetanalysis/write.alignment.R')
source('offtargetanalysis/plot.function.R')
source('offtargetanalysis/changeIRanges.R')
sgRNA <- readDNAStringSet('../sgRNA sequence.fasta')
sgRNA.alignment(cco.ann,sgRNA)
write.alignment2txt(cco.ann,
                    sgRNA,
                    path = '../align sgRNA target.txt')

# Artificial analysis of the most likely off target sites -----------------

cco.ann.co <- cco.ann[c(1,3,6,8,16,18,26)]

write.alignment2txt(cco.ann.co, 
                    sgRNA,
                    path = '../screen sgRNA off target.txt')

# plot the off target sites -----------------------------------------------

source('offtargetanalysis/plot.main.R')



# get the off target site sequences ---------------------------------------

cco.ann.co.extend <- changeIRanges(granges.obj = cco.ann.co,
                           upstream = 300, 
                           width = 600)

cco.ann.co.extend.seq <- getSeq(BS.hg19, cco.ann.co.extend)
names(cco.ann.co.extend.seq) <- cco.ann.co.extend$GENEID
names(cco.ann.co.extend.seq)[6] <- '144535'
writeXStringSet(cco.ann.co.extend.seq, 
                filepath = '../predict off target genes sequences')


