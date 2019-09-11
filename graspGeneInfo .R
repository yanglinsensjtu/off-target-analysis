
# Loader package ----------------------------------------------------------

if (!require(httr)) {
  install.packages('httr')
  library(httr)
} 

if (!require(rvest)) {
  install.packages('rvest')
  library(rvest)
} 

if (!require(stringr)) {
  install.packages('stringr')
  library(stringr)
} 
library(readr)


geneid <- read_csv('../geneid.txt', col_names = F)

# ncbi --------------------------------------------------------------------


geneinfo <- vector(mode = 'character',length = length(geneid$X1))
genenames <- vector(mode = 'character',length = length(geneid$X1))

for (i in seq_len(length(geneid$X1))) {
  
  weburl <- str_c('https://www.ncbi.nlm.nih.gov/gene/?term=',geneid$X1[i])
  NCBIgene <- read_html(weburl, encoding = 'UTF-8') 
  genesummary <- html_nodes(NCBIgene,'body div div form div div div div div div div div dl dd') %>% 
    html_text()
  geneinfo[i] <- genesummary[10]
  genename <- html_nodes(NCBIgene, 'body div div form h1 span') %>% 
    html_text()
  genenames[i] <- genename
  print(i)
  print(genename)
  print(genesummary[10]) 
  
  
}
off.target.genes<- data.frame(name = genenames, info = geneinfo, stringsAsFactors = F)

# uniprot -----------------------------------------------------------------

gene.disseases <- vector(mode = 'character',length = length(geneid$X1))
gene.functions <- vector(mode = 'character',length = length(geneid$X1))
for (i in seq_len(length(off.target.genes$name))) {
  uniproturl <- str_c('https://www.uniprot.org/uniprot/', off.target.genes$name[i],'_HUMAN')
  
  fit <- try(uniprot.gene <- read_html(uniproturl,encoding ='UTF-8'), silent = T)
  if("try-error" %in% class(fit))
  {
    next
  }
  else
  {
    uniprot.gene = fit
  }
  gene.dissease <- html_nodes(uniprot.gene,'.disseaseDescription') %>% 
    html_text()
  print(gene.dissease)
  if (length(gene.dissease) == 0 || nchar(gene.dissease) == 0) {
    gene.dissease <- NA
  }else{
    gene.dissease
  }
  gene.disseases[i] <- str_c(gene.dissease, sep = '', collapse = '')
  gene.function <- html_nodes(uniprot.gene,'div.annotation:nth-child(2)') %>% 
    html_text()
  print(gene.function)
  if (length(gene.function) == 0 || nchar(gene.function) == 0) {
    gene.function <- NA
  }else{
    gene.function
  }
  gene.functions[i] <- gene.function
  print(i)
  print(gene.dissease)
  print(gene.function)
  
}
off.target.genes$func <- gene.functions
off.target.genes$disease <- gene.disseases

oldpath <- getwd()
setwd(dir = '../')
write.csv(off.target.genes, file = 'offtargetgenes.txt',)
setwd(dir = oldpath)

