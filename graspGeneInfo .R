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
library(curl)

geneid <- read_csv('../geneid.txt', col_names = F)
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

oldpath <- getwd()
setwd(dir = '../')
write.csv(off.target.genes, file = 'offtargetgenes.txt',)
setwd(dir = oldpath)
i <- 3
genecardurl <- str_c('https://www.genecards.org/cgi-bin/carddisp.pl?gene=', off.target.genes$name[i])
genecard.gene <- read_html(curl(genecardurl, 
                                handle = curl::new_handle("useragent" = "Mozilla/5.0")),
                           encoding ='UTF-8')



