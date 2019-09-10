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
for (i in seq_len(length(geneid$X1))) {
  
  weburl <- str_c('https://www.ncbi.nlm.nih.gov/gene/?term=',geneid$X1[i])
  NCBIgene <- read_html(weburl, encoding = 'UTF-8') 
  genesummary <- html_nodes(NCBIgene,'body div div form div div div div div div div div dl dd') %>% 
    html_text()
  genename <- html_nodes(NCBIgene, 'body div div form h1 span') %>% 
    html_text()
  print(i)
  print(genename)
  print(genesummary[10]) 
  
}








