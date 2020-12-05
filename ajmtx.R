library(tidyverse)

# your file folder, with files named after 'network*x*'
# files include 
#('network*x*.txt' 'network*x*Labels.txt' 'network*x*summary.tsv' 'pi-title.tsv')
loc <- "/home/zaczou/projects/network_science/data/" 
netname <- "network407x407"

#-------------------construct adjacent matrix-------------------#
mtx <- read_delim(paste(loc,netname,".txt",sep = ""),delim = " ",col_names = F)
names <- read_delim(paste(loc,netname,"Labels.txt",sep = ""),delim = "\n",col_names = F)

trct <- function(string){
  tr <- str_trim(unlist(str_split(string,"\t"))[1])
  print(tr)
  return(tr)}

names <- map(names$X1,trct)
colnames(mtx) <- names
write_tsv(mtx,paste(loc,"_readable.tsv",sep = ""))

#------------------annotate the nodes with title---------------#
info.title <- read_tsv(paste(loc,"summary.tsv",sep = ""))
relation <- info.title[,c("Author","Title","Page")]
info.pi <- read_tsv(paste(loc,"pi-title.tsv",sep = ""))

# quality control
info.pi <- distinct(info.pi[info.pi$title %in% relation$Title,])
relation <- distinct(relation[relation$Title %in% info.pi$title,])

# combine the annotation
authors <- c()
for(title in relation$Title)
{
  tmp = unlist(info.pi[info.pi$title==title,"author"][1,1])
  names(tmp) <- NULL
  print(tmp)
  authors <- c(authors,tmp)
}
result <- mutate(relation,PI = authors)
colnames(result) <- c("Info","Title","Year","PI")

# write to file
write_tsv(result,paste(loc,"_title.tsv",sep = ""))
