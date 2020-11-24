library(tidyverse)

# your file folder, with files named after 'network*x*'
# files include ('network*x*.txt' 'network*x*Labels.txt' 'network*x*summary.tsv')
loc <- "/home/zaczou/projects/network_science/network82x82" 

#-------------------construct adjacent matrix-------------------#
mtx <- read_delim(paste(loc,".txt",sep = ""),delim = " ",col_names = F)
names <- read_delim(paste(loc,"Labels.txt",sep = ""),delim = "\n",col_names = F)

trct <- function(string){
  tr <- str_trim(unlist(str_split(string,"\t"))[1])
  print(tr)
  return(tr)}

names <- map(names$X1,trct)
colnames(mtx) <- names
write_tsv(mtx,paste(loc,"_readable.tsv",sep = ""))

#------------------annotate the nodes with title---------------#
info <- read_tsv(paste(loc,"summary.tsv",sep = ""))
relation <- info[,c("Author","Title")]
write_tsv(relation,paste(loc,"_title.tsv",sep = ""))
