
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(reshape2)

#--------------------------------------preprocessing-------------------------------------#

loc <- "/home/zaczou/projects/network_science/data/" 
netname <- "network927x927"

# construct adjacent matrix #
mtx <- read_delim(paste(loc,netname,".txt",sep = ""),delim = " ",col_names = F)
names <- read_delim(paste(loc,netname,"Labels.txt",sep = ""),delim = "\n",col_names = F)

trct <- function(string){
  tr <- str_trim(unlist(str_split(string,"\t"))[1])
  return(tr)}

names <- map(names$X1,trct)
colnames(mtx) <- names
write_tsv(mtx,paste(loc,netname,"_readable.tsv",sep = ""))

# annotate the nodes with title #
info.title <- read_tsv(paste(loc,netname,"summary.tsv",sep = ""))
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

# annotate the bulk lab
meta.lab <- read_tsv(paste(loc,"PI.tsv",sep = ""), col_names = F)
lab_list = c()
for(pi in result$PI)
{
  if(is.na(pi)) lab_list <- c(lab_list, "NA")
  else lab_list <- c(lab_list, unlist(meta.lab$X1[meta.lab$X4==pi]))
}
lab_list <- unlist(lab_list)
result <- result %>% mutate(Lab = lab_list)

# write to file
write_tsv(result,paste(loc,netname,"_title.tsv",sep = ""))


#--------------------------------------analytics------------------------------------------#

cres <- read_tsv("/home/zaczou/projects/network_science/data/cluster_result.tsv", col_names = F)
meta <- result
names <- unlist(names)

# filter out zero clusters
for(cols in colnames(cres)) if(sum(cres[,cols])<1) cres[,cols] = NULL
hc <- max.col(cres, ties.method = "first")
hc.anno <- tibble(info = names, cluster = hc)
write_tsv(hc.anno, "/home/zaczou/projects/network_science/data/hclust.tsv")

# filter out small clusters
cres.f <- cres
for(cols in colnames(cres)) if(sum(cres.f[,cols])<10) cres.f[,cols] = NULL

# pca analysis
pca <- function(data, Year)
{
  pr.out <- prcomp(data, scale. = T) #pca analysis
  pr.var <- pr.out$sdev^2      #calculate variance percentage
  prr <- as.data.frame(pr.out$rotation)
  pve <- 100*pr.var/sum(pr.var)
  res <- as_tibble(pr.out$x)   #ready to plot
  res %>%    #plot the result, PC1 and PC2
    ggplot(aes(PC1,PC2)) +
    geom_point(aes(color=Year), size = 2.5, alpha = 0.5) +
    geom_segment(aes(x=0,y=0,xend=2*prr[1,1], yend=2*prr[1,2]), 
                 arrow = arrow(length = unit(0.1, "cm")), show.legend = T) +
    geom_segment(aes(x=0,y=0,xend=2*prr[2,1], yend=2*prr[2,2]), 
                 arrow = arrow(length = unit(0.1, "cm")), show.legend = T) +
    geom_segment(aes(x=0,y=0,xend=2*prr[3,1], yend=2*prr[3,2]), 
                 arrow = arrow(length = unit(0.1, "cm")), show.legend = T) +
    labs(x = paste("PC1: ", round(pve[1]), "% variance", sep = ""), 
         y = paste("PC2: ", round(pve[2]), "% variance", sep = ""))
}

lab <- c()
for(tmp in names)
{
  lab <- c(lab, unlist(meta$Lab[meta$Info==tmp]))
}
pca(cres.f, lab)

# data filtering
meta <- dplyr::inner_join(hc.anno, meta, by = c("info"="Info"))
meta <- meta %>%
  dplyr::filter(Year>=2004)

# PI to cluster heatmaps
meta.pi <- meta %>% 
  group_by(PI) %>%
  summarize(cluster) %>%
  count(cluster)
df <- data.frame(row.names = 1:15)
for(ii in 1:nrow(meta.pi))
{
  if(!is.na(meta.pi[ii,"PI"]))
    df[unlist(meta.pi[ii,"cluster"]),unlist(meta.pi[ii,"PI"])] = unlist(meta.pi[ii,"n"])
}
df[is.na(df)] <- 0
for(pi in colnames(df))
{
  df[,pi] <- df[,pi]/sum(df[,pi])
}
pheatmap(df, colorRampPalette(rev(brewer.pal(n = 4, name ="RdYlBu")))(300),
         labels_col = colnames(df), labels_row = 1:15, cluster_rows = F)

# year to cluster heatmaps
meta.year <- meta %>% 
  group_by(Year) %>%
  summarize(cluster) %>%
  count(cluster) %>%
  arrange(cluster,year)
df <- data.frame(row.names = 2004:2019)
for(ii in 1:nrow(meta.year))
{
  df[as.character(unlist(meta.year[ii,"Year"])),as.character(unlist(meta.year[ii,"cluster"]))] = unlist(meta.year[ii,"n"])
}
df[is.na(df)] <- 0
df <- t(df)
for(year in colnames(df))
{
  df[,year] <- df[,year]/sum(df[,year])
}
df <- df[,-9]
df <- df[,-16]
pheatmap(df, colorRampPalette(rev(brewer.pal(n = 4, name ="RdYlBu")))(300),
         labels_col = colnames(df), labels_row = 1:15, cluster_cols = F, cluster_rows = F)

ent <- function(vec)
{
  vec <- vec[vec!=0]
  return(sum(-vec*log(vec)))
}
entp <- tibble(year = colnames(df), entropy = unname(apply(df,2,ent)))
entp %>% ggplot(aes(year,entropy,group=1)) +
  geom_point() +
  geom_smooth(formula = y~log(x))

# pi-trajectory
one_pi <- function(name)
{
  cres.tmp <- cres.f[names %in% unlist(meta[meta$PI==name,"info"]),]
  cres.tmp <- cres.tmp[,apply(cres.tmp,2,sum)>0]
  pca(cres.tmp, meta$Year[meta$info %in% names[names %in% unlist(meta[meta$PI==name,"info"])]])
}

# degree distribution
deg <- data.frame("0.005" = rowSums(mtx > 0.005),
                  "0.01" = rowSums(mtx > 0.01),
                  "0.05" = rowSums(mtx > 0.05),
                  "0.1" = rowSums(mtx > 0.1))
deg.m <- melt(deg)
colnames(deg.m) <- c("threshold", "value")
deg.m$threshold <- substring(deg.m$threshold, 2)
ggplot(deg.m, aes(x=value, fill=threshold)) + 
  geom_density(alpha=0.25) +
  xlab("k") + ylab("P(k)")
ggplot(deg.m, aes(x=value, fill=threshold)) +
  geom_density(alpha=0.25) +
  scale_x_log10() +
  scale_y_log10() +
  xlab("k") + ylab("P(k)")
print(paste(mean(deg),var(hdeg$counts)))
plot(log2(hdeg$breaks[-1]),log2(hdeg$counts))
