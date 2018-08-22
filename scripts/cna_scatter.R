setwd("~/Desktop/chrom_cna/Melanoma")
library(ggplot2)

generate_plot <- function(fname){
  df <- read.csv(fname, sep="\t")
  name <- basename(fname)
  name <- gsub(".xls","",name)
  p <- ggplot(data=df)+
    geom_point(aes(x = gene_middle, y = del_total), color = 'deepskyblue') +
    geom_point(aes(x = gene_middle, y = deep_del), color = 'blue3') +
    geom_point(aes(x = gene_middle, y = amp_total), color = 'darksalmon') +
    geom_point(aes(x = gene_middle, y = deep_amp), color = 'brown4') 
  q <-p + labs(x=name, y="Percentage", title="Melanoma")
  png(filename=paste(name, "_melanoma.png"))
  plot(q)
  dev.off()
}

generate_deep_plot <- function(fname){
  df <- read.csv(fname, sep="\t")
  name <- basename(fname)
  name <- gsub(".xls","",name)
  p <- ggplot(data=df)+
    geom_point(aes(x = gene_middle, y = deep_del), color = 'blue3') +
    geom_point(aes(x = gene_middle, y = deep_amp), color = 'brown4') 
  q <-p + labs(x=name, y="Percentage", title="Melanoma")
  png(filename=paste(name, "_melanoma_deep.png"))
  plot(q)
  dev.off()  
}

myFiles <- list.files(pattern="*.xls")
for (fname in myFiles){
  #generate_plot(fname)
  generate_deep_plot(fname)
}
