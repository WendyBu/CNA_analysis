# install.packages("pheatmap")
# install.packages("gridExtra")
library('ggplot2')
library('pheatmap')
library('grid')
library('gridExtra')


save_pheatmap <- function(x, filename, width=720, height=480) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename,width = width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}



myFiles <- list.files(pattern="*.xls")
for (fname in myFiles){
  name=basename(fname)
  name = gsub(".xls","", name)
  name = gsub("heatmap_","",name)
  my_data <- read.table(fname, sep="\t", stringsAsFactors = FALSE, header=TRUE)
  df <- as.matrix(my_data[ , c(8:ncol(my_data))])
  rnames <- my_data[ , 1]
  rownames(df) <- rnames
  q<-pheatmap(df, cluster_rows=FALSE,cluster_cols=TRUE, show_colnames =FALSE, show_rownames =FALSE, main=name)
  save_pheatmap(q, filename=name)
}


