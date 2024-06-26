############ Volcano plot advance #############

suppressPackageStartupMessages({
library(argparse)
library(tidyverse)
library(Seurat)
library(scRNAtoolVis)
library(cols4all)
library(dplyr)
})

optionParse <- ArgumentParser()
optionParse$add_argument("--input",help='Seurat obj in rds or Diffgene inclued each cell type diff gene in xls', required=T)
optionParse$add_argument("--marker_list",help='marker gene list, [default %(default)s]', default='')
optionParse$add_argument("--output",help='project output dir, [default %(default)s]', default='.')
optionParse$add_argument("--prefix",help='prefix of output plot, [default %(default)s]', default='')
optionParse$add_argument("--minpct",help='min.pct used for FindAllMarkers, [default %(default)s]', default=0.1)
optionParse$add_argument("--logfc",help='logfc threshold used for FindAllMarkers, [default %(default)s]', default=0.25)
optionParse$add_argument("--topgene",help='Numbers of top gene show on the plot, [default %(default)s]', default=3)
args <- optionParse$parse_args()

# check outdir exist
if (!dir.exists(args$output)) dir.create(args$output)

# check outdir, remove the last '/' if have
if (substr(args$output,nchar(args$output),nchar(args$output))=='/'){
      args$output <- substr(args$output,0,nchar(args$output)-1)
}

# check input format
if(substr(args$input,nchar(args$input)-2,nchar(args$input))=='rds'){
      rds <- readRDS(args$input)
      Idents(rds) <- rds$cluster.type
      data <- FindAllMarkers(rds, min.pct = args$minpct, logfc.threshold = args$logfc)
} else {
      xls <- list()
      merge_xls <- lapply(list.files(args$input,pattern = 'xls'), function(file){
            type <- gsub('^.*cluster','',file) %>% gsub('.MK.xls','',.)
            xls[[type]] <- read.table(paste0(args$input,'/',file),header = T) %>% 
            mutate(cluster = type)
      })
      data <- do.call(rbind,merge_xls)
      data$cluster <- as.factor(data$cluster)
}

run_marker <- function(data, marker, topgene, output, prefix){
      colours = c4a('dynamic',length(unique(data$cluster)))
      plot_1 <- jjVolcano(diffData = data,
            aesCol = c('navy','red'),
            tile.col = colours,
            pSize  = 0.5,
            fontface = 'italic',
            topGeneN = topgene,
            base_size = 10,
            flip = T,
            legend.position = c(0.9,0.1),
            polar = T)
      ggsave(paste0(output,'/',prefix,"_Advanced_volcano_circle.pdf"),plot_1)
    
      plot_2 <- jjVolcano(diffData = data,
            aesCol = c('navy','red'),
            tile.col = colours,
            pSize  = 0.5,
            fontface = 'italic',
            topGeneN = topgene,
            base_size = 10,
            flip = T,
            legend.position = c(0.9,0.1))
      ggsave(paste0(output,'/',prefix,"_Advanced_volcano_vertical.pdf"),plot_2)

      plot_3 <- jjVolcano(diffData = data,
            aesCol = c('navy','red'),
            tile.col = colours,
            pSize  = 0.5,
            fontface = 'italic',
            topGeneN = topgene,
            base_size = 10,
            legend.position = c(0.9,0.1))
      ggsave(paste0(output,'/',prefix,"_Advanced_volcano_horizontal.pdf"),plot_3)
}

run_marker(data,args$marker, args$topgene ,args$output, args$prefix)
















