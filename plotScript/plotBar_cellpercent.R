suppressPackageStartupMessages({
library(argparse)
library(tidyverse)
})

parser <- function() {
    parser <- ArgumentParser()
    parser$add_argument("--cellcountfile", help = "the cell count file, come from celllabel folder", required = T)
    parser$add_argument("--outdir", help = "the outdir of result", default = ".")
    parser$add_argument("--prefix", help = "the prefix of result", default = "Integration")
    parser$add_argument("--celltypeorder", help = "sort the celltype, split by ','")
    parser$add_argument("--sampleorder", help = "sort the sample, split by ','")
    args <- parser$parse_args()
    return(args)
}

create_dir <- function(dir){
    if(!dir.exists(dir)){
        dir.create(dir, recursive=T)
    }
}

get_celltype_count <- function(cellcountfile, celltypeorder, sampleorder) {
  cellcount <- read.table(cellcountfile, head = T, sep = "\t")
  cellcount <- cellcount %>% mutate(Cellname = cellcount[, 1]) %>% select(Cellname, everything()) %>% select(-2) %>%
    pivot_longer(cols = -1, names_to = "Sample", values_to = "Count") %>% group_by(Cellname) %>% mutate(Percent = Count/sum(Count))
  if(!is.null(celltypeorder)){
    cellcount$Cellname <- factor(cellcount$Cellname, levels = celltypeorder)
  }
  if(!is.null(sampleorder)){
    cellcount$Sample <- factor(cellcount$Sample, levels = sampleorder)
  }
  return(cellcount)
}

plot_stack_barplot <- function(cellcount) {
    p <- ggplot(cellcount, aes(x = Cellname, y = Percent, fill = Sample))+
        geom_bar(stat = "identity", alpha = 0.9) +
        geom_text(aes(label=Count),position = position_stack(vjust =0.5),size=3) +
        scale_fill_brewer(palette = "Set3") +
        scale_y_continuous(expand = c(0,0),label=scales::percent) +
        labs(y = "Percent") + theme_bw() + theme(axis.title.x =element_blank(), legend.key.size = unit(0.5,'cm'))
    return(p)
}

main <- function() {
    args <- parser()
    create_dir(args$outdir)
    celltypeorder <- NULL
    sampleorder <- NULL
    if(!is.null(args$celltypeorder)) {celltypeorder <- unlist(strsplit(args$celltypeorder, ","))}
    if(!is.null(args$sampleorder)) {sampleorder <- unlist(strsplit(args$sampleorder, ","))}
    cellcount <- get_celltype_count(args$cellcountfile, celltypeorder, sampleorder)
    p <- plot_stack_barplot(cellcount)
    ggsave(paste0(args$outdir, "/", args$prefix, ".celltype.CellsNumberPerSample.png"), p, width = length(cellcount$Cellname)/2, height= 8)
    ggsave(paste0(args$outdir, "/", args$prefix, ".celltype.CellsNumberPerSample.pdf"), p, width = length(cellcount$Cellname)/2, height= 8)
}

main()