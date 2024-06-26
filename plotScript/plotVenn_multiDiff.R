#' @name vennplot, up genes
#' @date 20230224
#' @param name, the name of celltype or sample, the number of name are less than 4 and more than 2, split by ","
#' @param degfile, the deg genes of Seurat::FindMarker, the number of degfiles are less than 4 and more than 2
#' @param outdir, the output of reslut
#' @param prefix, the prefix of result
#' @param top, the top genes need to compare


suppressPackageStartupMessages({
library(argparse)
library(VennDiagram)
library(RColorBrewer)
library(tidyverse)
})

Parser <- function() {
    parser <- ArgumentParser()
    parser$add_argument("--namefile", help = "the name of celltype or sample, split by ,", required = T)
    parser$add_argument("--degfile1", help = "the deg file contain deg genes", required = T)
    parser$add_argument("--degfile2", help = "the deg file contain deg genes", required = T)
    parser$add_argument("--degfile3", help = "the deg file contain deg genes")
    parser$add_argument("--degfile4", help = "the deg file contain deg genes")
    parser$add_argument("--outdir", help = "the output of result", default = ".")
    parser$add_argument("--prefix", help = "the prefix of result", default = "Integration")
    parser$add_argument("--top", help = "the number of top genes", default = "all")
    args <- parser$parse_args()
    return(args)
}

GetGenes <- function(degfile, args) {
    diffgenes <- read.table(degfile, header = T, sep = "\t", quote = "")
    diffgenes$genename <- diffgenes[, 1]
    diffgenes <- diffgenes %>% select(genename, everything()) %>% subset(avg_log2FC > 0)
    if (args$top != "all") {
        diffgenes <- diffgenes %>% group_by(genename) %>% top_n(n = as.numeric(args$top), wt = avg_log2FC)
    }
    return(diffgenes)
}

PlotVenn <- function(name, degfile1, degfile2, degfile3, degfile4, args) {
    if (!is.null(degfile3) && !is.null(degfile4)) {
        venn.diagram(x = list(degfile1$genename, degfile2$genename, degfile3$genename, degfile4$genename), fill = c(brewer.pal(7, "Set1")[1:4]), 
        category.names = c(name[[1]], name[[2]], name[[3]], name[[4]]), 
        alpha = c(0.5, 0.5, 0.5, 0.5), cex = 2, cat.cex = 3, cat.fontface = 4, lty = 2, fontfamily =3, resolution = 300, scaled = F, cat.default.pos = "outer", 
        filename = paste0(args$outdir, "/", args$prefix, ".", args$top, ".vennplot.tiff"))
    } else if (!is.null(degfile3) && is.null(degfile4)) {
        venn.diagram(x = list(degfile1$genename, degfile2$genename, degfile3$genename), fill = c(brewer.pal(7, "Set1")[1:3]), category.names = c(name[[1]], name[[2]], name[[3]]), 
        alpha = c(0.5, 0.5, 0.5), cex = 2, cat.cex = 3, cat.fontface = 4, lty = 2, fontfamily =3, resolution = 300, scaled = F, cat.default.pos = "outer",
        filename = paste0(args$outdir, "/", args$prefix, ".", args$top, ".vennplot.tiff"))
    } else {
        venn.diagram(x = list(degfile1$genename, degfile2$genename), fill = c(brewer.pal(7, "Set1")[1:2]), category.names = c(name[[1]], name[[2]]), 
        alpha = c(0.5, 0.5), cex = 2, cat.cex = 3, cat.fontface = 4, lty = 2, fontfamily =3, resolution = 300, scaled = F, cat.default.pos = "outer",
        filename = paste0(args$outdir, "/", args$prefix, ".", args$top, ".vennplot.tiff"))
    }
}

main <- function() {
    args <- Parser()
    name <- str_trim(unlist(read.table(args$namefile))) %>% str_split(pattern = ",") %>% unlist()
    degfile1 <- GetGenes(args$degfile1, args)
    degfile2 <- GetGenes(args$degfile2, args)
    degfile3 <- NULL
    degfile4 <- NULL
    if(!is.null(args$degfile3)) {degfile3 <- GetGenes(args$degfile3, args)}
    if(!is.null(args$degfile4)) {degfile4 <- GetGenes(args$degfile4, args)}
    PlotVenn(name, degfile1, degfile2, degfile3, degfile4, args)
}

main()
