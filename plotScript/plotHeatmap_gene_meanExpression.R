suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(tidyverse)
    library(argparse)
    library(Seurat)
})

Parser <- function(){
    parser <- ArgumentParser()
    parser$add_argument("--seurat_object", help = "the scRNA RDS file", required = T)
    parser$add_argument("--outdir", help = "the output of result", required = T, default = ".")
    parser$add_argument("--prefix", help = "the prefix of output", default = "Integration")
    parser$add_argument("--assay", help = "RNA or SCT", default = "RNA")
    parser$add_argument("--subsample", help = "sample used to analysis, split by,", default = "all")
    parser$add_argument("--subcluster", help = "celltype used to analysis, split by,", default = "all")
    parser$add_argument("--marker", help = "three col file, first col is genesetname, second col is genename, genename split by ,", required = T)
    parser$add_argument("--color", help = "one col, the color correspond with celltype, split by ,", required = T)
    parser$add_argument("--clusterorder", help = "the cluster order", default = "F")
    parser$add_argument("--ident", help = "col of meatadata")
    parser$add_argument("--rowsplit", help = "split the row", action = "store_true")
    args <- parser$parse_args()
    return(args)
}

SubsetRDS <- function(seuratObj, subset=list()){
    if(length(subset)!=0){
        meta <- seuratObj@meta.data
        for(name in names(subset)){
            if(subset[[name]][1]!='all'){
                meta <- meta[meta[, name] %in% subset[[name]], ]
            }
            print(unique(meta[, name]))
        }
        seuratObj <- subset(seuratObj, cells=rownames(meta))
    }
    return(seuratObj)
}

EditRds <- function(args){
    seurat_object <- readRDS(args$seurat_object)
    if(args$subsample!='all'){
        args$subsample <- unlist(strsplit(args$subsample, split=','))
    }
    if(args$subcluster!='all'){
        args$subcluster <- unlist(strsplit(args$subcluster, split=','))
        seurat_object <- subset(seurat_object, idents=args$subcluster)
    }
    subset_list <- list()
    if(!("sample" %in% colnames(seurat_object@meta.data))){
        seurat_object$sample <- seurat_object$orig.ident
    }
    subset_list$sample <- args$subsample

    seurat_object <- SubsetRDS(seurat_object, subset=subset_list)
    if(!is.null(args$ident)){
        Idents(seurat_object) <- seurat_object@meta.data[, args$ident]
    }
    if(args$clusterorder!='F'){
        args$clusterorder <- unlist(strsplit(args$clusterorder, split=','))
        levels(seurat_object) <- args$clusterorder
    }
    DefaultAssay(seurat_object) <- args$assay
    print(table(Idents(seurat_object)))
    return(seurat_object)
}

GetMarker <- function(i, seurat_object, marker){
    features <- c()
    row_split <- c()
    for (i in seq_along(marker[, 1])){
      marker_name <- marker[i, 2]
      if (length(marker_name) == 0)next
      marker_name = unlist(strsplit(marker_name, split = ","))
      error <- marker_name[!(marker_name %in% rownames(seurat_object))]
      if (length(error) > 0){
        print("error marker: ")
        print(paste(error, collapse = ", "))
      }
      marker_name <- marker_name[marker_name %in% rownames(seurat_object)] %>% unique()
      features <- c(features, marker_name)
      row_split <- c(row_split, rep(marker[i, 1], length(marker_name)))
    }
    row_split <- data.frame(times = row_split)
    row_split$times <- factor(row_split$times, levels = marker[, 1])
    return(list(features, row_split))
}

PlotMeanheatmap <- function(seurat_object, features, color, args, ...){
    p <- AverageHeatmap(seurat_object, markerGene = features, myanCol = color, annoCol = T, border = T, row_title = NULL, ...)
    pdf(paste0(args$outdir,"/",args$prefix,".meanheatmap.pdf"), width = 9.5/11*length(table(seurat_object@active.ident)), height = 10.5/41*length(features))
    print(p)
    dev.off()
    png(paste0(args$outdir,"/",args$prefix,".meanheatmap.png"), width = 760/11*length(table(seurat_object@active.ident)), height = 840/41*length(features))
    print(p)
    dev.off()
}

main <- function(){
    args <- Parser()
    if(!dir.exists(args$outdir)){dir.create(args$outdir, recursive=T)}
    seurat_object <- EditRds(args)
    marker <- read.table(args$marker, sep = "\t", header = T, stringsAsFactors = F)
    color <- str_trim(unlist(read.table(args$color, sep = ",")))
    features <- GetMarker(i, seurat_object, marker)[[1]]
    row_split <- GetMarker(i, seurat_object, marker)[[2]]
    if(args$rowsplit){
        PlotMeanheatmap(seurat_object, features, color, args, row_split = row_split$times, row_gap = unit(5, "mm"))
    } else{
        PlotMeanheatmap(seurat_object, features, color, args)
    }
    
}

main()
