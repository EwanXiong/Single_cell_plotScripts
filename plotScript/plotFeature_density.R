suppressMessages({
library(argparse)
})

parser <- ArgumentParser()
parser$add_argument("--seurat_object", help='RDS file')
parser$add_argument("--outdir", help='outdir of result, [default %(default)s]', default='.')
parser$add_argument("--subsample", help="sample used to analysis, split by , [default %(default)s]", default="all")
parser$add_argument("--subcluster", help="celltype used to analysis, split by , [default %(default)s]", default="all")
parser$add_argument("--marker", help = "marker gene file, split by", required=T)
parser$add_argument("--axis", action='store_true', help = "Whether to plot the x, y axis, [default %(default)s]", default=FALSE)
parser$add_argument("--reduction", help="reduction used for plot, tsne or umap , [default %(default)s]", default="umap")
parser$add_argument("--merge", action='store_true', help = "Whether to plot merged marker gene plot, [default %(default)s]", default=FALSE)
args <- parser$parse_args()

suppressMessages({
library(tidyverse)
library(Seurat)
library(Nebulosa)
})

createDir <- function(dir){
    if(!dir.exists(dir)){
        dir.create(dir, recursive=T)
    }
}

subsetRDS <- function(seuratObj, subset=list()){
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

editRDS <- function(args){
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

    seurat_object <- subsetRDS(seurat_object, subset=subset_list)

    print(table(Idents(seurat_object)))
    return(seurat_object)
}

density_plot <- function(seurat_obj, outdir, features, reduction, merge = FALSE, axis = TRUE){
    if (axis != T){
        axis <- theme(axis.line = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),
                      axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())
    } else {
        axis <- NULL
    }
    reductionUpper <- str_to_upper(reduction)
    if (merge != T){
        joint <- NULL
        for (gene in features){
            plots <- plot_density(seurat_obj, gene, reduction = reduction, joint = joint) + axis
            ggsave(paste0(outdir,'/', gene, '_', reductionUpper,'_density_plot.png'), plots, height = 5, width = 5)
            ggsave(paste0(outdir,'/', gene, '_', reductionUpper,'_density_plot.pdf'), plots)
        }
    } else {
        joint <- T
        plot_num <- length(features) + 1
        plots <- plot_density(seurat_obj, features, reduction = reduction, joint = joint) + axis
        merge_plot <- plots[[plot_num]]
        ggsave(paste0(outdir,'/',paste(features,sep = '', collapse = '+'), '_', reductionUpper, '_density_plot.png'), merge_plot, height = 5, width = 5)
        ggsave(paste0(outdir,'/',paste(features,sep = '', collapse = '+'), '_', reductionUpper, '_density_plot.pdf'), merge_plot)
    }
    
}

main <- function(){
    createDir(args$outdir)
    seurat_object <- editRDS(args)
    marker_name <- unlist(strsplit(args$marker, split=','))
    reduction <- unlist(args$reduction)
    if(args$axis){
        axis <- TRUE
    }else{
        axis <- FALSE
    }
    if(args$merge){
        merge <- TRUE
    }else{
        merge <- FALSE
    }
    density_plot(seurat_obj=seurat_object, outdir=args$outdir, features=marker_name, reduction=reduction, merge = merge, axis = axis)
}
main()









