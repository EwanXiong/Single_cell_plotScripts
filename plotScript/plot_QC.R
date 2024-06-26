suppressPackageStartupMessages({
library(argparse)
})

parser <- ArgumentParser()
parser$add_argument('--seurat_object', help='RDS file')
parser$add_argument('--outdir', help='outdir of result, [default %(default)s]', default='.')
parser$add_argument('--prefix', help='prefix of outdir, [default %(default)s]', default='program')
parser$add_argument('--assay', help='RNA or SCT, [default %(default)s]', default='RNA')
parser$add_argument("--subsample", help="sample used to analysis, split by , [default %(default)s]", default="all")
parser$add_argument("--subcluster", help="celltype used to analysis, split by , [default %(default)s]", default="all")
parser$add_argument("--ident", help="colname of metadata, used for celltype, cell.type cluster.type seurat_clusters, [default %(default)s]")
parser$add_argument("--clusterorder", help = "the cluster order , [default %(default)s]", default='F')
parser$add_argument("--vlnplot", action='store_true', help = "vlnplot")
parser$add_argument("--featureplot", action='store_true', help = "featureplot")
parser$add_argument("--features", help = "nCount_RNA, nFeature_RNA, percent.mt, percent.rb, split by ,", default="nCount_RNA")

parser$add_argument("--color", help="one col, color file, [default %(default)s]", default='F')
args <- parser$parse_args()

suppressPackageStartupMessages({
library(argparse)
library(Seurat)
library(purrr)
library(patchwork)
library(plyr)
library(tidyverse)
})

color <- c("#00CC99", "#F0E685", "#99CC00", "#C75127", "#5050FF", "#CDDEB7","#EE82EE",
"#924822","#FFC20A","#33CC00","#0099CC","#7A65A5","#AE1F63","#D595A7",
"#A9A9A9","#E4AF69","#990080","#14FFB1","#003399", "#FF1463","#D60047",
"#CE3D32","#996600","#809900","#008099","#0A47FF","#660099","#FFD147",
"#339900","#5A655E","#4775FF","#802268","#8B0000","#BBFFFF","#00868B",
"#00008B","#CDBA96","#FF69B4","#3B1B53","#D2B48C","#00FFFF","#836FFF",
"#BA6338","#991A00","#D58F5C","#E7C76F","#CC9900","#749B58","#6BD76B",
"#009966","#00D68F","#5DB1DD","#466983","#837B8D","#612A79","#990033"
)

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

# plotViolin <- function(feature, seurat_object, pt.size = 0, col = color, coord_flip = FALSE){
#     meta.data <- seurat_object@meta.data %>% select(sample, feature)
#     p <- ggplot(meta.data, aes(sample, feature)) +
#                       geom_violin(aes(fill = sample)) +
#                       geom_jitter(alpha = 0.5, size = pt.size) +
#                       theme_classic() +
#                       xlab('') +
#                       theme(legend.position = '')
#     pdf <- paste0(args$outdir, "/", args$prefix, "_", feature, "_violinplot.pdf")
#     plotFunc(p, pdf)
# }

plotFunc <- function(pic, pdf, pdfwidth=6, pdfheight=6, png=TRUE){
    pdf(pdf, width=pdfwidth, height=pdfheight, bg='white')
    print(pic)
    dev.off()
    if(png){
        png <- gsub('pdf$', 'png', pdf)
        png(png, width=pdfwidth, height=pdfheight, bg='white', res=300, units="in")
        print(pic)
        dev.off()
    }
}

featureplot <- function(feature, seurat_object, dir){
    color <- c('lightgrey', '#FC4E07')
    p <- FeaturePlot(seurat_object, cols=color, features=feature, reduction="umap") + theme(aspect.ratio=1)
    pdf <- paste0(dir, "/", args$prefix, '_', feature,"_featureplot.pdf")
    plotFunc(p, pdf)
}

plotViolin <- function(feature, seurat_object, dir, ...){
    p <- VlnPlot(seurat_object, features = feature, ...) +
        labs(x=NULL, y=feature, title=NULL)+
        theme(
            legend.position="none",
            axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
        )
    width <- max(6, length(levels(seurat_object))/3)
    pdf <- paste0(dir, "/", args$prefix, '_', feature,"_violinplot.pdf")
    plotFunc(p, pdf, pdfwidth=width, pdfheight=2)
}

main <- function(){
    createDir(args$outdir)
    seurat_object <- editRDS(args)

    if(args$color!='F'){
        color <- read_lines(args$color, skip_empty_rows=T)
    }

    features <- unlist(strsplit(args$features, split=","))
    for(f in features){
        if(!f %in% colnames(seurat_object@meta.data) & f=="percent.rb"){
            seurat_object[["percent.rb"]] <- PercentageFeatureSet(seurat_object, pattern = "^RPL|^RPS|^MRPL|^MRPS")
        }
    }
    if(args$vlnplot){
        llply(features, plotViolin, seurat_object, dir=args$outdir, pt.size = 0, col = color)
    }
    if(args$featureplot){
        llply(features, featureplot, seurat_object, dir=args$outdir)
    }
    
}
main()