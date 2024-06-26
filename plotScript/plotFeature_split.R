suppressMessages({
library(argparse)
})

parser <- ArgumentParser()
parser$add_argument('--seurat_object', help='RDS file')
parser$add_argument('--outdir', help='outdir of result, [default %(default)s]', default='.')
parser$add_argument('--prefix', help='prefix of outdir, [default %(default)s]', default='project')
parser$add_argument("--subsample", help="sample used to analysis, split by , [default %(default)s]", default="all")
parser$add_argument("--subcluster", help="celltype used to analysis, split by , [default %(default)s]", default="all")
parser$add_argument('--splitby', help='sample or column name of meta.data, [default %(default)s]', default='sample')
parser$add_argument("--onePic", help="plot all picture together, [default %(default)s]", default=FALSE, action='store_true')
parser$add_argument("--color", help="one col, color file, [default %(default)s]", default='F')
parser$add_argument("--marker", help = "marker gene file, split by", required=T)
parser$add_argument("--label", action='store_true', help = "Whether to label the clusters")
args <- parser$parse_args()

suppressMessages({
library(tidyverse)
library(Seurat)
library(Hmisc)
})

createDir <- function(dir){
    if(!dir.exists(dir)){
        dir.create(dir, recursive=T)
    }
}

color <- c("#00CC99", "#F0E685", "#99CC00", "#C75127", "#5050FF", "#CDDEB7","#EE82EE",
"#924822","#FFC20A","#33CC00","#0099CC","#7A65A5","#AE1F63","#D595A7",
"#A9A9A9","#E4AF69","#990080","#14FFB1","#003399", "#FF1463","#D60047",
"#CE3D32","#996600","#809900","#008099","#0A47FF","#660099","#FFD147",
"#339900","#5A655E","#4775FF","#802268","#8B0000","#BBFFFF","#00868B",
"#00008B","#CDBA96","#FF69B4","#3B1B53","#D2B48C","#00FFFF","#836FFF",
"#BA6338","#991A00","#D58F5C","#E7C76F","#CC9900","#749B58","#6BD76B",
"#009966","#00D68F","#5DB1DD","#466983","#837B8D","#612A79","#990033"
)
################# perpare marker #####################



getMarker <- function(seurat_object, marker){
    marker_name <- unlist(strsplit(marker, split=','))
    if(length(marker_name)>0)next
    #marker_name <- unlist(strsplit(marker_name, split=","))
    error <- marker_name[!(marker_name %in% rownames(seurat_object))]
    if(length(error)>0){
        print("error marker: ")
        print(paste(error, collapse=", "))
    }
    marker_name <- marker_name[marker_name %in% rownames(seurat_object)] %>% unique()
    return(marker_name)
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

featureParams <- function(scRNA){
    cluster <- levels(scRNA)
    cluster_num <- length(cluster)
    if(NA %in% (as.numeric(cluster))){
        if(cluster_num>20){
            ncol <- ceiling(cluster_num/20)
            ratio <- 1.8
        }else{
            ncol <- 1
            ratio <- 1.6
        }
    }else{
        if(cluster_num>20){
            ncol <- ceiling(cluster_num/20)
            ratio <- 12/9
        }else{
            ncol <- 1
            ratio <- 12/10
        }
    }
    return(list(ncol=ncol, ratio=ratio))
}

featureSplitCluster <- function(scRNA, marker, metadata, dir, prefix, color,label =F){
    if(length(unique(scRNA@meta.data[, metadata])) > 1){
        params <- featureParams(scRNA)
        cluster <- levels(scRNA)
        clusterColor <- color[1:length(cluster)]
        names(clusterColor) <- cluster
        meta <- unique(scRNA@meta.data[, metadata])
        df <- data.frame(barcode=colnames(scRNA), splitby=scRNA@meta.data[[metadata]])
        if (label){
        show <- TRUE
        } else {
        show <- FALSE
        }
        for(m in meta){
            df_sub <- df %>% filter(splitby==m)
            scRNA_sub <- subset(scRNA, cells=df_sub$barcode)
            color_split <- clusterColor[levels(scRNA_sub)]
            for(gene in marker){
                dir_sub <- paste0(dir, '/split', capitalize(metadata))
                if(!dir.exists(dir_sub)) dir.create(dir_sub)
                for(method in c('umap','tsne')){
                    x_max <- max(Embeddings(scRNA, reduction=method)[,1])
                    x_min <- min(Embeddings(scRNA, reduction=method)[,1])
                    y_max <- max(Embeddings(scRNA, reduction=method)[,2])
                    y_min <- min(Embeddings(scRNA, reduction=method)[,2])
                    methodUpper <- str_to_upper(method)
                    p <- FeaturePlot(scRNA_sub, features = gene, label = show, reduction = method, cols = c('#CECECE','#FF5E00')) + 
                    theme(aspect.ratio=1) + 
                    guides(color = guide_legend(override.aes = list(size=4), ncol=params$ncol)) + 
                    lims(x=c(x_min, x_max), y=c(y_min, y_max)) + theme(legend.position = "right")  
                    ggsave(p, path = dir_sub, filename = paste0(prefix, '.', capitalize(metadata), '_', m, '_', methodUpper, '_FeaturePlot_', gene, '.png'), device = 'png', width=7*params$ratio)
                    ggsave(p, path = dir_sub, filename = paste0(prefix, '.', capitalize(metadata), '_', m, '_', methodUpper, '_FeaturePlot_', gene, '.pdf'), device = 'pdf', width=7*params$ratio)
                }
            }
        }
    }
}

featureSplitClusterOnePic <- function(scRNA, marker, metadata, dir, prefix, color,label=F){
    params <- featureParams(scRNA)
    if (label){
        show <- TRUE
    } else {
        show <- FALSE
    }
    for(gene in marker){
        for (method in c('umap','tsne')){
            methodUpper <- str_to_upper(method)
            p <- FeaturePlot(scRNA, features = gene, label = show, reduction = method, split.by = metadata, cols = c('#CECECE','#FF5E00')) +theme(legend.position = "right")
            num <- length(unique(scRNA@meta.data[, metadata]))
            ggsave(p, path = dir, filename = paste0(prefix, '.', capitalize(metadata),'_', methodUpper , '_FeaturePlot_', gene, '.png'), device = 'png', width=7*max(params$ratio, num/1), limitsize = FALSE)
            ggsave(p, path = dir, filename = paste0(prefix, '.', capitalize(metadata),'_', methodUpper , '_FeaturePlot_', gene, '.pdf'), device = 'pdf', width=7*max(params$ratio, num/1), limitsize = FALSE)
        }
    }
}

main <- function(){
    createDir(args$outdir)
    seurat_object <- editRDS(args)
    marker_name <- unlist(strsplit(args$marker, split=','))
    if(args$color!='F'){
        color <- read_lines(args$color, skip_empty_rows=T)
    }
    if(args$onePic){
        featureSplitClusterOnePic(seurat_object, marker_name, args$splitby, args$outdir, args$prefix, color, args$label)
    }else{
        featureSplitCluster(seurat_object, marker_name, args$splitby, args$outdir, args$prefix, color, args$label)
    }
}
main()
