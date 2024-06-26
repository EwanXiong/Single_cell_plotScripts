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
parser$add_argument("--marker", help = "marker gene file", required=T)

parser$add_argument("--preanno", action='store_true', help = "plot marker before annotaion")
parser$add_argument("--featureplot", action='store_true', help = "feature plot")
parser$add_argument("--vlnplot", action='store_true', help = "vlnplot")
parser$add_argument("--dotplot", action='store_true', help = "dotplot")
parser$add_argument("--heatmap", action='store_true', help = "heatmap")
parser$add_argument("--split_celltype", action='store_true', help = "output featureplot for each celltype dir")

parser$add_argument("--color", help="one col, color file, [default %(default)s]", default='F')
args <- parser$parse_args()

suppressPackageStartupMessages({
library(argparse)
library(Seurat)
library(ggplot2)
library(purrr)
library(patchwork)
library(plyr)
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

getMarker <- function(i, seurat_object, marker){
    marker_name <- marker[marker[,1]==i, 2]
    if(length(marker_name)==0)next
    marker_name <- unlist(strsplit(marker_name, split=","))
    error <- marker_name[!(marker_name %in% rownames(seurat_object))]
    if(length(error)>0){
        print("error marker: ")
        print(paste(error, collapse=", "))
    }
    marker_name <- marker_name[marker_name %in% rownames(seurat_object)] %>% unique()
    return(marker_name)
}

parserMarker <- function(seurat_object, marker, preanno=FALSE, split=FALSE){
    celltype <- levels(seurat_object)
    if(preanno){
        celltype <- marker[,1]
    }
    features <- c()
    features_list <- list()

    for(i in celltype){
        if(!(i %in% marker[,1])){
            print(paste0("not find celltype: ", i))
            next
        }
        marker_name <- getMarker(i, seurat_object, marker)
        features <- c(features, marker_name)
        features_list[[as.character(i)]] <- marker_name
    }

    return(list(features_list, unique(features)))
}

modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75,0,-0.75,0),"cm"), coord_flip = FALSE, legend = FALSE,...){
    p <- VlnPlot(obj, features = feature, pt.size = pt.size,...) +
        labs(x=NULL, y=feature, title=NULL)+
        theme(
            legend.position="none",
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_text(size = rel(1),angle = 0,vjust = 0.5),
            plot.margin = plot.margin)
    if(coord_flip){
        p <- p + theme(
            axis.title.x = element_text(size = rel(0.8), angle = 0),
            plot.margin = plot.margin)+
            scale_y_continuous(position = "right")+
            coord_flip()
    }
    return(p)
}

StackVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(0,0,0,0),"cm"), coord_flip = FALSE,...){
    plot_list <- purrr::map(features,function(x) modify_vlnplot(obj, feature = x, coord_flip = coord_flip,...))
    if(!coord_flip){
        plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),axis.ticks.x = element_line())
        p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
    }else{
        plot_list[[1]] <- plot_list[[1]] + theme(axis.text.y = element_text(angle=0), axis.ticks.y = element_line())
        p <- patchwork::wrap_plots(plotlist = plot_list,ncol = length(features))
    }

    h <- max(2, length(features)/2+max(nchar(levels(obj)))/10)
    w <- max(2, length(levels(obj))/12*8)

    pdf <- paste0(args$outdir, "/", args$prefix, "_stackviolin.pdf")
    ggsave(p, filename = gsub("pdf$", "png", pdf), width=w, height=h, bg="white", dpi=400, units="in")
    ggsave(p, filename = pdf, width=w, height=h, bg="white")
}

dotplotTop <- function(seurat_object, features, exp=FALSE){
    if(length(features)>0){
        # color <- c('#000080', '#FAEBD7', '#B22222')
        color <-  c("#0099CC", "white", "#CC0033")
        len <- length(features)
        dot_size <- ifelse(len>30, 2.5, ifelse(len>25, 3, ifelse(len>20, 3.5, ifelse(len>15, 4, ifelse(len>10, 4.5, 5)))))
        p <- DotPlot(seurat_object, features=features, dot.scale=6) + RotatedAxis() + 
            scale_color_gradientn(colors=color) + 
            theme(panel.border=element_rect(color='black', fill=NA, size=0.5, linetype="solid"), 
                axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), 
                axis.line=element_blank()) + 
            labs(x=NULL, y=NULL)
        pdfwidth <- ifelse(len < 15, 6, len*6/15)
        pdfheight <- ifelse(length(levels(seurat_object)) < 10, 6, length(levels(seurat_object))*3/6)
        pdf <- paste0(args$outdir, "/", args$prefix, "_dotplot.pdf")
        if(len>=35){
            p <- p+theme(legend.position='top', legend.justification=c(0.5, 1))
        }
        if(exp){
            xls <- gsub('pdf$', 'xls', pdf)
            df <- p$data[, -5]
            write.table(df, xls, row.names=F, col.names=T, quote=F, sep='\t')
        }
        ggsave(p, filename = gsub("pdf$", "png", pdf), width=pdfwidth, height=pdfheight, dpi=400, units="in", limitsize=F, bg="white")
        ggsave(p, filename = pdf, width=pdfwidth, height=pdfheight, limitsize=F, bg="white")
    }
}

heatmapTop <- function(seurat_object, features){
    if(length(features)>0){
        # fill_color <- c('#000080', '#FAEBD7', '#B22222')
        fill_color <- c("#0099CC", "white", "#CC0033")
        seurat_object <- ScaleData(seurat_object, features=rownames(seurat_object))
        seurat_object <- subset(seurat_object, downsample=500)
        p <- DoHeatmap(seurat_object, features=features, label=F, group.colors=color, raster=T)+
            scale_fill_gradientn(colors=fill_color)+
            # NoLegend()+
            #scale_fill_viridis()+
            theme(axis.text.y = element_text(size = 8))
        if(length(features) <= 60){
            pdfheight <- 7.5
            pngheight <- 550
        }else{
            pdfheight <- length(features)/10+2
            pngheight <- length(features)*10
        }
        pdf <- paste0(args$outdir, "/", args$prefix, "_heatmap.pdf")
        ggsave(p, filename = gsub("pdf$", "png", pdf), width=6, height=pdfheight, bg="white", dpi=400, units="in")
        ggsave(p, filename = pdf, width=6, height=pdfheight, bg="white")
    }
}

plotFunc <- function(pic, pdf, params=list('pdfwidth'=6, 'pdfheight'=6, 'pngwidth'=480, 'pngheight'=480), png=TRUE){
    pdf(pdf, width=params$pdfwidth, height=params$pdfheight, bg='white')
    print(pic)
    dev.off()
    if(png){
        png <- gsub('pdf$', 'png', pdf)
        png(png, width=params$pngwidth, height=params$pngheight)
        print(pic)
        dev.off()
    }
}

featureplot <- function(seurat_object, features, dir){
    color <- c('lightgrey', '#FC4E07')
    p <- FeaturePlot(seurat_object, cols=color, features=features, reduction="umap") + theme(aspect.ratio=1)
    pdf <- paste0(dir, "/", args$prefix, '_', features,"_featureplot.pdf")
    plotFunc(p, pdf)
}

main <- function(){
    createDir(args$outdir)
    seurat_object <- editRDS(args)

    marker <- read.table(args$marker, head=T, sep = "\t", stringsAsFactors=F)
    marker_list <- parserMarker(seurat_object, marker, args$preanno)
    feature <- marker_list[[2]]
    feature_list <- marker_list[[1]]

    if(args$color!='F'){
        color <- read_lines(args$color, skip_empty_rows=T)
    }
    if(args$vlnplot){
        StackVlnPlot(seurat_object, feature = feature, pt.size = 0, col = color, coord_flip = FALSE)
    }
    if(args$dotplot){
        dotplotTop(seurat_object, feature)
    }
    if(args$heatmap){
        heatmapTop(seurat_object, feature)
    }
    if(args$featureplot & args$split_celltype){
        for(c in names(feature_list)){
            dir <- file.path(args$outdir, c)
            createDir(dir)
            f <- feature_list[[c]]
            llply(f, featureplot, seurat_object=seurat_object, dir=dir)
        }
    }else if(args$featureplot){
        llply(feature, featureplot, seurat_object=seurat_object, dir=args$outdir)
    }
}
main()