suppressMessages({
library(argparse)
})

parser <- ArgumentParser()
parser$add_argument('--seurat_object', help='RDS file')
parser$add_argument('--outdir', help='outdir of result, [default %(default)s]', default='.')
parser$add_argument('--prefix', help='prefix of outdir, [default %(default)s]', default='project')
parser$add_argument("--subsample", help="sample used to analysis, split by , [default %(default)s]", default="all")
parser$add_argument("--subcluster", help="celltype used to analysis, split by , [default %(default)s]", default="all")
parser$add_argument("--color", help="one col, color file, [default %(default)s]", default='F')
parser$add_argument("--marker", help = "marker gene file, split by", required=T)
args <- parser$parse_args()

suppressMessages({
library(tidyverse)
library(Seurat)
library(ggplot2)
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

violin_data <- function(seurat_object,marker){
    seurat_object$cell_type <- seurat_object@active.ident
    vln.dat <- FetchData(seurat_object,c(marker,'cell_type'))
    vln.dat.plot<-vln.dat %>% 
            reshape2::melt(,marker) %>%
            rename("Gene"="variable") %>%
            group_by(cell_type,Gene)
    vln.dat.plot$Gene <- factor(vln.dat.plot$Gene, levels = marker)
    return(vln.dat.plot)
}

plot_violin <- function(violin_data,outdir,prefix){
    p <- ggplot(data = violin_data,
    aes(x = value, y = cell_type, fill = cell_type)) +
    geom_violin(scale = 'width',
    color= 'black',size= 0.45,alpha= 0.8)+
    scale_fill_manual(values = color)+
    facet_grid(cols = vars(Gene), scales = 'free_x') +
    scale_x_continuous(breaks = seq(0, 10, by = 1)) + theme_bw()+
    theme(panel.grid = element_blank(),
    axis.title.y = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),strip.background = element_blank())
    marker_num <- length(unique(violin_data$Gene))
    type_num <- length(unique(violin_data$cell_type))
    width_fix <- marker_num*1.3
    height_fix <- type_num*0.8
    ggsave(p,filename = paste0(outdir,'/',prefix,'_violin_plot.png'),width=width_fix, height=height_fix, bg="white",dpi = 600)
    ggsave(p,filename = paste0(outdir,'/',prefix,'_violin_plot.pdf'),width=width_fix, height=height_fix, bg="white",dpi = 600)
}


main <- function(){
    createDir(args$outdir)
    seurat_object <- editRDS(args)
    marker_name <- unlist(strsplit(args$marker, split=','))
    vdata <- violin_data(seurat_object, marker_name)
    plot_violin(vdata, args$outdir, args$prefix)
}
main()