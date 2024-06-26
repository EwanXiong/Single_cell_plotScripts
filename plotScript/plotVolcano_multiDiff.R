suppressPackageStartupMessages({
library(tidyverse)
library(reshape2)
library(argparse)
library(ggrepel)
})

parser <- ArgumentParser()
parser$add_argument("--outdir", help="",default=".")
parser$add_argument("--prefix", help="",default="multi")
parser$add_argument("--diff", help="diff file, split by , contains columns gene_id, avg_log2FC and p_val_adj")
parser$add_argument("--cluster", help="cluster name, split by , number equal to diff")
parser$add_argument("--col", help="col number of gene_id, avg_log2FC, p_val_adj. [default %(default)s]", default='1,3,6')
parser$add_argument("--type", help="up or down or both. [default %(default)s]", default='both')
parser$add_argument("--top", help="top diff gene label, 0 was not label. [default %(default)s]", default=5)
parser$add_argument("--labelsize", help="label size of gene. [default %(default)s]", type="double", default=4)
parser$add_argument("--avg_log2FC", help="avg_log2FC. [default %(default)s]", type="double", default=0.25)
parser$add_argument("--p_val_adj", help="p_val_adj. [default %(default)s]", type="double", default=0.05)
args <- parser$parse_args()

color <- c("#00CC99", "#F0E685", "#99CC00", "#C75127", "#5050FF", "#CDDEB7","#EE82EE",
"#924822","#FFC20A","#33CC00","#0099CC","#7A65A5","#AE1F63","#D595A7",
"#A9A9A9","#E4AF69","#990080","#14FFB1","#003399", "#FF1463","#D60047",
"#CE3D32","#996600","#809900","#008099","#0A47FF","#660099","#FFD147",
"#339900","#5A655E","#4775FF","#802268","#8B0000","#BBFFFF","#00868B",
"#00008B","#CDBA96","#FF69B4","#3B1B53","#D2B48C","#00FFFF","#836FFF",
"#BA6338","#991A00","#D58F5C","#E7C76F","#CC9900","#749B58","#6BD76B",
"#009966","#00D68F","#5DB1DD","#466983","#837B8D","#612A79","#990033"
)

parseArgs <- function(args){
    args$diff <- unlist(strsplit(args$diff, split=","))
    args$cluster <- unlist(strsplit(args$cluster, split=","))
    args$col <- as.numeric(unlist(strsplit(args$col, split=",")))
    args$outprefix <- file.path(args$outdir, args$prefix)
    return(args)
}

createDir <- function(dir){
    if(!dir.exists(dir)){
        dir.create(dir, recursive=T)
    }
}

plotFunc <- function(pic, file, width=7, height=7, res=300){
    pdf(file, width=width, height=height )
    print(pic)
    dev.off()
    png <- gsub('pdf$', 'png', file)
    png(png, width=width, height=height, res=res, units='in')
    print(pic)
    dev.off()
}

plotVolcano <- function(df, top, limit, args){
    pointColor <- c('#FC4E07', '#00AFBB', '#696969')
    # pointColor <- c("Firebrick", "NavyBlue", "#696969")
    if(args$type=='up'){
        pointColor <- pointColor[c(1, 3)]
    }
    if(args$type=='down'){
        pointColor <- pointColor[c(2, 3)]
    }
    if(nrow(top)==0){
        pointColor <- pointColor[3]
    }
    file <- paste0(args$outprefix, '_volcano.pdf')

    filterDF <- filter(df, top==0)
    filterDF$cluster <- factor(filterDF$cluster, levels=args$cluster)
    limit$x <- factor(limit$x, levels=args$cluster)
    p <- ggplot()
    if(args$type=='up'){p <- p + geom_col(data=limit, aes(limit[,1], limit[,4]), fill = "#D3D3D3")}
    if(args$type=='down'){p <- p + geom_col(data=limit, aes(limit[,1], limit[,3]), fill = "#D3D3D3")}
    if(args$type=='both'){p <- p + geom_col(data=limit, aes(limit[,1], limit[,3]), fill = "#D3D3D3") + geom_col(data=limit, aes(limit[,1], limit[,4]), fill = "#D3D3D3")}

    p <- p + geom_jitter(data=filterDF, aes(cluster, avg_log2FC, color=type), shape=19, size=0.5, width=0.4)
    if(args$top!=0 & nrow(top)!=0){
        p <- p + geom_jitter(data=top, aes(cluster, avg_log2FC, color=type), shape=19, size=2, width=0.4) + 
            geom_text_repel(data=top, aes(cluster, avg_log2FC, label=gene_id), size=args$labelsize, force=1)
    }

    p <- p + geom_tile(data=limit, aes(limit[,1], 0, fill=limit[,1]), height=args$avg_log2FC*1.8) + 
        geom_text(data=limit, aes(limit[,1], 0, label=1:length(limit[,1])), size=5, color ="white")

    p <- p + scale_color_manual(values = pointColor) + 
        scale_fill_manual(values = as.character(limit$color)) +
        labs(x=NULL, y="average log2FC", color=NULL, fill=NULL)+ 
        theme_minimal()+  
        theme(axis.line.y=element_line(color = "black", size = 0.8), axis.title.y=element_text(size = 18), axis.text.y=element_text(size = 13), 
            axis.line.x=element_blank(), axis.text.x = element_blank(), panel.grid = element_blank(), legend.position = "right", legend.direction = "vertical", 
            legend.text=element_text(size = 15), aspect.ratio=1) + 
        guides(color = guide_legend(override.aes = list(size = 6)), fill = guide_legend(override.aes = list(size = 6)))

    width <- ifelse(nrow(limit)<4, 7, (nrow(limit)-4)/5+7)
    plotFunc(p, file, width=width, height=7)
}

main <- function(args){
    args <- parseArgs(args)
    createDir(args$outdir)
    df <- data.frame()
    top <- data.frame()
    ymin <-c()
    ymax <- c()
    for(i in 1:length(args$diff)){
        diff <- read.table(args$diff[i], head=T, sep='\t', quote="", stringsAsFactors=F)[, args$col]
        colnames(diff) <- c("gene_id", "avg_log2FC", "p_val_adj")
        if(args$type=="up"){
            diff <- diff[diff$avg_log2FC>0, ]
        }
        if(args$type=="down"){
            diff <- diff[diff$avg_log2FC<0, ]
        }
        ymin <- c(ymin, min(diff$avg_log2FC))
        ymax <- c(ymax, max(diff$avg_log2FC))
        diff <- diff[abs(diff$avg_log2FC) >= args$avg_log2FC, ]

        diff$type <- ifelse(diff$avg_log2FC >= args$avg_log2FC & diff$p_val_adj < args$p_val_adj, 'UP', ifelse(diff$avg_log2FC <= (-args$avg_log2FC) & diff$p_val_adj < args$p_val_adj, 'DOWN', 'Not Significant'))
        # print(unique(diff$type))
        diff$cluster <- args$cluster[i]
        diff$top <- 0

        if(args$top!=0){
            topClusterUp <- diff[diff$type=='UP', ]
            if(nrow(topClusterUp)>0){
                topClusterUp <- top_n(topClusterUp, args$top, abs(topClusterUp$avg_log2FC))
            }else{
                topClusterUp <- data.frame()
            }

            topClusterDown <- diff[diff$type=='DOWN', ]
            if(nrow(topClusterDown)>0){
                topClusterDown <- top_n(topClusterDown, args$top, abs(topClusterDown$avg_log2FC))
            }else{
                topClusterDown <- data.frame()
            }

            topCluster <- rbind(topClusterUp, topClusterDown)
            if(nrow(topCluster)!=0){
                diff$top <- case_when(!(diff$gene_id %in% topCluster$gene_id) ~ 0, diff$gene_id %in% topCluster$gene_id ~ 1)
                topCluster$top <- 1
            }
            top <- rbind(top, topCluster)
        }
        df <- rbind(df, diff)
    }
    level <- c('UP', 'DOWN', 'Not Significant')
    df$type <- factor(df$type, levels=level[level %in% df$type])
    df$cluster <- factor(df$cluster, levels=args$cluster)
    top$type <- factor(top$type, levels=level[level %in% top$type])

    limit <- data.frame(x=args$cluster, color=color[1:length(args$cluster)], ymin=ymin, ymax=ymax)
    plotVolcano(df, top, limit, args)
}
main(args)