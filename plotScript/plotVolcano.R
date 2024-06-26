suppressMessages({
library(argparse)
library(ggpubr)
library(ggrepel)
library(tidyverse)
})

parser <- ArgumentParser()
parser$add_argument('--outdir', help='outdir of result, [default %(default)s]', default='.')
parser$add_argument('--prefix', help='prefix of outdir, [default %(default)s]', default='project')
parser$add_argument("--method", help="p_val or p_val_adj or pct. [default %(default)s]", default="p_val_adj")
parser$add_argument("--diff", help="diff expressed genes file from seurat")
parser$add_argument("--top", help="top genes to show. [default %(default)s]", default=5)
parser$add_argument("--gene", help="genes to show, split by ,")
parser$add_argument("--title", help="plot title")
args <- parser$parse_args()

createDir <- function(dir){
    if(!dir.exists(dir)){
        dir.create(dir, recursive=T)
    }
}

plotVolcano <- function(diff, args){
    x_lim <- max(abs(diff$avg_log2FC))
    palette <- c('#FC4E07', '#D3D3D3', '#00AFBB')
    names(palette) <- c('Up', 'Not significant', 'Down')
    title <- args$title
    ylab <- ifelse(args$method=='p_val', "-log10(pval)", ifelse(args$method=='p_val_adj', "-log10(p_val_adj)", "abs(pct.1 - pct.2)"))
    up <- diff[!is.na(diff$label) & diff$type=="Up", ]
    down <- diff[!is.na(diff$label) & diff$type=="Down", ]
    if(sum(diff[,'p_val'])==0)return(NULL)

    p <- ggplot(diff, aes(avg_log2FC, value, color = type)) +
        geom_point(size = 1) + labs(y=ylab, title=args$title) + 
        scale_color_manual(values = palette)+
        theme_classic()+
        geom_vline(xintercept = c(-0.25, 0.25), linetype = 'dashed')+
        geom_text_repel(aes(label = label), color = "black", data=diff, size=4)+
        xlim(-x_lim, x_lim)+
        theme(legend.title = element_blank(),
            legend.text = element_text(size = 12),
            plot.title = element_text(face = 'bold', size = 14, vjust = .5, hjust = .5),
            axis.title = element_text(face = 'bold', size = 12, vjust = .5, hjust = .5),
            plot.background=element_blank(),
            panel.background=element_rect(fill='transparent', color='black'))+
            guides(color = guide_legend(override.aes = list(size = 4), ncol=1))+
            geom_point(aes(avg_log2FC, value), data = up, color = '#FC4E07', size = 2)+
            geom_point(aes(avg_log2FC, value), data = down, color = '#00AFBB', size = 2)

    pdf(paste0(args$outdir, '/', args$prefix, '.volcano.pdf'))
    print(p)
    dev.off()

    png(paste0(args$outdir, '/', args$prefix, '.volcano.png'))
    print(p)
    dev.off()
}

main <- function(){
    createDir(args$outdir)
    diff <- read.table(args$diff, header=T, row.names=1, stringsAsFactors=F)
    diff$geneID <- rownames(diff)
    diff$type <- 'Not significant'
    if(args$method=="p_val"){
        diff$type[which((diff$avg_log2FC > 0.25) & (diff$p_val < 0.05  ))] <- 'Up'
        diff$type[which((diff$avg_log2FC < -0.25) & (diff$p_val < 0.05) )] <- 'Down'
    }else{
        diff$type[which((diff$avg_log2FC > 0.25) & (diff$p_val_adj < 0.05  ))] <- 'Up'
        diff$type[which((diff$avg_log2FC < -0.25) & (diff$p_val_adj < 0.05) )] <- 'Down'
    }

    if(args$method=="p_val" | args$method=='p_val_adj'){
        diff$value <- -log10(diff[, args$method])
    }else{
        diff$value <- abs(diff$pct.1 - diff$pct.2)
    }

    if(!is.null(args$gene)){
        gene <- unlist(strsplit(args$gene, split = ',')) %>% 
            intersect(diff$geneID) %>% 
            unique()
    }else{
        diff <- diff %>% arrange(avg_log2FC)
        up <- diff %>% filter(type=="Up") %>% select(geneID) %>% tail(as.numeric(args$top))
        down <- diff %>% filter(type=="Down") %>% select(geneID) %>% head(as.numeric(args$top))
        up_down <- rbind(up, down)
        if(nrow(up_down)>0){
            gene <- unique(up_down$geneID)
        }else{
            gene <- NULL
            print("there is no significant genes")
        }
    }

    diff$label <- NA
    diff$label[match(gene, diff$geneID)] <- gene

    plotVolcano(diff, args)
}
main()
