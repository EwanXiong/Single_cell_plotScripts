suppressMessages({
library(tidyverse)
library(argparse)
})


parser <- ArgumentParser()
parser$add_argument("--enrich", help="sig enrich files, split by ,")
parser$add_argument("--name", help="name for each enrich file")
parser$add_argument("--top", help="top sig. [default %(default)s]", type="double", default=5)
parser$add_argument("--prefix", help="the prefix of result")
parser$add_argument("--outdir", help="the output dir")
parser$add_argument("--width", help="width", type="double")
parser$add_argument("--height", help="height", type="double")
args <- parser$parse_args()

parseArgs <- function(args){
    args$enrich <- unlist(strsplit(args$enrich,split=","))
    args$name <- unlist(strsplit(args$name,split=","))
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

plotMerge <- function(df, args){
    color <- c("#0099CC", "white", "#CC0033")
    p <- ggplot(df, aes(factor(name, levels=args$name), factor(Description, levels=rev(unique(df$Description))), color=pvalue, size=GeneRatio))+
        geom_point()+
        scale_color_gradientn(colors=color)+ 
        guides(color=guide_colorbar(title='-log10(pvalue)', title.position='top', label=T, raster=FALSE, order = 0),
                size = guide_legend(title = 'GeneRatio', title.position='top',
                title.theme = element_text(size=9, colour='black'),
                label.position='right', label.theme=element_text(size=7), order=1))+
        theme_bw()+
        labs(x = "")+
        theme(axis.title.y = element_blank(),
                axis.text.x = element_text(size=21, angle=90, hjust=1, vjust=0.5),
                axis.text.y = element_text(size=21),
                panel.grid = element_blank(),
                panel.border = element_rect(color = 'black', size=0.8, fill=NA),
                legend.position = 'right',
                # legend.position = 'top',
                # legend.direction = 'horizontal',
                legend.direction = 'vertical', 
                legend.justification = c(0.5, 0.5))

    len <- max(nchar(df$Description))
    width <- ifelse(!is.null(args$width), args$width, max(7, len/10+length(args$name)/4))
    height <- ifelse(!is.null(args$height), args$height, max(6, length(df$Description)/10+max(nchar(args$name)/5)))
    file <- paste(args$outprefix, '.enrich.buble.pdf', sep='')
    print(width)
    print(height)
    plotFunc(p, file, width=width, height=height)
}

main <- function(args){
    args <- parseArgs(args)
    createDir(args$outdir)

    df <- data.frame()
    for(i in 1:length(args$enrich)){
        data <- read.table(args$enrich[i], sep='\t', quote="", head=T)
        # print(head(data))
        data <- data %>% 
            arrange(desc(Count), pvalue) %>% 
            slice_max(Count, n=args$top, with_ties=F) %>% 
            select(Description, GeneRatio, pvalue)

        # data <- head(data, args$top)[, c("Description", "GeneRatio", "pvalue")]
        geneRatio <- sapply(data$GeneRatio, function(x){tmp <- as.numeric(unlist(strsplit(as.character(x), split="/")));tmp[1]/tmp[2]})
        data$GeneRatio <- geneRatio
        data$pvalue <- -(log10(data$pvalue))
        data$name <- args$name[i]
        df <- rbind(df, data)
    }
    plotMerge(df, args)
}
main(args)