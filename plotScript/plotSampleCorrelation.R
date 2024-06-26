suppressPackageStartupMessages({
library(Seurat)
library(corrplot)
library(argparse)
})

parser <- ArgumentParser()
parser$add_argument('--seurat_object', help='RDS file')
parser$add_argument('--outdir', help='outdir of result, [default %(default)s]', default='.')
parser$add_argument('--prefix', help='prefix of outdir, [default %(default)s]', default='program')
args <- parser$parse_args()

createDir <- function(dir){
    if(!dir.exists(dir)){
        dir.create(dir, recursive=T)
    }
}

plotCor <- function(data, prefix){
    pdf(paste(prefix,'.corrplot.pdf',sep=''))
    corrcol <- colorRampPalette(c("red","orange","blue","white","white"))
    corrplot(data, cl.lim = c(0, 1),tl.col = "black",col = rev(corrcol(50)))
    dev.off()
    png(paste(prefix,'.corrplot.png',sep=''))
    corrplot(data, cl.lim = c(0, 1),tl.col = "black",col = rev(corrcol(50)))
    dev.off()
    write.table(data, paste0(prefix,'.corrplot.xls'), sep="\t", quote=F, col.names=NA)
}

corSample <- function(obj, prefix){
    variable_gene <- VariableFeatures(obj)

    Idents(obj) <- obj$sample
    average <- AverageExpression(object = obj)$RNA
    M <- cor(average)
    order.hc <- corrMatOrder(M, order = "hclust")
    M.hc <- M[order.hc, order.hc]
    print(M.hc)
    plotCor(M.hc, paste0(prefix, "_allGene"))

    M <- cor(average[variable_gene, ])
    order.hc <- corrMatOrder(M, order = "hclust")
    M.hc <- M[order.hc, order.hc]
    print(M.hc)
    plotCor(M.hc, paste0(prefix, "_variableGene"))
}

main <- function(){
    createDir(args$outdir)
    seurat_object <- readRDS(args$seurat_object)
    corSample(seurat_object, file.path(args$outdir, args$prefix))
}
main()