# seurat对象中，meta.data中需要有sample列，clusterorder参数不要使用，如果需要更改细胞顺序，在R中levels(seurat_object) <- c(顺序)进行修改
suppressPackageStartupMessages({
library(argparse)
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
})

parser <- function() {
    parser <- ArgumentParser()
    parser$add_argument("--seurat_object", help = "the RDS file", required = T)
    parser$add_argument("--outdir", help = "the outdir of result", default = ".")
    parser$add_argument("--prefix", help = "the prefix of result", default = "Integration")
    parser$add_argument("--features", help = "the genes need to plot, split by ','")
    parser$add_argument("--sampleorder", help = "sort the sample")
    parser$add_argument("--clusterorder", help = "sort the cluster, split by ','")
    args <- parser$parse_args()
}

create_dir <- function(dir) {
    if(!dir.exists(dir)) {
        dir.create(dir, recursive=T)
    }
}

get_eachsample_obj <- function(seurat_object) {
    meta <- seurat_object@meta.data
    if (!"sample" %in% colnames(meta)) {print("there is no col of sample in meta")}
    obj_list = list()
    for (i in unique(meta$sample)) {
        obj_list[[i]] <- subset(seurat_object, sample %in% i)
    } 
    return(obj_list)
}

get_mean_expr <- function(feature, object, clusterorder) {
    mean_expr <- data.frame(AverageExpression(object, features = feature, assays = "RNA", slot = "data")$RNA)
    if(!is.null(clusterorder)) {
        mean_expr <- mean_expr %>% select(clusterorder)
    }
    mean_expr <- t(as.matrix(mean_expr))
    return(mean_expr)
}

get_eachgene_mean <- function(feature, obj_list, sampleorder, clusterorder) {
  mean_gene_expr <- as.data.frame(map_dfr(obj_list, get_mean_expr, feature = feature, clusterorder = clusterorder))
  name_mean_expr <- data.frame(AverageExpression(obj_list[[1]], features = feature, assays = "RNA", slot = "data")$RNA)
  if(!is.null(clusterorder)) {name_mean_expr <- name_mean_expr %>% select(clusterorder)}
  rownames(mean_gene_expr) <- colnames(name_mean_expr)
  if(!is.null(sampleorder)) {
    mean_gene_expr <- mean_gene_expr %>% select(sampleorder)
  }
  return(mean_gene_expr)
}

average_heatmap <- function(mean_gene_expr, feature) {
    htdf <- scale(mean_gene_expr, scale = T, center = T)
    htCol = c("#0099CC", "white", "#CC0033")
    htRange = c(-2, 0, 2)
    col_fun <- circlize::colorRamp2(htRange, htCol)
    p <- Heatmap(htdf, name = "Z-score", cluster_columns = FALSE, cluster_rows = FALSE, row_title = feature, right_annotation = NULL,
             show_row_names = TRUE, row_names_gp = grid::gpar(fontface = "italic",fontsize = 10), row_names_side = "left", column_names_side = "bottom",
             column_names_rot = 45, col = col_fun, rect_gp = gpar(col= "white"))
}

main <- function() {
    args <- parser()
    create_dir(args$outdir)
    seurat_object <- readRDS(args$seurat_object)
    obj_list <- get_eachsample_obj(seurat_object)
    features <- unlist(strsplit(args$features, ","))
    sampleorder = NULL
    clusterorder = NULL
    if(!is.null(args$sampleorder)) {sampleorder = unlist(strsplit(args$sampleorder, ","))}
    if(!is.null(args$clusterorder)) {sampleorder = unlist(strsplit(args$clusterorder, ","))}
    for (feature in features) {
        mean_gene_expr <- get_eachgene_mean(feature, obj_list, sampleorder, clusterorder)
        write.csv(mean_gene_expr, paste0(args$outdir,"/",args$prefix, feature, ".meanexpr.csv"))
        p <- average_heatmap(mean_gene_expr, feature)
        pdf(paste0(args$outdir,"/",args$prefix, feature, ".meanheatmap.pdf"), width = 3.8*length(colnames(mean_gene_expr)), height = length(rownames(mean_gene_expr)))
        print(p)
        dev.off()
        png(paste0(args$outdir,"/",args$prefix, feature, ".meanheatmap.png"), width = 380*length(colnames(mean_gene_expr)), height = 100*length(rownames(mean_gene_expr)))
        print(p)
        dev.off()
    }
}

main()