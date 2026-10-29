#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-i", "--input_file"), type="character", default=NULL),
  make_option(c("-o", "--output_file"), type="character", default=NULL),
  make_option(c("--BC_ranking"), type="character", default=NULL),
  make_option(c("-n", "--n_clusters"), type="integer", default=NULL),
  make_option(c("--n_bcs"), type="integer", default=8),
  make_option(c("--BC_output"), type="character", default=NULL)
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

input_file <- opt$input_file
output_file <- opt$output_file
n_clust <- opt$n_clusters
n_bcs <- opt$n_bcs

suppressPackageStartupMessages({library(clue)})

calc_entropy <- function(u) {
  p <- u[u>0]
  p <- p/sum(p)
  -sum(p*log(p))
}

align_classes <- function(d, ref) {
  bcs <- as.data.frame(d)
  for(j in 1:ncol(bcs)){
    bcs[,j] <- factor(bcs[,j], levels=as.factor(as.numeric(unique(bcs[,j]))))
  }
  levs <- lapply(bcs, levels)
  n_levs <- sapply(levs, length)
  cols_to_change <- setdiff(colnames(bcs), ref)
  for(i in cols_to_change) {
    if(n_levs[i] != n_levs[ref]) {
      message(paste0(i," has ",n_levs[i]," levels; reference (",ref,") has ",n_levs[ref]," (not modifying)"))
      next
    }
    hung <- clue::solve_LSAP(table(bcs[,ref], bcs[,i]), maximum=TRUE)
    lookup <- cbind(seq_along(hung), levs[[i]][hung])
    levels(bcs[,i]) <- lookup[order(as.integer(lookup[,2])),1]
    bcs[,i] <- as.factor(as.character(bcs[,i]))
  }
  return(bcs)
}

cross_meth_entropy <- function(df){
  cm_ent_list <- lapply(colnames(df), function(c){
    df <- align_classes(d=df, ref=c)
    lv <- levels(df[[c]])
    ps <- apply(df, 1, function(u) table(factor(u, levels=lv)))
    ps <- t(ps)
    ent <- apply(ps, 1, calc_entropy)
    names(ent) <- row.names(df)
    return(ent)
  })
  mean_ent <- rowMeans(as.data.frame(cm_ent_list))
  return(mean_ent)
}

label_df <- read.delim(input_file, stringsAsFactors=FALSE, row.names=1, numerals="no.loss")

# n_clusters <- apply(label_df, 2, function(u){length(unique(u))})  ← COMMENTED OUT
# label_lists <- split(names(n_clusters), n_clusters)               ← COMMENTED OUT

message("Note: cluster number split is disabled - using first n_bcs columns")
columns <- colnames(label_df)[1:min(n_bcs, ncol(label_df))]
message(paste0("Using ", length(columns), " columns without cluster number split"))

result_ent <- cross_meth_entropy(label_df[, columns])
result_df <- as.data.frame(result_ent)
colnames(result_df) <- c("entropy")

dir.create(dirname(output_file), showWarnings=FALSE, recursive=TRUE)
write.table(result_df, file=output_file, sep="\t", col.names=NA, quote=FALSE)
message("Done!")
