#!/usr/bin/env Rscript
library(data.table)
library(DecodeR)
library(optparse)
library(parallel)

option_list = list(
  make_option(c("-I", "--fast5"), type = "character", default = NULL,
              help = "path to a fast5 file[default = %default]", metavar = "character"),
  
  make_option(c("-M", "--model"), type = "character", default = NULL,
              help = "The path of pretrained model[default = %default]", metavar = "character"),
  
  make_option(c("--cores"), type = "integer", default = 1,
              help = "the number of cores to use, i.e. at most how many child processes will be run simultaneously [default = %default]", metavar = "integer"),
  
  make_option(c("-t", "--cutoff"), type = "double", default = 0,
              help = "the cutoff of minimum probability for classified reads, high for accuracy, low for recovery [default = %default]", metavar = "double"),
  
  make_option(c("-l", "--include.lowest"), type = "logical", default = FALSE,
              help = "include.lowest of cutoff[default = %default]", metavar = "logical"),
  
  make_option(c("-O", "--output_file"), type = "character", default = NULL,
              help = "The path and file name of output[default = %default]", metavar = "integer")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

load(opt$model)
files <- list.files(opt$dir_input, full.names = T)

pred <- mclapply(files, FUN = function(x) DecodeR(fast5 = x, model = Fit1, NT = opt$cores, cutoff = opt$cutoff))
pred <- do.call(rbind, pred)
if(!dir.exists(dirname(opt$output_file))) dir.create(dirname(opt$output_file), recursive = TRUE)
fwrite(pred, opt$output_file, sep = "\t", row.names = F, quote = F)
