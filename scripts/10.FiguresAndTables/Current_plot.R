library(PorexploreR)
library(data.table)
library(rhdf5)
library(dplyr)
library(ggplot2)
library(smoother)
library(parallel)
library(spatstat)
library(Cairo)
library(changepoint)

readFast5 <- function(file, NT) {
  reads <- rhdf5::h5ls(file = file, recursive = FALSE)$name
  parallel::mclapply(X = reads, FUN = function(read) {
    h5_signals <- rhdf5::h5read(file = file, name = paste0(read, "/Raw"), read.attributes = TRUE)
    signal_meta <- rhdf5::h5read(file, paste0(read, "/channel_id"), read.attributes = TRUE)
    
    range <- attr(signal_meta, "range")
    
    digitisation <- attr(signal_meta, "digitisation")
    scaling <- range/digitisation
    offset <- attr(signal_meta, "offset")
    sampling_rate <- attr(signal_meta, "sampling_rate")
    
    new("Squiggle", raw_signal = as.integer(h5_signals[[1]]), 
        range = as.numeric(range), digitisation = as.integer(digitisation), 
        offset = as.integer(offset), sampling_rate = as.integer(sampling_rate), 
        scaling = as.numeric(scaling))
  }, mc.cores = NT)
}

normalize_signal <- function(sig) {
  med = median(sig)
  mad = median(abs(sig - med))
  (sig - med) / max(0.01, (mad * 1.4826))
}

barcode_scale <- function(s, len = 100) {
  if(length(s) < len) {
    times <- ceiling(len/length(s))
    s <- rep(s, each = times)
  }
  
  b <- as.numeric(cut(seq_along(s), breaks = len, include.lowest = T))
  se <- as.numeric(by(s, b, FUN = median))
  return(se)
}

BinCollapse <- function(sig, BinFun = "median") {
  steps <- diff(c(sig[1], sig))
  ir <- data.table::as.data.table(sort(c(IRanges::IRanges(steps > 0), IRanges::IRanges(steps <= 0))))
  data.table::setnames(ir, "width", "width_up")
  ir$width_down <- c(ir[-1, width_up], 0)
  
  ir[, left := end - width_up/2]
  ir[, right := end + width_down/2]
  
  bin <- cut(seq_along(sig), breaks = c(0, ir$left, length(sig)), labels = FALSE)
  dre <- c(1, ifelse(as.numeric(S4Vectors::runValue(S4Vectors::Rle(steps > 0))) == 1, 1, -1))
  
  if(BinFun == "median") {
    res <- as.numeric(do.call(c, by(sig, bin, FUN = function(x) rep(median(x), length(x)))))
  }
  
  if(BinFun == "mean") {
    res <- as.numeric(do.call(c, by(sig, bin, FUN = function(x) rep(mean(x), length(x)))))
  }
  
  if(BinFun == "max") {
    binmax <- as.numeric(by(sig, bin, FUN = function(x) max(x)))
    binmin <- as.numeric(by(sig, bin, FUN = function(x) min(x)))
    binmax[which(dre != 1)] <- binmin[which(dre != 1)]
    res <- rep(binmax, S4Vectors::runLength(S4Vectors::Rle(bin)))
  }
  return(res)
}

GetBarcode3 <- function(read, length = 20000, plot = TRUE, col = "red") {
  raw_sig <- PorexploreR::signal(read)
  if(length(raw_sig) > length) raw_sig <- raw_sig[seq_len(length)]
  
  Mat <- data.table::data.table(x = seq_along(raw_sig), norm = normalize_signal(raw_sig))
  Mat[, dema := smoother::smth(norm, method = "dema", n = 20)]
  
  cp2 <- Mat[!is.na(dema), suppressWarnings(changepoint::cpt.meanvar(dema, class = FALSE))[[1]]] + Mat[!is.na(dema), min(x)]
  
  ansmean <- suppressWarnings(changepoint::cpt.meanvar(Mat[cp2:nrow(Mat), dema], penalty = "Asymptotic", pen.value = 1e-10, method = "PELT"))
  whichBin <- which.max(diff(c(0, head(ansmean@cpts, 4))))
  polyA_Pos <- ifelse(whichBin == 1, cp2, ansmean@cpts[whichBin - 1] + cp2)
  polyA_Sig <- ansmean@param.est$mean[whichBin]
  
  if(polyA_Pos < 2500) {
    BCSs <- NULL
  } else {
    if(mean(Mat[round(polyA_Pos*0.75):polyA_Pos, dema] < polyA_Sig) > 0.9) {
      # Mat[1:polyA_Pos,  smth := smoother::smth(norm, method = "dema", n = round(polyA_Pos/400), v = 1)]
      # BCSs <- Mat[!is.na(smth), smth]
      BCSs <- Mat[1:polyA_Pos, norm]
      if(plot) {
        Mat[1:polyA_Pos,  smth := smoother::smth(norm, method = "dema", n = round(polyA_Pos/400), v = 1)]
        Mat[cp2:nrow(Mat), smth := NA]
        return(Mat)
      }
    } else {
      if(mean(Mat[round(cp2*0.75):cp2, dema] < polyA_Sig) > 0.9) {
        BCSs <- Mat[1:cp2, norm]
        if(plot) {
          Mat[1:polyA_Pos, smth := smoother::smth(norm, method = "dema", n = round(polyA_Pos/400), v = 1)]
          Mat[cp2:nrow(Mat), smth := NA]
          return(Mat)
        }
      } else {
        BCSs <- NULL
      }
    }
  }
  
  if(!is.null(BCSs)) {
    if(length(BCSs) < 2100) {
      BCSs <- NULL
    }
  }
  return(BCSs)
}


read2gene <- fread("/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/analysis/02.BamReadsSplit/read2gene_20201127.txt")

dir_f5 <- "/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/Rawdata/20201127/20201126/no_sample/20201126_1357_MN26652_FAO61922_e7fbd202/fast5_pass"
file1 <- file.path(dir_f5, list.files(dir_f5)[1])
fast5_1 <- readFast5(file = file1, NT = 2)
names(fast5_1) <- gsub("read_", "", h5ls(file = file1, recursive = FALSE)$name)

read2gene[read %in% names(fast5_1), .N, Barcode]
read2gene2U <- read2gene[read %in% names(fast5_1), ]


M1 <- GetBarcode3(read = fast5_1[[read2gene2U[Barcode == "RTA-03", read][2]]], plot = T, length = 20000)
M2 <- GetBarcode3(read = fast5_1[[read2gene2U[Barcode == "RTA-08", read][1]]], plot = T, length = 20000)
M3 <- GetBarcode3(read = fast5_1[[read2gene2U[Barcode == "RTA-10", read][22]]], plot = T, length = 20000)
M4 <- GetBarcode3(read = fast5_1[[read2gene2U[Barcode == "RTA-21", read][5]]], plot = T, length = 20000)
M5 <- GetBarcode3(read = fast5_1[[read2gene2U[Barcode == "RTA-27", read][3]]], plot = T, length = 20000)
M6 <- GetBarcode3(read = fast5_1[[read2gene2U[Barcode == "RTA-33", read][18]]], plot = T, length = 20000)
M7 <- GetBarcode3(read = fast5_1[[read2gene2U[Barcode == "RTA-37", read][6]]], plot = T, length = 20000)

save(M1, M2, M3, M4, M5, M6, M7, file = "/mnt/raid61/Personal_data/tangchao/Temp/Ms.RData")


library(ggplot2)
library(data.table)
load("/mnt/raid61/Personal_data/tangchao/Temp/Ms.RData")

ggplot(M1) + 
  geom_line(aes(x = x, y = norm)) + 
  geom_line(aes(x = x, y = smth), color = "#E41A1C") + 
  theme_light() + 
  labs(y = "Scaled current") + 
  theme(axis.title.x = element_blank()) -> p1


ggplot(M2) + 
  geom_line(aes(x = x, y = norm)) + 
  geom_line(aes(x = x, y = smth), color = "#377EB8") + 
  theme_light() + 
  labs(y = "Scaled current") + 
  theme(axis.title.x = element_blank()) -> p2


ggplot(M3) + 
  geom_line(aes(x = x, y = norm)) + 
  geom_line(aes(x = x, y = smth), color = "#4DAF4A") + 
  theme_light() + 
  labs(y = "Scaled current") + 
  theme(axis.title.x = element_blank()) -> p3


ggplot(M4) + 
  geom_line(aes(x = x, y = norm)) + 
  geom_line(aes(x = x, y = smth), color = "#984EA3") + 
  theme_light() + 
  labs(y = "Scaled current") + 
  theme(axis.title.x = element_blank()) -> p4


ggplot(M5) + 
  geom_line(aes(x = x, y = norm)) + 
  geom_line(aes(x = x, y = smth), color = "#984EA3") + 
  theme_light() + 
  labs(y = "Scaled current") + 
  theme(axis.title.x = element_blank()) -> p5


ggplot() + 
  geom_line(data = M6, aes(x = x, y = norm)) + 
  geom_line(data = M6[!is.na(smth)], aes(x = x, y = smth), color = "#FF7F00") + 
  theme_light() + 
  labs(y = "Scaled current") + 
  theme(axis.title.x = element_blank()) -> p6


ggplot(M7) + 
  geom_line(aes(x = x, y = norm)) + 
  geom_line(aes(x = x, y = smth), color = "#A65628") + 
  theme_light() + 
  labs(y = "Scaled current") + 
  theme(axis.title.x = element_blank()) -> p7


library(patchwork)
p1 / p2 / p3 / p5 / p6 / p7
ggsave("/mnt/raid61/Personal_data/tangchao/Temp/Ms.pdf", width = 10, height = 8)

