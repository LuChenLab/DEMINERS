#' Internal functions
#'
#'
#' @rdname InterFuns
#' @name BarcodeSegment
#'
#' @importFrom rhdf5 h5ls
#' @importFrom rhdf5 h5read
#' @importFrom parallel mclapply
#' @importFrom data.table data.table
#' @importFrom smoother smth
#' @importFrom changepoint cpt.meanvar
#' @importFrom data.table :=
#'
#' @param file The file of fast5
#' @param NT The number of cores to use, i.e. at most how many child processes will be run simultaneously.

BarcodeSegment <- function(file, NT, length = 20000, MinLength = 10, ChangePoints = 98, StateStat = "Mean") {
  reads <- rhdf5::h5ls(file = file, recursive = FALSE)$name
  Fun1 <- function(read) {
    # Step1: read fast5
    h5_signals <- suppressWarnings(rhdf5::h5read(file = file,
                                                 name = paste0(read, "/Raw"),
                                                 read.attributes = TRUE))
    signal_meta <- suppressWarnings(rhdf5::h5read(file, paste0(read, "/channel_id"),
                                                  read.attributes = TRUE))
    range <- attr(signal_meta, "range")
    digitisation <- attr(signal_meta, "digitisation")
    scaling <- range/digitisation
    offset <- attr(signal_meta, "offset")
    raw_sig <- as.numeric(scaling * (as.integer(h5_signals[[1]]) + offset))

    # Step2: Get barcode signal
    if(length(raw_sig) > length) raw_sig <- raw_sig[seq_len(length)]

    Mat <- data.table::data.table(x = seq_along(raw_sig), norm = normalize_signal(raw_sig))
    Mat[, dema := smoother::smth(norm, method = "dema", n = 20)]

    cp2 <- Mat[!is.na(dema), suppressWarnings(changepoint::cpt.meanvar(dema, class = FALSE))[[1]]] + Mat[!is.na(dema), min(x)]
    ansmean <- tryCatch(suppressWarnings(changepoint::cpt.meanvar(Mat[cp2:nrow(Mat), dema], penalty = "Asymptotic", pen.value = 1e-10, method = "PELT")), error = function(e) NULL)
    if(is.null(ansmean)) return(NULL)

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
      } else {
        if(mean(Mat[round(cp2*0.75):cp2, dema] < polyA_Sig) > 0.9) {
          BCSs <- Mat[1:cp2, norm]
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
    if(is.null(BCSs)) return(NULL)

    # Step3: Segment barcode signal
    MyChangePoint(sig = BCSs, ChangePoints = ChangePoints, MinLength = MinLength, StateStat = StateStat)
  }
  parallel::mclapply(X = reads, FUN = function(x) tryCatch(Fun1(read = x), error = function(e) NULL), mc.cores = NT) -> BarcodeSegs
  names(BarcodeSegs) <- reads
  BarcodeSegs <- do.call(rbind, BarcodeSegs)
  colnames(BarcodeSegs) <- paste0("BIN", sprintf("%03d", seq_len(ncol(BarcodeSegs))))
  return(BarcodeSegs)
}


#' @rdname InterFuns
#' @name normalize_signal
normalize_signal <- function(sig) {
  med = median(sig)
  mad = median(abs(sig - med))
  (sig - med) / max(0.01, (mad * 1.4826))
}


#' @rdname InterFuns
#' @name MyChangePoint
#' @importFrom changepoint cpt.meanvar
MyChangePoint <- function(sig, MinLength = 10, ChangePoints = 68, StateStat = "Mean") {
  if(is.null(StateStat) | is.na(StateStat)) {
    stop("StateStat must be one of Mean or Median")
  }

  if(length(StateStat) != 1) {
    stop("StateStat must be one of Mean or Median")
  }

  if(!is.element(StateStat, c("Mean", "Median"))) {
    stop("StateStat must be one of Mean or Median")
  }

  cp0 <- suppressWarnings(changepoint::cpt.meanvar(data = sig,
                                                   Q = ChangePoints,
                                                   penalty = "Manual",
                                                   method = "BinSeg",
                                                   class = FALSE,
                                                   minseglen = MinLength,
                                                   param.estimates = FALSE,
                                                   pen.value = 0.0001)) - 0.5
  bins <- cut(seq_along(sig), c(0, cp0, length(sig)), include.lowest = T, labels = FALSE)

  if(StateStat == "Mean") {
    bin_sig <- as.numeric(by(sig, bins, mean))
  } else {
    bin_sig <- as.numeric(by(sig, bins, median))
  }
  return(bin_sig)
}
