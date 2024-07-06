#' Barcode Prediction
#'
#' @rdname DecodeR
#' @name DecodeR
#'
#' @param fast5 path to a fast5 file
#' @param NT the number of cores to use, i.e. at most how many child processes will be run simultaneously.
#' @param model a trained model for prediction
#' @param cutoff the cutoff of minimum probability for classified reads, high for accuracy, low for recovery
#' @param include.lowest include.lowest of cutoff
#' @param clear TRUE for a clear output
#'
#' @importFrom data.table as.data.table
#' @importFrom data.table melt.data.table
#' @importFrom data.table .SD
#' @importFrom data.table :=
#' @export
#'
#' @examples
#' fast5file <- system.file("extdata/demo2_0.fast5", package = "DecodeR")
#' data("Model_2barcodes")
#'
#' pred <- DecodeR(fast5 = fast5file, model = Model_2barcodes)
#' hist(pred$Probability, xlab = "Probability", main = "Histogram of Probability")
#' table(pred$Barcode)
#'
#' pred2 <- DecodeR(fast5 = fast5file, model = Model_2barcodes, cutoff = 0.8)
#' table(pred2$Barcode)

DecodeR <- function(fast5, model = NULL, NT = 1, cutoff = 0, include.lowest = FALSE, clear = FALSE) {
  barcodeSigsBinMat <- BarcodeSegment(file = fast5, NT = NT)
  if(is.null(model)) return(barcodeSigsBinMat)
  pred0 <- predict(model, barcodeSigsBinMat)
  pred1 <- predict(model, barcodeSigsBinMat, type = "prob")
  # pred <- as.data.table(pred, keep.rownames = "Read")
  # PPs <- melt.data.table(pred, id.vars = "Read", variable.name = "Barcode", value.name = "Probability")
  # PPs <- PPs[, .SD[which.max(Probability), ], Read]
  PPs <- data.table::data.table(Read = rownames(barcodeSigsBinMat), Barcode = pred0, Probability = apply(pred1, 1, max))
  if(include.lowest) {
    PPs[Probability < cutoff, Barcode := "unclassified"]
  } else {
    PPs[Probability <= cutoff, Barcode := "unclassified"]
  }
  if(clear) {
    return(PPs)
  } else {
    return(data.table::data.table(PPs, data.table(pred1)))
  }
}

# DecodeR <- function(fast5, model, NT = 1, cutoff = 0, include.lowest = FALSE, clear = FALSE) {
#   barcodeSigsBinMat <- BarcodeSegment(file = fast5, NT = NT)
#   pred <- predict(model, barcodeSigsBinMat, type = "prob")
#   pred <- as.data.table(pred, keep.rownames = "Read")
#   PPs <- melt.data.table(pred, id.vars = "Read", variable.name = "Barcode", value.name = "Probability")
#   PPs <- PPs[, .SD[which.max(Probability), ], Read]
#   if(include.lowest) {
#     PPs[Probability < cutoff, Barcode := "unclassified"]
#   } else {
#     PPs[Probability <= cutoff, Barcode := "unclassified"]
#   }
#   if(clear) {
#     return(PPs)
#   } else {
#     return(merge(pred, PPs, by = "Read"))
#   }
# }
