SJFromBam <- function(bamFile, mapqFilter = NA_integer_, txdb, genome, cores = 1, ...) {
  introns <- unique(unlist(GenomicFeatures::intronsByTranscript(txdb)))
  names(introns) <- NULL
  
  sbp <- Rsamtools::ScanBamParam(mapqFilter = mapqFilter)
  
  galp0 <- GenomicAlignments::readGAlignments(Rsamtools::BamFile(bamFile), use.names = TRUE, param = sbp)
  galp0 <- galp0[as.character(GenomicRanges::seqnames(galp0)) %in% names(genome)]
  GenomeInfoDb::seqlevels(galp0) <- levels(droplevels(GenomeInfoDb::seqnames(galp0)))
  
  galp1 <- galp0[GenomicAlignments::njunc(galp0) > 0]
  RC <- GenomicAlignments::cigar(galp1)
  
  ROverH <- OverHang(RC, cores = cores)
  ROverH <- data.table::data.table(Read = names(galp1), OverHang = ROverH)
  data.table::setkey(ROverH, Read)
  ROverH <- unique(ROverH)
  ROverH <- ROverH[, .(OverHang = min(OverHang)), Read]
  
  SJ0 <- unlist(GenomicAlignments::junctions(galp1))
  
  mmr <- data.table::data.table(ReadID = names(galp0))
  data.table::setkey(mmr, ReadID)
  mmr <- mmr[, .N, ReadID]
  mmp <- mmr[, .(Rs = .N), N]
  mmp$P <- prop.table(mmp$Rs)
  
  galp1 <- galp0[!names(galp0) %in% mmr[N == mmp[P > 0.5, min(N)], ReadID]]
  galp0 <- galp0[names(galp0) %in% mmr[N == mmp[P > 0.5, min(N)], ReadID]]
  
  galp0 <- galp0[GenomicAlignments::njunc(galp0) > 0]
  galp1 <- galp1[GenomicAlignments::njunc(galp1) > 0]
  
  # uniquely-mapping reads
  juncs <- data.table::data.table(as.data.frame(unlist(GenomicAlignments::junctions(galp0), use.names = FALSE)))
  juncs <- juncs[, .N, by = c("seqnames", "start", "end")]
  introns <- data.table::data.table(as.data.frame(introns))
  juncs1 <- merge(juncs, introns, by = c("seqnames", "start", "end"))
  
  juncs0 <- juncs[!paste0(seqnames, ":", start, "-", end) %in% juncs1[, paste0(seqnames, ":", start, "-", end)]]
  juncs0_1 <- juncs0[, GenomicRanges::GRanges(seqnames = seqnames, ranges = IRanges::IRanges(start = start, end = end), strand = "+")]
  juncs0$Motif1 <- SpliceMotif(x = juncs0_1, fa = genome)
  
  juncs0[Motif1 %in% c("GTAG", "GCAG", "ATAC"), strand := "+"]
  juncs0[Motif1 %in% c("CTAC", "CTGC", "GTAT"), strand := "-"]
  juncs0[is.na(strand), strand := "*"]
  
  juncs <- rbind(juncs0[, .(seqnames, start, end, strand, N, annotation = 0)], juncs1[, .(seqnames, start, end, strand, N, annotation = 1)])
  juncs$Motif <- SpliceMotif(x = juncs[, GenomicRanges::GRanges(seqnames = seqnames, ranges = IRanges::IRanges(start = start, end = end), strand = "+")], fa = genome)
  juncs$Motif <- plyr::mapvalues(juncs$Motif, c("GTAG", "CTAC", "GCAG", "CTGC", "ATAC", "GTAT"), 1:6)
  juncs[!Motif %in% 1:6, Motif := 0]
  juncs[, strand := plyr::mapvalues(strand, c("*", "+", "-"), 0:2)]
  data.table::setnames(juncs, "N", "umr")
  juncs_um <- data.table::copy(juncs)
  
  # multi-mapping reads
  juncs <- data.table::data.table(as.data.frame(unlist(GenomicAlignments::junctions(galp1), use.names = FALSE)))
  juncs <- juncs[, .N, by = c("seqnames", "start", "end")]
  introns <- data.table::data.table(as.data.frame(introns))
  juncs1 <- merge(juncs, introns, by = c("seqnames", "start", "end"))
  
  juncs0 <- juncs[!paste0(seqnames, ":", start, "-", end) %in% juncs1[, paste0(seqnames, ":", start, "-", end)]]
  juncs0_1 <- juncs0[, GenomicRanges::GRanges(seqnames = seqnames, ranges = IRanges::IRanges(start = start, end = end), strand = "+")]
  juncs0$Motif1 <- SpliceMotif(x = juncs0_1, fa = genome)
  
  juncs0[Motif1 %in% c("GTAG", "GCAG", "ATAC"), strand := "+"]
  juncs0[Motif1 %in% c("CTAC", "CTGC", "GTAT"), strand := "-"]
  juncs0[is.na(strand), strand := "*"]
  
  juncs <- rbind(juncs0[, .(seqnames, start, end, strand, N, annotation = 0)], juncs1[, .(seqnames, start, end, strand, N, annotation = 1)])
  juncs$Motif <- SpliceMotif(x = juncs[, GenomicRanges::GRanges(seqnames = seqnames, ranges = IRanges::IRanges(start = start, end = end), strand = "+")], fa = genome)
  juncs$Motif <- plyr::mapvalues(juncs$Motif, c("GTAG", "CTAC", "GCAG", "CTGC", "ATAC", "GTAT"), 1:6)
  juncs[!Motif %in% 1:6, Motif := 0]
  juncs[, strand := plyr::mapvalues(strand, c("*", "+", "-"), 0:2)]
  data.table::setnames(juncs, "N", "mmr")
  juncs_mm <- data.table::copy(juncs)
  
  juncs <- merge(juncs_um, juncs_mm, by = c("seqnames", "start", "end", "strand", "annotation", "Motif"), all = TRUE)
  juncs[is.na(juncs)] <- 0
  
  gr <- juncs[, GenomicRanges::GRanges(seqnames, IRanges::IRanges(start, end))]
  ov <- GenomicRanges::findOverlaps(gr, SJ0, type = "equal")
  ov <- data.table::data.table(sj = S4Vectors::queryHits(ov), Read = names(SJ0)[S4Vectors::subjectHits(ov)])
  ov <- merge(ov, ROverH, by = "Read")
  data.table::setkey(ov, sj)
  ov <- ov[, .(OverHang = max(OverHang)), sj]
  
  juncs[ov$sj, OverHang := ov$OverHang]
  juncs <- juncs[, .(seqnames, start, end, strand, Motif, annotation, umr, mmr, OverHang)]
  return(juncs)
}


SpliceMotif <- function(x, fa) {
  donorSite <- GenomicRanges::resize(x, width = 2, fix = "start")
  donorMotif <- suppressWarnings(BSgenome::getSeq(fa, donorSite))
  
  acceptorSite <- GenomicRanges::resize(x, width = 2, fix = "end")
  acceptorMotif <- suppressWarnings(BSgenome::getSeq(fa, acceptorSite))
  
  paste0(donorMotif, acceptorMotif)
}

OverHang <- function(x, cores = 1) {
  CRL <- stringr::str_extract_all(x, "([[:digit:]]+)M")
  parallel::mcmapply(CRL, FUN = function(x) {
    min(as.numeric(gsub("M", "",x[c(1, length(x))])))
  }, mc.cores = cores)
}

library(GenomicFeatures)
library(Biostrings)
library(data.table)
gtfFile <- "/mnt/raid61/Personal_data/songjunwei/reference/plasmoDB/PlasmoDB-53_PbergheiANKA.gff"
txdb <- GenomicFeatures::makeTxDbFromGFF(gtfFile)
genome <- readDNAStringSet("/mnt/raid61/Personal_data/songjunwei/reference/plasmoDB/PlasmoDB-53_PbergheiANKA_Genome.fasta")
names(genome) <- mapply(function(x) x[1], strsplit(names(genome), " \\| "))

RTA03 <- SJFromBam(bamFile = "/mnt/raid61/Personal_data/tangchao/Temp/20211025/05.TranscriptClean/RTA-03.minimap2genome.bam", txdb = txdb, genome = genome)
bams <- list.files("/mnt/raid61/Personal_data/tangchao/Temp/20211025/05.TranscriptClean", "bam$", full.names = T)

lapply(bams, function(x) {
  SJ <- SJFromBam(bamFile = x, txdb = txdb, genome = genome)
  fwrite(SJ, gsub("bam$", "SJ.out.tab", x), sep = "\t", quote = F, col.names = F, row.names = F)
})

library(ICAS)

SJs <- list.files("/mnt/raid61/Personal_data/tangchao/Temp/20211025/05.TranscriptClean", "SJ.out.tab$", full.names = T)

meta <- data.frame(Stage = c("Trophozoite", "Trophozoite", "Trophozoite", "Schizont", "Schizont", "Schizont"), 
                   row.names = c("RTA-03", "RTA-10", "RTA-16", "RTA-17", "RTA-24", "RTA-32"))

MyObj <- ICASDataSetFromSJFile(SJFiles = SJs, 
                               postfix = ".minimap2genome.SJ.out.tab",
                               colData = meta,
                               uniqMapOnly = TRUE,
                               annotationOnly = FALSE,
                               minSJRowSum = 10,
                               minSJRead = 1,
                               minOverhang = 3,
                               design = "Stage")

MyObj <- PSICalculater(object = MyObj, MMJF = 0.01, MinSumSpliceSite = 3)
psi(MyObj)[1:4, ]

DE1 <- FindAllMarkers(object = MyObj, 
                      only.pos = T, 
                      delta.threshold = 0,
                      return.thresh = 1,
                      test.use = "tobit",
                      print.bar = FALSE,
                      NT = 20)


DE2 <- FindAllMarkers(object = MyObj, 
                      only.pos = T, 
                      delta.threshold = 0,
                      return.thresh = 1,
                      test.use = "wilcox",
                      print.bar = FALSE,
                      NT = 20)


