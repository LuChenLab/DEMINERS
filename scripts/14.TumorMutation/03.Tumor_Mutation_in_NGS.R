library(data.table)
library(rtracklayer)
library(GenomicAlignments)

Candi <- openxlsx::read.xlsx("/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/analysis/14.TumorMutation/Tumor_specific_mutation.xlsx")
setDT(Candi)
Candi[, seq := NULL]

chain <- suppressWarnings(import.chain(con = "/mnt/raid61/Personal_data/tangchao/Document/UCSC/goldenPath/hg19/liftOver/hg19ToHg38.over.chain"))

gr0 <- Candi[, GRanges(seqnames = CHROM, IRanges(start = POS, end = POS))]
gr1 <- suppressWarnings(rtracklayer::liftOver(x = gr0, chain = chain))
gr1 <- unlist(gr1)
names(gr1) <- as.character(gr0)


# DMG-1


bamfile <- file.path("/mnt/raid64/visium/analysis_ONT/merge_fq/NGS_spaceOut/B1_DMG_1_hs/outs/possorted_genome_bam.bam")
file.exists(bamfile)

DMG_1 <- lapply(X = seq_along(gr1), function(i) {
  print(i)
  gr <- gr1[i]
  sbp <- Rsamtools::ScanBamParam(mapqFilter = 255, which = gr, what = c("flag", "mapq", "seq"))
  galp0 <- GenomicAlignments::readGAlignments(Rsamtools::BamFile(bamfile), use.names = TRUE, param = sbp)
  galp0 <- galp0[mcols(galp0)$flag %in% c(0, 16)]
  galp0 <- galp0[mcols(galp0)$mapq == 255]
  
  ReferenceSpace <- GenomicAlignments::cigarRangesAlongReferenceSpace(GenomicAlignments::cigar(galp0), ops = "M", with.ops = TRUE, pos = GenomicAlignments::start(galp0))
  seqs <- S4Vectors ::mcols(galp0)$seq
  QuerySpace <- GenomicAlignments::cigarRangesAlongQuerySpace(GenomicAlignments::cigar(galp0), ops = "M", with.ops = TRUE)
  
  # Get base sequence of gr ----
  parallel::mclapply(seq_along(seqs), FUN = function(x){
    if(length(QuerySpace[[x]]) > 0) {
      XVector::subseq(seqs[rep(x, length(QuerySpace[[x]]))], QuerySpace[[x]])
    } else {
      NULL
    }
  }, mc.cores = 10) -> seqi
  
  lapply(seq_along(ReferenceSpace), function(x) {
    GenomicRanges::GRanges(seqnames = 1, ranges = ReferenceSpace[[x]], seq = seqi[[x]])
  }) -> GenoSeq
  names(GenoSeq) <- S4Vectors ::mcols(galp0)$qname
  
  mcmapply(seq_along(GenoSeq), FUN =  function(i) {
    x <- GenoSeq[[i]]
    ov <- IRanges::findOverlaps(IRanges::ranges(gr), IRanges::ranges(x))
    if(length(ov) == 0) {
      res <- NA
    } else {
      si <- x[S4Vectors::subjectHits(ov)]
      st <- GenomicRanges::start(gr) - GenomicRanges::start(si) + 1
      if(st == 0) st = 1
      res <- as.character(XVector::subseq(si$seq, start = st, width = 1))
    }
    return(res)
  }, mc.cores = 10) -> bases
  bases <- bases[!is.na(bases)]
  alphabetFrequency(DNAStringSet(paste0(bases, collapse = "")))[, 1:4]
})


# DMG-2

bamfile <- file.path("/mnt/raid64/visium/analysis_ONT/merge_fq/NGS_spaceOut/EK_4_1_DMG_2_hs/outs/possorted_genome_bam.bam")
file.exists(bamfile)

DMG_2 <- lapply(X = seq_along(gr1), function(i) {
  print(i)
  gr <- gr1[i]
  sbp <- Rsamtools::ScanBamParam(mapqFilter = 255, which = gr, what = c("flag", "mapq", "seq"))
  galp0 <- GenomicAlignments::readGAlignments(Rsamtools::BamFile(bamfile), use.names = TRUE, param = sbp)
  galp0 <- galp0[mcols(galp0)$flag %in% c(0, 16)]
  galp0 <- galp0[mcols(galp0)$mapq == 255]
  
  ReferenceSpace <- GenomicAlignments::cigarRangesAlongReferenceSpace(GenomicAlignments::cigar(galp0), ops = "M", with.ops = TRUE, pos = GenomicAlignments::start(galp0))
  seqs <- S4Vectors ::mcols(galp0)$seq
  QuerySpace <- GenomicAlignments::cigarRangesAlongQuerySpace(GenomicAlignments::cigar(galp0), ops = "M", with.ops = TRUE)
  
  # Get base sequence of gr ----
  parallel::mclapply(seq_along(seqs), FUN = function(x){
    if(length(QuerySpace[[x]]) > 0) {
      XVector::subseq(seqs[rep(x, length(QuerySpace[[x]]))], QuerySpace[[x]])
    } else {
      NULL
    }
  }, mc.cores = 10) -> seqi
  
  lapply(seq_along(ReferenceSpace), function(x) {
    GenomicRanges::GRanges(seqnames = 1, ranges = ReferenceSpace[[x]], seq = seqi[[x]])
  }) -> GenoSeq
  names(GenoSeq) <- S4Vectors ::mcols(galp0)$qname
  
  mcmapply(seq_along(GenoSeq), FUN =  function(i) {
    x <- GenoSeq[[i]]
    ov <- IRanges::findOverlaps(IRanges::ranges(gr), IRanges::ranges(x))
    if(length(ov) == 0) {
      res <- NA
    } else {
      si <- x[S4Vectors::subjectHits(ov)]
      st <- GenomicRanges::start(gr) - GenomicRanges::start(si) + 1
      if(st == 0) st = 1
      res <- as.character(XVector::subseq(si$seq, start = st, width = 1))
    }
    return(res)
  }, mc.cores = 10) -> bases
  bases <- bases[!is.na(bases)]
  alphabetFrequency(DNAStringSet(paste0(bases, collapse = "")))[, 1:4]
})



# DMG-3


bamfile <- file.path("/mnt/raid64/visium/analysis_ONT/merge_fq/NGS_spaceOut/FK_6_1_DMG_3_hs/outs/possorted_genome_bam.bam")
file.exists(bamfile)

DMG_3 <- lapply(X = seq_along(gr1), function(i) {
  print(i)
  gr <- gr1[i]
  sbp <- Rsamtools::ScanBamParam(mapqFilter = 255, which = gr, what = c("flag", "mapq", "seq"))
  galp0 <- GenomicAlignments::readGAlignments(Rsamtools::BamFile(bamfile), use.names = TRUE, param = sbp)
  galp0 <- galp0[mcols(galp0)$flag %in% c(0, 16)]
  galp0 <- galp0[mcols(galp0)$mapq == 255]
  
  ReferenceSpace <- GenomicAlignments::cigarRangesAlongReferenceSpace(GenomicAlignments::cigar(galp0), ops = "M", with.ops = TRUE, pos = GenomicAlignments::start(galp0))
  seqs <- S4Vectors ::mcols(galp0)$seq
  QuerySpace <- GenomicAlignments::cigarRangesAlongQuerySpace(GenomicAlignments::cigar(galp0), ops = "M", with.ops = TRUE)
  
  # Get base sequence of gr ----
  parallel::mclapply(seq_along(seqs), FUN = function(x){
    if(length(QuerySpace[[x]]) > 0) {
      XVector::subseq(seqs[rep(x, length(QuerySpace[[x]]))], QuerySpace[[x]])
    } else {
      NULL
    }
  }, mc.cores = 10) -> seqi
  
  lapply(seq_along(ReferenceSpace), function(x) {
    GenomicRanges::GRanges(seqnames = 1, ranges = ReferenceSpace[[x]], seq = seqi[[x]])
  }) -> GenoSeq
  names(GenoSeq) <- S4Vectors ::mcols(galp0)$qname
  
  mcmapply(seq_along(GenoSeq), FUN =  function(i) {
    x <- GenoSeq[[i]]
    ov <- IRanges::findOverlaps(IRanges::ranges(gr), IRanges::ranges(x))
    if(length(ov) == 0) {
      res <- NA
    } else {
      si <- x[S4Vectors::subjectHits(ov)]
      st <- GenomicRanges::start(gr) - GenomicRanges::start(si) + 1
      if(st == 0) st = 1
      res <- as.character(XVector::subseq(si$seq, start = st, width = 1))
    }
    return(res)
  }, mc.cores = 10) -> bases
  bases <- bases[!is.na(bases)]
  alphabetFrequency(DNAStringSet(paste0(bases, collapse = "")))[, 1:4]
})



# GBM-6


bamfile <- file.path("/mnt/raid64/visium/analysis_ONT/merge_fq/NGS_spaceOut/HW3_GBM_6_hs/outs/possorted_genome_bam.bam")
file.exists(bamfile)

GBM_6 <- lapply(X = seq_along(gr1), function(i) {
  print(i)
  gr <- gr1[i]
  sbp <- Rsamtools::ScanBamParam(mapqFilter = 255, which = gr, what = c("flag", "mapq", "seq"))
  galp0 <- GenomicAlignments::readGAlignments(Rsamtools::BamFile(bamfile), use.names = TRUE, param = sbp)
  galp0 <- galp0[mcols(galp0)$flag %in% c(0, 16)]
  galp0 <- galp0[mcols(galp0)$mapq == 255]
  
  ReferenceSpace <- GenomicAlignments::cigarRangesAlongReferenceSpace(GenomicAlignments::cigar(galp0), ops = "M", with.ops = TRUE, pos = GenomicAlignments::start(galp0))
  seqs <- S4Vectors ::mcols(galp0)$seq
  QuerySpace <- GenomicAlignments::cigarRangesAlongQuerySpace(GenomicAlignments::cigar(galp0), ops = "M", with.ops = TRUE)
  
  # Get base sequence of gr ----
  parallel::mclapply(seq_along(seqs), FUN = function(x){
    if(length(QuerySpace[[x]]) > 0) {
      XVector::subseq(seqs[rep(x, length(QuerySpace[[x]]))], QuerySpace[[x]])
    } else {
      NULL
    }
  }, mc.cores = 10) -> seqi
  
  lapply(seq_along(ReferenceSpace), function(x) {
    GenomicRanges::GRanges(seqnames = 1, ranges = ReferenceSpace[[x]], seq = seqi[[x]])
  }) -> GenoSeq
  names(GenoSeq) <- S4Vectors ::mcols(galp0)$qname
  
  mcmapply(seq_along(GenoSeq), FUN =  function(i) {
    x <- GenoSeq[[i]]
    ov <- IRanges::findOverlaps(IRanges::ranges(gr), IRanges::ranges(x))
    if(length(ov) == 0) {
      res <- NA
    } else {
      si <- x[S4Vectors::subjectHits(ov)]
      st <- GenomicRanges::start(gr) - GenomicRanges::start(si) + 1
      if(st == 0) st = 1
      res <- as.character(XVector::subseq(si$seq, start = st, width = 1))
    }
    return(res)
  }, mc.cores = 10) -> bases
  bases <- bases[!is.na(bases)]
  alphabetFrequency(DNAStringSet(paste0(bases, collapse = "")))[, 1:4]
})



GBM_6 <- data.table(ID = names(gr1), do.call(rbind, GBM_6))
colnames(GBM_6)[2:5] <- paste0("GBM_6:", colnames(GBM_6)[2:5])

DMG_3 <- data.table(ID = names(gr1), do.call(rbind, DMG_3))
colnames(DMG_3)[2:5] <- paste0("DMG_3:", colnames(DMG_3)[2:5])

DMG_2 <- data.table(ID = names(gr1), do.call(rbind, DMG_2))
colnames(DMG_2)[2:5] <- paste0("DMG_2:", colnames(DMG_2)[2:5])

DMG_1 <- data.table(ID = names(gr1), do.call(rbind, DMG_1))
colnames(DMG_1)[2:5] <- paste0("DMG_1:", colnames(DMG_1)[2:5])

NGS <- merge(merge(merge(DMG_1, DMG_2), DMG_3), GBM_6)

openxlsx::write.xlsx(NGS, "/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/analysis/14.TumorMutation/Tumor_specific_mutation_in_NGS.xlsx")




