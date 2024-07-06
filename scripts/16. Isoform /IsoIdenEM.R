#!/usr/bin/env Rscript
#author:Chao Tang
#modify:Junwei Song(2022.6.4)

options(warn=-1)
if (!"optparse" %in% installed.packages()){
  stop('There is no package called "optparse"', call. = FALSE)
}

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(S4Vectors))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(GenomicRanges))


arguments <- parse_args(
  OptionParser(usage = "%prog [options] bam", 
               description="New transcript reconstruction from noisy reads.", 
               option_list = list(
                 make_option(opt_str = c("-r","--reference"), 
                             type    = "character", 
                             default = "NULL", 
                             help    = "The prefix and name for reference genome fasta file [default %default]"),
                 make_option(opt_str = c("-g","--gtf"), 
                             type    = "character", 
                             default = "NULL", 
                             help    = "The prefix and name for gtf/gff3 file [default %default]"), 
                 make_option(opt_str = c("-o", "--output"), 
                             type    = "character", 
                             default = NULL, 
                             help    = "The prefix and name of new transcripts fasta file [default prefix of genomeBam]"), 
                 make_option(opt_str = c("-Q", "--MinMapQ"), 
                             type    = "integer", 
                             default = 60, 
                             help    = "The minimum mapq for bam file [default %default]"), 
                 make_option(opt_str = c("-s", "--ignoreStrand"), 
                             type    = "logical", 
                             default = FALSE, 
                             help    = "whether for strand-specific analysis [default %default]"),
                 make_option(opt_str = c("--MinReads1"), 
                             type    = "integer", 
                             default = 3, 
                             help    = "The minimum reads for NewJunction [default %default]"), 
                 make_option(opt_str = c("--MinReads2"), 
                             type    = "integer", 
                             default = 10, 
                             help    = "The minimum reads for NewSpliceSite [default %default]"), 
                 make_option(opt_str = c("--MinReads3"), 
                             type    = "integer", 
                             default = 10, 
                             help    = "The minimum reads for Truncated [default %default]"), 
                 make_option(opt_str = c("--MinReads4"), 
                             type    = "integer", 
                             default = 10, 
                             help    = "The minimum reads for NewSplicedGene [default %default]"), 
                 make_option(opt_str = c("--MinReads5"), 
                             type    = "integer", 
                             default = 10, 
                             help    = "The minimum reads for unspliced new transcript [default %default]"), 
                 make_option(opt_str = c("-L", "--MinL"), 
                             type    = "integer", 
                             default = 10, 
                             help    = "The minimum length of insertion or deletion for transcriptomeBam [default %default]"), 
                 make_option(opt_str = c("-c", "--cores"), 
                             type    = "integer", 
                             default = 10, 
                             help    = "The minimum cores for parallel [default %default]"))),
  positional_arguments = 1
)

# gBam <- "/mnt/raid61/Personal_data/songjunwei/Project/Map4k4_bkONT/target_pcr/AA_clean.bam"
# output <- "/mnt/raid61/Personal_data/songjunwei/Project/Map4k4_bkONT/target_pcr/Isoform_make/AA_"
# reference <- "/mnt/raid61/Personal_data/songjunwei/reference/mouse/GRCm38/release101/Mus_musculus.GRCm38.dna.primary_assembly.fa"
# gtf <- "/mnt/raid61/Personal_data/songjunwei/reference/mouse/GRCm38/release101/Mus_musculus.GRCm38.101.sorted.gtf"
# MinMapQ <- 60
ignoreStrand <- F
MR_NewJunction <- 3
MR_NewSpliceSite <- 5
MR_Truncated <- 5
MR_NewSplicedGene <- 5
MR_Unspliced <- 5
cores <- 5

gBam <- arguments$args[1]
arguments <- arguments$opt
reference <- arguments$reference
gtf <- arguments$gtf
MinMapQ <- arguments$MinMapQ
ignoreStrand <- arguments$ignoreStrand
MR_NewJunction <- arguments$MinReads1
MR_NewSpliceSite <- arguments$MinReads2
MR_Truncated <- arguments$MinReads3
MR_NewSplicedGene <- arguments$MinReads4
MR_Unspliced <- arguments$MinReads5
cores <- arguments$cores
output <- arguments$output

if(is.null(output)) {
  output <- gsub(".bam$", "", gBam)
}

logFile <- paste0(output, "NewTtanscriptReconstructionFromGenomeMapping.log")

zz <- file(logFile, open = "wt") # "w" or "wt": Open for writing in text mode.
sink(zz, type = c("message"), append = TRUE)



#### Functions #####
DonorSiteSeq <- function (x, Genome, exon = 3, intron = 20) 
{
  if (!is(x, "GRanges")) {
    stop("x must be a GRanges object of splice junction")
  }
  if (!any(is.element(is(Genome), c("TwoBitFile", "FaFile", 
                                    "BSgenome", "FaFileList", "XStringSet")))) {
    stop(paste0("unable to find an inherited method for function 'getSeq' for signature '", 
                is(Genome)[1], "'"))
  }
  donorSite <- GenomicRanges::flank(x = x, width = max(exon, 
                                                       intron), both = T, start = TRUE)
  donorMotif <- BSgenome::getSeq(Genome, donorSite)
  if (exon == intron) {
    Seq <- donorMotif
  }
  else {
    if (exon < intron) {
      Seq <- subseq(donorMotif, start = (intron - exon + 
                                           1))
    }
    else {
      Seq <- subseq(donorMotif, end = (exon + intron))
    }
  }
  x$DonorMotif <- Seq
  return(x)
}


AcceptorSiteSeq <- function (x, Genome, intron = 20, exon = 3) 
{
  if (!is(x, "GRanges")) {
    stop("x must be a GRanges object of splice junction")
  }
  if (!any(is.element(is(Genome), c("TwoBitFile", "FaFile", 
                                    "BSgenome", "FaFileList", "XStringSet")))) {
    stop(paste0("unable to find an inherited method for function 'getSeq' for signature '", 
                is(Genome)[1], "'"))
  }
  acceptorSite <- GenomicRanges::flank(x = x, width = max(exon, 
                                                          intron), both = T, start = FALSE)
  acceptorMotif <- BSgenome::getSeq(Genome, acceptorSite)
  if (exon == intron) {
    Seq <- acceptorMotif
  }
  else {
    if (exon < intron) {
      Seq <- subseq(acceptorMotif, end = (intron + exon))
    }
    else {
      Seq <- subseq(acceptorMotif, start = (exon - intron + 
                                              1))
    }
  }
  x$AcceptorMotif <- Seq
  return(x)
}


Clu2L <- function(x) {
  pi <- unlist(base::strsplit(x, "[:;-]"))
  paste(pi[1], pi[-1], sep = ":")
}

Clu2SJ <- function(x) {
  pi <- unlist(base::strsplit(x, "[:;]"))
  paste(pi[1], pi[-1], sep = ":")
}

SpliceMotif <- function(x, Genome) {
  if(!is(x, "GRanges")) {
    g <- as(x, "GRanges")
  } else {
    g <- x
  }
  g <- DonorSiteSeq(x = g, Genome = genome, exon = 0, intron = 2)
  g <- AcceptorSiteSeq(x = g, Genome = genome, exon = 0, intron = 2)
  paste0(g$DonorMotif, g$AcceptorMotif) %in% c("GTAG", "CTAC", "GCAG", "CTGC", "ATAC", "GTAT")
}

CanonicalSite <- function(x, Genome) {
  if(!is(x, "GRanges")) {
    g <- as(x, "GRanges")
  } else {
    g <- x
  }
  g <- DonorSiteSeq(x = g, Genome = genome, exon = 0, intron = 2)
  g <- AcceptorSiteSeq(x = g, Genome = genome, exon = 0, intron = 2)
  paste0(g$DonorMotif, g$AcceptorMotif) %in% c("GTAG", "GCAG", "ATAC")
}

bam2GRange <- function(map, sj, genome) {
  ft <- min(start(map)):max(end(map))
  sjsi <- as(sj, "GRanges")
  gap <- do.call(c, lapply(seq_along(sjsi), function(y) start(sjsi[y]):end(sjsi[y])))
  exongr <- GRanges(seqnames = as.character(runValue(seqnames(map))), 
                    ranges = IRanges::shift(IRanges(!ft %in% gap), shift = min(ft) - 1), 
                    strand = as.character(runValue(strand(sjsi))))
  exongr$type <- "exon"
  txgr <- range(exongr)
  txgr$type <- "transcript"
  res <- c(txgr, exongr)
  
  if(all(strand(as(sj, "GRanges")) == "*")) {
    strand(res) <- ifelse(sum(CanonicalSite(x = paste0(sj, ":+"), Genome = genome)) >= sum(CanonicalSite(x = paste0(sj, ":-"), Genome = genome)), "+", "-")
  }
  return(res)
}


#######  Host and EM ####### 
HostGeneIso <- function(from_vec,gtf,other=NULL,FromName="transcript_name",ToName="gene_name",check2=F) {
  gtf.dataframe <- as.data.frame(gtf)
  covert <- unique(data.frame(from=gtf.dataframe[,FromName],to=gtf.dataframe[,ToName]))
  covert <- na.omit(covert)
  covert$from <- as.character(covert$from)
  covert$to <- as.character(covert$to)

  #anno
  to_vec <- plyr::mapvalues(from = covert$from, to = covert$to,
                            x = as.character(from_vec), warn_missing = F)
  
  #novel
  if (length(grep("[0-9A-Za-z]:[0-9]",to_vec))!=0) {
    #helped by SJs
    if (!is.null(other)) {
      #sj info
      spl <- stringr::str_split_fixed(rownames(other),"@",3)
      tmp1 <- data.table(gene=spl[,1],iso=spl[,2],sj=spl[,3])
      tmp2 <- tmp1[grep("[0-9A-Za-z]:[0-9]",tmp1$iso)]
      novel_iso <- data.frame(iso=grep("[0-9A-Za-z]:[0-9]",to_vec,value = T))
      to_vec <- plyr::mapvalues(from = tmp2$iso, to = tmp2$gene,
                                x = novel_iso$iso, warn_missing = F)
    }
    #findoverlap
    novel_iso <- data.frame(iso=grep("[0-9A-Za-z]:[0-9]",to_vec,value = T))
    tmp3 <- stringr::str_split_fixed(novel_iso$iso,":",3)
    tmp4 <- ifelse(tmp3[,3]=="+"| tmp3[,3]=="-", 
                   paste0(tmp3[,1],":",tmp3[,2],":",tmp3[,3]) ,paste0(tmp3[,1],":",tmp3[,2],":*")) 
    tmp4 <- data.frame(sj=tmp4,iso=novel_iso)
    
    # tmp5 <- stringr::str_split_fixed(tmp3[,2],"-",2)
    # table(tmp5[,1]<tmp5[,2])
    
    qry <- tmp4$sj %>% as(., "GRanges")
    qry$iso <- tmp4$iso
    qry$sj <-  paste0(qry@seqnames, ":", qry@ranges, ":", qry@strand)
    names(qry) <- qry$sj
    gtftmp <- gtf.dataframe[, c("seqnames", "start", "end", "strand",ToName)]
    colnames(gtftmp) <- c("seqnames", "start", "end", "strand","Toname")
    obj <- makeGRangesFromDataFrame(gtftmp, keep.extra.columns = TRUE)
    names(obj) <- obj$Toname
    hts = findOverlaps(qry, obj)
    find <- unique(data.frame(sj = names(qry)[queryHits(hts)], gene = names(obj)[subjectHits(hts)]))
    Hostgene <- merge(find,mcols(qry))
    tmp5 <- as.data.table(Hostgene)[,.N,by=c("sj")]
    hostgene <- Hostgene[Hostgene$sj %in% tmp5[tmp5$N==1]$sj,]
    to_vec <- plyr::mapvalues(from = hostgene$iso, to = hostgene$gene,
                              x = to_vec, warn_missing = F) 
  }
  #novel double check
  if (length(grep("[0-9A-Za-z]:[0-9]",to_vec))!=0 & check2==T) {
    gtf_gene = gtf.dataframe[gtf.dataframe$type=="gene",][,c("seqnames","start","end","gene_id","gene_name")]
    gtf_gene$seqnames=gsub("chr","",gtf_gene$seqnames)
    
    check_gene_belong <- function(x,gtf_gene){
      tmp = unlist(strsplit(x, ":"))
      tmp_start=unlist(strsplit(tmp[2], "-"))[1]
      tmp_end=unlist(strsplit(tmp[length(tmp)], "-"))[2]
      chr=gsub("chr","",tmp[1])
      tmp3=as.numeric(c(chr,tmp_start,tmp_end))
      gtf_gene_sel_chr = gtf_gene[gtf_gene$seqnames == tmp3[1], ]
      start = unique(gtf_gene_sel_chr$start)
      end = unique(gtf_gene_sel_chr$end)
      my_end = head(end[which(end > tmp3[2])], n = 1)
      my_start = tail(start[which(start < tmp3[2])], n = 1)
      gene_name = gtf_gene_sel_chr[which(gtf_gene_sel_chr$start == my_start | gtf_gene_sel_chr$end == my_end), "gene_name"]
      gene_name = ifelse(gene_name[1]=="","_",paste(gene_name,collapse = "_"))
      return(gene_name)
    }
    novel_iso <- data.frame(iso=grep("[0-9A-Za-z]:[0-9]",to_vec,value = T))
    novel_iso$gene <- unlist(lapply(novel_iso$iso, function(x){   check_gene_belong(x,gtf_gene)   }))
    to_vec <- plyr::mapvalues(from = novel_iso$iso, to = novel_iso$gene,
                              x = to_vec, warn_missing = F)
  }
  return(to_vec)
}

filter_transcripts <- function(data,cores=5) {
  data$suffix <- as.numeric(sub(".*-", "", data$transcript_name))
  data$suffix[grep(":",data$transcript_name)] <- "NA"
  # 计算每个transcript_name的计数, 找到每个read_id的最大count
  counts <- merge(data,data.table(data)[,.(count=.N),by=transcript_name],by="transcript_name",sort=F)
  max_counts <- aggregate(count ~ read_id, counts, max)
  max_data <- merge(counts, max_counts, by = c("read_id", "count"),sort=F)
  # 对于每个Read ID，如果有多个transcript_name具有相同的最大计数，选择后缀最小的
  result_data <- do.call(rbind, mclapply(split(max_data, max_data$read_id), function(df) {
    if (nrow(df) > 1) {
      df[df$suffix == min(df$suffix), ]
    } else {
      df
    }
  },mc.cores = cores))
  final_data <- data.table(result_data[, c("read_id", "transcript_id", "transcript_name", "gene_name")])
  rownames(final_data) <- NULL
  return(final_data)
}

convert_to_data_list <- function(isoTab) {
  setDT(isoTab)
  unique_transcripts <- isoTab[, .(ref = list(unique(transcript_name))), by = read_id]
  data_list <- setNames(lapply(unique_transcripts$ref, function(x) list(ref = x)), unique_transcripts$read_id)
  return(data_list)
}

get_compatibility <- function(data_list) {
  compatibility_dict <- list()
  
  for (k in 1:length(data_list)) {
    read_name <- names(data_list[k])
    
    for (alignment in data_list[[k]]$ref) {
      score <- 1.0 / length(data_list[[k]]$ref)
      compatibility_dict[paste(read_name, alignment, sep = "_")] <- score
    }
  }
  return(compatibility_dict)
}


calculate_abundance <- function(abundance_dict,compatibility_dict) {
  last_abundance_dict <- abundance_dict
  abundance_dict <- list()
  total <- 0
  tmp_convergence <- 0
  abundance_dict <- lapply(abundance_dict, as.numeric)
  
  for (pair_name in names(compatibility_dict)) {
    #拆分
    split_string <- strsplit(pair_name, "_")[[1]]
    ref_name <- split_string[length(split_string)]
    #计算ref_score
    score <- as.numeric(unlist(compatibility_dict[[pair_name]]))
    if (is.null(abundance_dict[[ref_name]])) {
      abundance_dict[[ref_name]] <- 0
    }
    abundance_dict[[ref_name]] <- abundance_dict[[ref_name]] + score
    total <- total + score
  }
  
  for (ref_name in names(abundance_dict)) {
    #计算丰度
    abundance_dict[[ref_name]] <- as.numeric(unlist(abundance_dict[[ref_name]])) / total
    
    if (em_round > 1) {
      tmp_convergence <- tmp_convergence + abs(last_abundance_dict[[ref_name]] - abundance_dict[[ref_name]])
    }
  }
  if (em_round == 1) {
    convergence <- 1
  } else {
    convergence <<- tmp_convergence  #global
  }
  return(abundance_dict)
}

update_compatibility <- function(abundance_dict,compatibility_dict) {
  total <- list()
  abundance_dict <- lapply(abundance_dict, as.numeric)
  total <- lapply(total, as.numeric)
  
  for (pair_name in names(compatibility_dict)) {
    #拆分
    split_string <- strsplit(pair_name, "_")[[1]]
    ref_name <- split_string[length(split_string)]
    read_name <- split_string[1]
    if (is.null(total[[read_name]])) {
      total[[read_name]] <- 0
    }
    total[[read_name]] <- total[[read_name]] + abundance_dict[[ref_name]]
  }
  
  for (pair_name in names(compatibility_dict)) {
    #拆分
    split_string <- strsplit(pair_name, "_")[[1]]
    ref_name <- split_string[length(split_string)]
    read_name <- split_string[1]
    
    compatibility_dict[[pair_name]] <- abundance_dict[[ref_name]] / total[[read_name]]
  }
  
  return(compatibility_dict)
}

## EM main Function
EM_calculate <- function(input_data,gtf,min_gene_reads=5,cores=5){

  input_dataf <- filter_transcripts(data = input_data,cores = cores)
  input_dataff <- input_dataf[gene_name %in% input_dataf[,.N,by="gene_name"][N>=min_gene_reads]$gene_name]
  
  output_df <- mclapply(split(input_dataff,f = input_dataff$gene_name), function(x){
    data_list <- convert_to_data_list(x) # convert data
    compatibility_dict <- get_compatibility(data_list) #直接计算每条reads的ref占比
    
    em_round <- 0
    convergence <- 1 #初始值
    convergence_target <- 0.005 #变化值
    max_em_rounds <- 300
    abundance_dict <- list()
    last_abundance_dict <- list()
    # Iterate until convergence threshold or max EM round are reached
    while (convergence > convergence_target && em_round < max_em_rounds) {
      em_round <- em_round + 1
      # Calculate abundance from compatibility assignments    计算丰度
      abundance_dict <- calculate_abundance(abundance_dict,compatibility_dict)
      # Update compatibility assignments  更新兼容性字典
      compatibility_dict <- update_compatibility(abundance_dict,compatibility_dict)
      cat("EM Round:", em_round, "/ Convergence value:", convergence ,"\n")
    }
    print(paste("Exit EM loop after", em_round, "rounds"))
    print(paste("Convergence value:", convergence))
    if (!(convergence <= convergence_target)) {
      print(paste("Convergence target (", convergence_target, ") could not be reached after", max_em_rounds, "rounds"))
    }
    
    # Compute estimated counts and TPM
    count_df <- data.frame(transcript_name = names(abundance_dict), em_ratio = unlist(abundance_dict), stringsAsFactors = FALSE)
    rownames(count_df) <- count_df$transcript_name
    count_df$est_count <- count_df$em_ratio * length(data_list)
    return(count_df)
  },mc.cores = cores) %>% rbindlist()
  return(output_df)
}




# Read genome ----
message("Read genome")

genome <- readDNAStringSet(reference)
names(genome) <- mapply(function(x) x[1], base::strsplit(names(genome), " "))


# Phase gtf ----
message("Phase gtf")

TxDb <- suppressMessages(suppressWarnings(makeTxDbFromGFF(gtf)))

GeneRegion <- genes(TxDb)

Tx2Gene <- data.table(Tx = mcols(unlist(transcriptsBy(TxDb, "gene")))$tx_name, Gene = names(unlist(transcriptsBy(TxDb, "gene"))))
AllIntrons <- unlist(intronsByTranscript(TxDb, use.names = T))
AllIntrons <- data.table(Tx = names(AllIntrons), introns = paste0(seqnames(AllIntrons), ":", start(AllIntrons), "-", end(AllIntrons)), strand = as.character(strand(AllIntrons)))
AllIntrons <- merge(Tx2Gene, AllIntrons, by = "Tx")

site4Gene <- unlist(intronsByTranscript(TxDb, use.names = TRUE))
site4Gene$Tx <- names(site4Gene)
site4Gene <- as.data.table(site4Gene)
site4Gene <- merge(Tx2Gene, site4Gene, by = "Tx")

if(ignoreStrand) {
  site4Gene[, start := paste0(seqnames, ":", start)]
  site4Gene[, end := paste0(seqnames, ":", end)]
} else {
  site4Gene[, start := paste0(seqnames, ":", start, ":", strand)]
  site4Gene[, end := paste0(seqnames, ":", end, ":", strand)]
}

site4Gene <- site4Gene[, .(seqnames = unique(seqnames), strand = unique(strand), Site = list(c(start, end))), by = Gene]


sj4Tx <- unlist(intronsByTranscript(TxDb, use.names = TRUE))
sj4Tx$Tx <- names(sj4Tx)
sj4Tx <- as.data.table(sj4Tx)
sj4Tx[, range := paste(start, end, sep = "-")]
sj4Tx[, seqnames := as.character(seqnames)]
sj4Tx[, intron := paste0(seqnames, ":", range)]

setkey(sj4Tx, Tx, start, end)
sj4Tx <- split(sj4Tx, sj4Tx$seqnames)

sj4Tx <- mclapply(sj4Tx, FUN = function(x) {
  x[, .SD[, .(strand = unique(strand), SJ = paste(unique(seqnames), paste(range, collapse = ":"), sep = ":"), introns = list(intron))], by = Tx]
}, mc.cores = cores)
sj4Tx <- do.call(rbind, sj4Tx)


# Read bam ----
message("Read bam")

params <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isNotPassingQualityControls = FALSE), 
                                  what = c("mapq", "flag"), 
                                  mapqFilter = MinMapQ)
map0 <- GenomicAlignments::readGAlignments(file = gBam, param = params, use.names = TRUE)
map0 <- map0[S4Vectors::mcols(map0)$flag %in% c(0, 16)] # remove secondary and supplementary alignment
map1 <- map0[njunc(map0) > 0]
map0 <- map0[njunc(map0) == 0]
gc()


# Junction of reads ----
message("Junction of reads")

read2sj <- unlist(junctions(map1))
read2sj$Read <- names(read2sj)
read2sj <- as.data.table(read2sj)
read2sj[, range := paste(start, end, sep = "-")]

read2sj[, seqnames := as.character(seqnames)]
setkey(read2sj, Read, start, end)
read2sj <- split(read2sj, read2sj$seqnames)

read2sj <- mclapply(read2sj, FUN = function(x) {
  x[, .SD[, .(SJ = paste(unique(seqnames), paste(range, collapse = ":"), sep = ":"), start = min(start), end = max(end), strand = unique(strand))], by = Read]
}, mc.cores = cores)

read2sj <- cbind(do.call(rbind, read2sj), seqnames = rep(names(read2sj), mapply(nrow, read2sj)))

if(ignoreStrand) {
  sj2read <- read2sj[, .(Reads = list(Read)), by = list(SJ, seqnames, start, end)]
} else {
  sj2read <- read2sj[, .(Reads = list(Read)), by = list(SJ, seqnames, start, end, strand)]
}


# known ----
message("known transcripts")

## known spliced transcripts ----

if(ignoreStrand) {
  Equal2Tx <- merge(sj4Tx, sj2read, by = "SJ")
} else {
  Equal2Tx <- merge(sj4Tx, sj2read, by = c("SJ", "strand"))
}

## known unspliced transcripts ----

UnsplicedTx <- intronsByTranscript(TxDb, use.names = TRUE)
UnsplicedTx <- transcripts(TxDb, use.names = TRUE)[setdiff(names(UnsplicedTx), names(unlist(UnsplicedTx)))]

map02Tx <- as.data.table(GenomicRanges::findOverlaps(GRanges(map0), UnsplicedTx, ignore.strand = ignoreStrand))
map02Tx_uniq <- map02Tx[queryHits %in% map02Tx[, .N, queryHits][N == 1, queryHits]]
read2Tx_uniq <- data.table(Read = names(map0)[map02Tx_uniq[, queryHits]], Tx = names(UnsplicedTx)[map02Tx_uniq[, subjectHits]])

map02Tx2 <- as.data.table(GenomicRanges::findOverlaps(GRanges(map0), UnsplicedTx, type = "within", ignore.strand = ignoreStrand))
map02Tx_uniq <- map02Tx2[queryHits %in% map02Tx2[, .N, queryHits][N == 1, queryHits]]
read2Tx_uniq2 <- data.table(Read = names(map0)[map02Tx_uniq[, queryHits]], Tx = names(UnsplicedTx)[map02Tx_uniq[, subjectHits]])

read2Tx_uniq <- unique(rbind(read2Tx_uniq, read2Tx_uniq2))
read2Tx_uniq <- read2Tx_uniq[, .(Reads = list(Read)), Tx]

Reads2KnownTx <- c(Equal2Tx$Reads, read2Tx_uniq$Reads)
names(Reads2KnownTx) <- c(Equal2Tx$Tx, read2Tx_uniq$Tx)


# Novel ----
message("Novel transcript")

## Novel spliced transcript ----

if(ignoreStrand) {
  sj2read_novel <- sj2read[!SJ %in% sj4Tx$SJ]
  sj2read_novel <- sj2read_novel[, .(Reads = list(unlist(Reads))), by = list(SJ, seqnames, start, end)]
  sj2read_novel$N <- mapply(length, sj2read_novel$Reads)
} else {
  sj2read_novel <- sj2read[!sj2read[, paste0(SJ, ":", strand)] %in% sj4Tx[, paste0(SJ, ":", strand)]]
  sj2read_novel$N <- mapply(length, sj2read_novel$Reads)
  sj2read_novel <- sj2read_novel[, .SD[which.max(N),], SJ]
}

sj2read_novel <- sj2read_novel[N >= min(MR_NewJunction, MR_NewSpliceSite, MR_Truncated, MR_NewSplicedGene)]


if(ignoreStrand) {
  Int4Tx_List <- split(AllIntrons[, introns], AllIntrons[, Tx])
} else {
  Int4Tx_List <- split(AllIntrons[, paste0(introns, ":", strand)], AllIntrons[, Tx])
}

site4Gene_List <- site4Gene$Site
names(site4Gene_List) <- site4Gene$Gene



NewSplicedTxGRanges <- mclapply(seq_len(nrow(sj2read_novel)), function(i) {
  if(ignoreStrand) {
    sj <- sj2read_novel[i, Clu2SJ(SJ)]
  } else {
    sj <- paste0(sj2read_novel[i, Clu2SJ(SJ)], ":", sj2read_novel[i, strand])
  }
  
  if(ignoreStrand) {
    ss <- sj2read_novel[i, Clu2L(SJ)]
  } else {
    ss <- paste0(sj2read_novel[i, Clu2L(SJ)], ":", sj2read_novel[i, strand])
  }
  
  isTruncatedTx <- mapply(function(x) all(sj %in% x), Int4Tx_List)
  
  Truncated <- NULL
  if(sj2read_novel[i, N] >= MR_Truncated) {
    if(any(isTruncatedTx)) {
      mapi <- map1[unlist(sj2read_novel[i, Reads])]
      res <- bam2GRange(map = mapi, sj = sj, genome = genome)
      res$gene_id <- AllIntrons[Tx %in% names(which(isTruncatedTx)), paste0(unique(Gene), collapse = "|")]
      res$transcript_id <- sj2read_novel[i, SJ]
      res$source <- "Truncated"
      res$NReads <- sj2read_novel[i, N]
      Truncated <- res
      if(ignoreStrand) {
        strand(Truncated) <- ifelse(sum(CanonicalSite(x = paste0(sj, ":+"), Genome = genome)) >= sum(CanonicalSite(x = paste0(sj, ":-"), Genome = genome)), "+", "-")
      }
    } else {
      Truncated <- NULL
    }
  } else {
    Truncated <- NULL
  }
  
  NewJunction <- NULL
  if(is.null(Truncated)) {
    if(sj2read_novel[i, N] >= MR_NewJunction) {
      if(any(isTruncatedTx)) {
        NewJunction <- NULL
      } else {
        NewJunction <- mapply(function(x) all(ss %in% x), site4Gene_List)
        if(any(NewJunction)) {
          mapi <- map1[unlist(sj2read_novel[i, Reads])]
          res <- bam2GRange(map = mapi, sj = sj, genome = genome)
          res$gene_id <- paste0(unique(names(which(NewJunction))), collapse = "|")
          res$transcript_id <- sj2read_novel[i, SJ]
          res$source <- "NewJunction"
          res$NReads <- sj2read_novel[i, N]
          NewJunction <- res
        } else {
          NewJunction <- NULL
        }
      }
    } else {
      NewJunction <- NULL
    }
  }
  
  NewSpliceSite <- NULL
  if(is.null(NewJunction)) {
    if(sj2read_novel[i, N] >= MR_NewSpliceSite) {
      if(any(isTruncatedTx)) {
        NewSpliceSite <- NULL
      } else {
        if(any(ss %in% unlist(site4Gene_List, use.names = FALSE)) & any(!ss %in% unlist(site4Gene_List, use.names = FALSE))) {
          NewSpliceSite <- mapply(function(x) any(ss %in% x), site4Gene_List)
          mapi <- map1[unlist(sj2read_novel[i, Reads])]
          res <- bam2GRange(map = mapi, sj = sj, genome = genome)
          res$gene_id <- paste0(unique(names(which(NewSpliceSite))), collapse = "|")
          res$transcript_id <- sj2read_novel[i, SJ]
          res$source <- "NewSpliceSite"
          res$NReads <- sj2read_novel[i, N]
          NewSpliceSite <- res
        } else {
          NewSpliceSite <- NULL
        }
      }
    } else {
      NewSpliceSite <- NULL
    }
  }
  
  NewSplicedGene <- NULL
  if(is.null(NewSpliceSite)) {
    if(sj2read_novel[i, N] >= MR_NewSplicedGene) {
      if(any(isTruncatedTx)) {
        NewSplicedGene <- NULL
      } else {
        if(any(ss %in% unlist(site4Gene_List, use.names = FALSE))) {
          NewSplicedGene <- NULL
        } else {
          mapi <- map1[unlist(sj2read_novel[i, Reads])]
          res <- bam2GRange(map = mapi, sj = sj, genome = genome)
          res$gene_id <- sj2read_novel[i, SJ]
          res$transcript_id <- paste0("transcript:", sj2read_novel[i, SJ])
          res$source <- "NewSplicedGene"
          res$NReads <- sj2read_novel[i, N]
          NewSplicedGene <- res
        }
      }
    } else {
      NewSplicedGene <- NULL
    }
  }
  
  res <- c(Truncated, NewJunction, NewSpliceSite, NewSplicedGene)
  if(is.null(res)) {
    res <- NULL
  } else {
    if(is.list(res)) {
      res <- res[[1]]
    } else {
      res
    }
  }
  return(res)
}, mc.cores = cores)

NewSplicedTxGRanges <- NewSplicedTxGRanges[!mapply(is.null, NewSplicedTxGRanges)]



## Novel unspliced transcript  ----

if(ignoreStrand) {
  NewUnsplicedTxGRanges <- NULL
} else {
  map0 <- map0[countOverlaps(GRanges(map0), UnsplicedTx) == 0]
  ebg <- unlist(exonsBy(TxDb, "gene"))
  map0 <- map0[countOverlaps(GRanges(map0), ebg, type = "within") == 0] # bu shu yu ren he yi ge exon
  map0Cov <- coverage(GRanges(map0))
  
  map0Cov <- lapply(map0Cov, function(x) IRanges(x >= MR_Unspliced))
  map0Cov <- map0Cov[mapply(length, map0Cov) > 0]
  map0Cov <- do.call(c, lapply(seq_along(map0Cov), function(x) GRanges(seqnames = names(map0Cov)[x], ranges = map0Cov[[x]])))
  
  cluster <- split(names(map0)[subjectHits(findOverlaps(map0Cov, GRanges(map0)))], queryHits(findOverlaps(map0Cov, GRanges(map0))))
  
  mclapply(cluster, function(x) {
    i <- map0[x]
    i <- i[width(i)/width(reduce(ranges(i))) > 0.8]
    if(length(i) > 1) {
      i <- i[strand(i) == names(which.max(table(strand(i))))]
    } else {
      return(NULL)
    }
    
    if(length(i) > 0) {
      i <- reduce(GRanges(i))
    } else {
      return(NULL)
    }
    res <- c(i, i)
    res$type <- c("transcript", "exon")
    res$gene_id <- as.character(i)
    res$transcript_id <- paste0("transcript:", as.character(i))
    res$source <- "NewUnsplicedGene"
    res$NReads <- length(x)
    return(res)
  }, mc.cores = cores) -> NewUnsplicedTxGRanges
  
  NewUnsplicedReads2Tx <- cluster[!mapply(is.null, NewUnsplicedTxGRanges)]
  
  NewUnsplicedTxGRanges <- NewUnsplicedTxGRanges[!mapply(is.null, NewUnsplicedTxGRanges)]
  names(NewUnsplicedTxGRanges) <- NULL
  
  names(NewUnsplicedReads2Tx) <- mapply(NewUnsplicedTxGRanges, FUN = function(x) unique(as.character(x)))
}


# New transcript sequence  ----
message("New transcript sequence")


if(!is.null(NewUnsplicedTxGRanges)) {
  NewTxGRanges <- c(NewSplicedTxGRanges, NewUnsplicedTxGRanges)
} else {
  NewTxGRanges <- NewSplicedTxGRanges
}

NewTxSeq <- mclapply(NewTxGRanges, function(x) {
  TxSeq <- DNAStringSet(do.call(c, BSgenome::getSeq(genome, x[x$type == "exon"])))
  names(TxSeq) <- unique(x$transcript_id)
  return(TxSeq)
}, mc.cores = cores)
NewTxSeq <- do.call(c, NewTxSeq)

writeXStringSet(NewTxSeq, filepath = paste0(output, "NewIsoform.fa"))  



# New transcript reads  ----
message("New transcript reads")


if(!is.null(NewUnsplicedTxGRanges)) {
  NewSplicedTxReads <- sj2read_novel[, Reads]
  names(NewSplicedTxReads) <- sj2read_novel[, SJ]
  NewSplicedTxReads <- NewSplicedTxReads[names(NewSplicedTxReads) %in% names(NewTxSeq)]
  NewTxReads <- c(NewSplicedTxReads, NewUnsplicedReads2Tx)
} else {
  NewTxReads <- sj2read_novel[, Reads]
  names(NewTxReads) <- sj2read_novel[, SJ]
  NewTxReads <- NewTxReads[names(NewTxReads) %in% names(NewTxSeq)]
}

Reads2Tx <- c(Reads2KnownTx, NewTxReads)
jsonlite::write_json(Reads2Tx, paste0(output, "Reads2transcript.json"))



# New transcript GTF   ----
message("New transcript GTF")

gtf_Tab <- fread(cmd = paste('zgrep -v "#"', gtf))
gtf_Tab <- gtf_Tab[V3 == "gene"]
gtf_Tab$ID <- gsub('\\"', "", gsub("gene_id ", "", mapply(function(x) x[1], strsplit(gtf_Tab[, V9], ";"))))

NewTxGRangesTab <- as.data.table(unlist(GRangesList(NewTxGRanges)))

NewTxKnownGeneGTF <- unique(gtf_Tab[gtf_Tab$ID %in% NewTxGRangesTab$gene_id, 1:9])

NewTxGRangesTabNewGene <- NewTxGRangesTab[!gene_id %in% gtf_Tab$ID, ]
NewTxNewGeneGTF <- NewTxGRangesTabNewGene[, .(V1 = unique(seqnames), 
                                              V2 = paste0(unique(source), collapse = "|"), 
                                              V3 = "gene", 
                                              V4 = min(start), 
                                              V5 = max(end), 
                                              V6 = ".", 
                                              V7 = unique(strand), 
                                              V8 = ".", 
                                              V9 = paste0('gene_id "', unique(gene_id), '";')), by = gene_id]

NewTxGTF <- NewTxGRangesTab[type == "transcript", .(V1 = seqnames, 
                                                    V2 = paste0(unique(source), collapse = "|"), 
                                                    V3 = "transcript", 
                                                    V4 = start, 
                                                    V5 = end, 
                                                    V6 = ".", 
                                                    V7 = strand, 
                                                    V8 = ".", 
                                                    V9 = paste0('gene_id "', gene_id, '"; transcript_id "', transcript_id, '";')), by = transcript_id]

NewExonGTF <- NewTxGRangesTab[type == "exon", .(V1 = seqnames, 
                                                V2 = paste0(unique(source), collapse = "|"), 
                                                V3 = "exon", 
                                                V4 = start, 
                                                V5 = end, 
                                                V6 = ".", 
                                                V7 = strand, 
                                                V8 = ".", 
                                                V9 = paste0('gene_id "', gene_id, '"; transcript_id "', transcript_id, '";')), by = list(transcript_id, seqnames, start, end, strand)]


GTFTab <- rbind(NewTxKnownGeneGTF[, paste0("V", 1:9)], NewTxNewGeneGTF[, paste0("V", 1:9)], NewTxGTF[, paste0("V", 1:9)], NewExonGTF[, paste0("V", 1:9)])

GTFTabByGene <- split(GTFTab, mapply(function(x) grep("gene_id", x, value = TRUE), strsplit(GTFTab$V9, ";")))

New_GTF <- do.call(rbind, mclapply(GTFTabByGene, function(x) {
  Tf <- x[V3 == "gene"]
  Ts <- x[V3 != "gene"]
  Ts$TxID <- mapply(function(x) grep("transcript_id", x, value = TRUE), strsplit(Ts$V9, ";"))
  setkey(Ts, TxID, V4)
  rbind(Tf, Ts[, 1:9])
}, mc.cores = cores))

New_GTF <- unique(New_GTF)
fwrite(New_GTF, paste0(output, "NewIsoform.gtf"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


system(command = paste0("cp ", gtf, " ", output, "New.gtf"))
fwrite(x = New_GTF[V2 %in% c("NewJunction", "NewSplicedGene", "NewSpliceSite", "NewUnsplicedGene", "Truncated")], 
       file = paste0(output, "New.gtf"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)




# EM calculate   ----
new_gtf <- rtracklayer::import(paste0(output, "New.gtf"))

#Reads2Tx <- jsonlite::read_json(paste0(output, "Reads2transcript.json"))
# iso_tab <- do.call(rbind, lapply(names(Reads2Tx), function(tx) {
#   data.frame(transcript_id = tx, read_id = Reads2Tx[[tx]])
# }))

iso_tab <- data.frame(read_id   =  unlist(Reads2Tx, use.names = F),
                  transcript_id = rep(names(Reads2Tx),sapply(Reads2Tx, length)))
iso_tab$transcript_id <- gsub("\\.[0-9]*$","",iso_tab$transcript_id)
iso_tab$transcript_name <- HostGeneIso(from_vec = iso_tab$transcript_id,gtf = new_gtf,
                                       FromName = "transcript_id",ToName = "transcript_name")
iso_tab$gene_name <- HostGeneIso(from_vec = iso_tab$transcript_id,gtf = new_gtf,
                                       FromName = "transcript_id",ToName = "gene_name")
fwrite(iso_tab,paste0(output, "Reads2transcript.txt"))


iso_mat <- EM_calculate(input_data = iso_tab,gtf = new_gtf,min_gene_reads = 5,cores = 20)
fwrite(iso_mat,paste0(output, "isoform_mat.txt"))

sink(type = "message")
close(zz)
## end
