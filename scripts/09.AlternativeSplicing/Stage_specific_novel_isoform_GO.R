library(tidyr)
library(ggplot2)
library(data.table)
library(clusterProfiler)

library(ggplot2)
library(data.table)
novel_iso <- readRDS("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/00analysis/rep4_pb_analysis/01isoform_rawtab.rds")
novel_iso
novel_iso$Stage <- mapply(function(x) tail(x, 1), strsplit(novel_iso$sample, "\\."))
known_iso <- novel_iso[grepl("^PB", isoform)]
novel_iso <- novel_iso[!grepl("^PB", isoform)]

Tx_N <- rbind(known_iso[, .(Annotation = "Annotated", N = length(unique(isoform))), Stage],
              novel_iso[, .(Annotation = "Novel", N = length(unique(isoform))), Stage])

Tx_N <- Tx_N[, .(Annotation, N, P = N / sum(N) * 100), Stage]
Tx_N[, Stage := factor(Stage, levels = c("Trophozoite", "Schizont"))]
Tx_N[, L := paste0(N, "\n(", round(P, 1), "%)")]

ggplot(Tx_N, aes(x = Stage, y = N, fill = Annotation)) + 
  geom_col() + 
  geom_text(aes(label = L), position = position_stack(vjust = 0.5)) +
  theme_bw(base_size = 15) +
  # guides(fill = guide_legend(title = "")) +
  theme(legend.position = "top", legend.title = element_blank()) +
  labs(y = "Number of isoforms")
# ggsave("/mnt/raid61/Personal_data/tangchao/Temp/20211025/05.TranscriptClean/No.Isoform.pdf", width = 3, height = 4)
ggsave("/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/analysis/09.AlternativeSplicing/01.Plasmodium/20211025/No.Isoform.pdf", width = 3, height = 4)



novel_iso <- readRDS("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/00analysis/rep4_pb_analysis/rep4_pb_6sample_isoform.rds")
novel_iso <- readRDS("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/00analysis/rep4_pb_analysis/01isoform_rawtab.rds")
novel_iso <- do.call(rbind, novel_iso)
novel_iso$Stage <- mapply(function(x) tail(x, 1), strsplit(novel_iso$sample, "\\."))
novel_iso <- novel_iso[!grepl("^PB", isoform)]

novel_iso[, length(unique(isoform)), Stage]

Schizont_508 <- setdiff(novel_iso[Stage == "Schizont", unique(isoform)], novel_iso[Stage != "Schizont", unique(isoform)])
Schizont_508_hostgene <- novel_iso[isoform %in% Schizont_508, unique(gene)]

Trophozoite_72 <- setdiff(novel_iso[Stage != "Schizont", unique(isoform)], novel_iso[Stage == "Schizont", unique(isoform)])
Trophozoite_72_hostgene <- novel_iso[isoform %in% Trophozoite_72, unique(gene)]

intersect(Schizont_508_hostgene, Trophozoite_72_hostgene)


ID_transform <- readRDS("/mnt/raid61/Personal_data/yangqingxin/data/PlasDB_Orthology/Science_genealign_inDiiffstrain.Rds")%>%as.data.table()
setkey(ID_transform,berghei)

transID <- function(x){
  y <- ID_transform[x]$falciparum
  y <- subset(y,!is.na(y))
  return(y)
}

Schizont_508_hostgene_2 <- transID(Schizont_508_hostgene)

library(org.Pf.plasmo.db)
egoBP1 <- enrichGO(Schizont_508_hostgene_2, 
                   OrgDb = org.Pf.plasmo.db,
                   keyType = "SYMBOL",
                   # ont = "BP",
                   pAdjustMethod = "BH")



egoBP2 <- enrichGO(transID(Trophozoite_72_hostgene), 
                   OrgDb = org.Pf.plasmo.db,
                   keyType = "SYMBOL",
                   # ont = "BP",
                   pAdjustMethod = "BH")

JW <- openxlsx::read.xlsx("/mnt/raid61/Personal_data/tangchao/Temp/Supplementary Data 6.xlsx", sheet = 3, colNames = T)
JW <- as.data.table(JW)
colnames(JW) <- as.character(JW[1, ])
JW <- JW[-1, ]

transID(JW[`differential object` == "gene expression", `Gene Name`])
transID(JW[`differential object` != "gene expression", `Gene Name`])

egoBP <- rbind(data.table(Class = "Schizont specific isoform", egoBP1@result),
               data.table(Class = "Trophozoite specific isoform", egoBP2@result))
openxlsx::write.xlsx(egoBP, "/mnt/raid61/Personal_data/tangchao/Temp/Stage_specific_isoform_GO.xlsx")


egoBP_Sig <- egoBP[pvalue < 0.05]
egoBP_Sig <- egoBP_Sig[, c(1, 3:6)]
egoBP_Sig$GeneRatio <- mapply(function(x) eval(parse(text = x)), egoBP_Sig$GeneRatio)
egoBP_Sig[, Stage := gsub(" specific isoform", "", Class)]
setkey(egoBP_Sig, Stage, GeneRatio)
egoBP_Sig$Description <- factor(egoBP_Sig$Description, levels = egoBP_Sig$Description)

ggplot(egoBP_Sig, aes(x = Stage, y = Description, size = GeneRatio, colour = -log10(pvalue))) + 
  geom_point() + 
  scale_colour_gradient2(high = "red", low = "grey") + 
  theme_classic(base_size = 15)
ggsave("/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/analysis/09.AlternativeSplicing/01.Plasmodium/20211025/Stage_specific_novel_isoform_GO.pdf", width = 9.5, height = 4)

