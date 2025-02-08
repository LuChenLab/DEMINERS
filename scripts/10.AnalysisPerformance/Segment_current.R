---
title: "Using scaled positional signal to train a classifier"
author: "Chao Tang"
date: 'Report created: `r Sys.Date()`'
output: 
  html_document: 
    code_folding: "hide"
    toc: true
    toc_depth: 4
    toc_float: 
      collapsed: false
    number_sections: true
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, 
                      warning = FALSE, 
                      message = FALSE, 
                      fig.align = 'center')
knitr::opts_knit$set(root.dir = "/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex")
```

```{r required packages}
library(PorexploreR)
library(data.table)
library(rhdf5)
library(ggplot2)
library(parallel)
library(changepoint)
library(caret)
library(pbapply)
library(GenomicAlignments)
library(ggbeeswarm)
library(patchwork)
```

```{r}
MyChangePoint2 <- function(sig, MinLength = 10, ChangePoints = 68) {
  cp0 <- suppressWarnings(changepoint::cpt.meanvar(data = sig, 
                                                   Q = ChangePoints, 
                                                   penalty = "Manual", 
                                                   method = "BinSeg", 
                                                   class = FALSE, 
                                                   minseglen = MinLength, 
                                                   param.estimates = FALSE, 
                                                   pen.value = 0.0001)) - 0.5
  bins <- cut(seq_along(sig), c(0, cp0, length(sig)), include.lowest = T, labels = FALSE)
  # bins <- c(0, cp0, length(sig))
  return(as.numeric(bins))
}
```

```{r}
read2gene <- fread("/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/analysis/02.BamReadsSplit/read2gene_20201127.txt")
read2gene[, read := paste0("read_", read)]

dir_bs <- "/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/analysis/03.BarcodeProcessing/01.NormalBarcodeSignal/20201127"
files <- list.files(path = dir_bs, pattern = ".signal", full.names = TRUE)

load(files[1])

read2gene[read %in% names(barcode), .N, Barcode]
read2gene2U <- read2gene[read %in% names(barcode), ]
```

```{r}
sigi_1 <- barcode[[read2gene2U[Barcode == "RTA-03", read][2]]]
sigiM_1 <- data.table(x = seq_along(sigi_1), Current = sigi_1)
sigiM_1$Bin <- MyChangePoint2(sig = sigi_1, ChangePoints = 98, MinLength = 10)
binM_1 <- sigiM_1[, .(BinP = mean(Current)), Bin]
```

```{r}
sigi_2 <- barcode[[read2gene2U[Barcode == "RTA-08", read][1]]]
sigiM_2 <- data.table(x = seq_along(sigi_2), Current = sigi_2)
sigiM_2$Bin <- MyChangePoint2(sig = sigi_2, ChangePoints = 98, MinLength = 10)
binM_2 <- sigiM_2[, .(BinP = mean(Current)), Bin]
```

```{r}
sigi_3 <- barcode[[read2gene2U[Barcode == "RTA-10", read][22]]]
sigiM_3 <- data.table(x = seq_along(sigi_3), Current = sigi_3)
sigiM_3$Bin <- MyChangePoint2(sig = sigi_3, ChangePoints = 98, MinLength = 10)
binM_3 <- sigiM_3[, .(BinP = mean(Current)), Bin]
```

```{r}
sigi_4 <- barcode[[read2gene2U[Barcode == "RTA-21", read][5]]]
sigiM_4 <- data.table(x = seq_along(sigi_4), Current = sigi_4)
sigiM_4$Bin <- MyChangePoint2(sig = sigi_4, ChangePoints = 98, MinLength = 10)
binM_4 <- sigiM_4[, .(BinP = mean(Current)), Bin]
```

```{r}
sigi_5 <- barcode[[read2gene2U[Barcode == "RTA-27", read][3]]]
sigiM_5 <- data.table(x = seq_along(sigi_5), Current = sigi_5)
sigiM_5$Bin <- MyChangePoint2(sig = sigi_5, ChangePoints = 98, MinLength = 10)
binM_5 <- sigiM_5[, .(BinP = mean(Current)), Bin]
```

```{r}
sigi_6 <- barcode[[read2gene2U[Barcode == "RTA-33", read][18]]]
sigiM_6 <- data.table(x = seq_along(sigi_6), Current = sigi_6)
sigiM_6$Bin <- MyChangePoint2(sig = sigi_6, ChangePoints = 98, MinLength = 10)
binM_6 <- sigiM_6[, .(BinP = mean(Current)), Bin]
```

```{r}
sigi_7 <- barcode[[read2gene2U[Barcode == "RTA-37", read][6]]]
sigiM_7 <- data.table(x = seq_along(sigi_7), Current = sigi_7)
sigiM_7$Bin <- MyChangePoint2(sig = sigi_7, ChangePoints = 98, MinLength = 10)
binM_7 <- sigiM_7[, .(BinP = mean(Current)), Bin]
```












```{r fig.width=10, fig.height=2}
ggplot(sigiM_1, aes(x, Current)) + 
  geom_line(color = "#E41A1C") + 
  theme_bw(base_size = 15) + 
  geom_vline(xintercept = cumsum(runLength(Rle(sigiM_1[, Bin])))) + 
  theme(axis.text = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid = element_blank()) + 
  labs(y = "Scaled current")
```

```{r fig.width=10, fig.height=2}
ggplot(binM_1, aes(x = Bin, y = BinP)) + 
  geom_step(color = "#E41A1C") + 
  theme_bw(base_size = 15) + 
  theme(axis.text = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid = element_blank()) + 
  labs(y = "Segment current") -> p1
p1
```


```{r fig.width=10, fig.height=2}
ggplot(binM_2, aes(x = Bin, y = BinP)) + 
  geom_step(color = "#377EB8") + 
  theme_bw(base_size = 15) + 
  theme(axis.text = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid = element_blank()) + 
  labs(y = "Segment current") -> p2
p2
```


```{r fig.width=10, fig.height=2}
ggplot(binM_3, aes(x = Bin, y = BinP)) + 
  geom_step(color = "#4DAF4A") + 
  theme_bw(base_size = 15) + 
  theme(axis.text = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid = element_blank()) + 
  labs(y = "Segment current") -> p3
p3
```


```{r fig.width=10, fig.height=2}
ggplot(binM_4, aes(x = Bin, y = BinP)) + 
  geom_step(color = "#984EA3") + 
  theme_bw(base_size = 15) + 
  theme(axis.text = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid = element_blank()) + 
  labs(y = "Segment current") -> p4
p4
```


```{r fig.width=10, fig.height=2}
ggplot(binM_5, aes(x = Bin, y = BinP)) + 
  geom_step(color = "#984EA3") + 
  theme_bw(base_size = 15) + 
  theme(axis.text = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid = element_blank()) + 
  labs(y = "Segment current") -> p5
p5
```


```{r fig.width=10, fig.height=2}
ggplot(binM_6, aes(x = Bin, y = BinP)) + 
  geom_step(color = "#FF7F00") + 
  theme_bw(base_size = 15) + 
  theme(axis.text = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid = element_blank()) + 
  labs(y = "Segment current") -> p6
p6
```


```{r fig.width=10, fig.height=2}
ggplot(binM_7, aes(x = Bin, y = BinP)) + 
  geom_step(color = "#A65628") + 
  theme_bw(base_size = 15) + 
  theme(axis.text = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid = element_blank()) + 
  labs(y = "Segment current") -> p7
p7
```


```{r fig.width=10, fig.height=5}
p1 / p2 / p3
ggsave("/mnt/raid61/Personal_data/tangchao/Temp/SigCurrent_1.pdf", width = 10, height = 5)
```

```{r fig.width=10, fig.height=5}
p5 / p6 / p7
ggsave("/mnt/raid61/Personal_data/tangchao/Temp/SigCurrent_2.pdf", width = 10, height = 5)
```

```{r}
Mat <- rbind(binM_1[[2]], 
             binM_2[[2]], 
             binM_3[[2]], 
             binM_4[[2]], 
             binM_5[[2]], 
             binM_6[[2]], 
             binM_7[[2]])
row.names(Mat) <- c("Barcode01", "Barcode02", "Barcode03", "Barcode04", "Barcode22", "Barcode23", "Barcode24")
colnames(Mat) <- paste0("BIN", sprintf("%03d", seq_len(ncol(Mat))))
```

