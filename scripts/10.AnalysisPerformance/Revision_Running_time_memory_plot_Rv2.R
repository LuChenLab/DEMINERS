# pidstat -r -u -p 1361 1 > deeplexicon_pidstat2.log
# pidstat -r -u -p 13231 1 > DecodeR_pidstat2.log
library(data.table)
library(ggplot2)

DecodeR_CPU <- lapply(list.files(path = "~", "DecodeR_pidstat", full.names = T), function(x) {
  dat1 <- fread(x, header = F, fill = T)[V4 != "PID" & V10 == "R"][, .SD[1, ], .(V1)]
  colnames(dat1) <- c("Date", "AM", "UID", "PID", "%usr", "%system", "%guest", "%CPU", "CPU", "Command")
  data.table(dat1, Time = seq_len(nrow(dat1)))
})
DecodeR_CPU <- data.table(Method = "DEMIENRS", Rep = rep(seq_along(DecodeR_CPU), mapply(nrow, DecodeR_CPU)), do.call(rbind, DecodeR_CPU))
DecodeR_CPU[, `%CPU` := as.numeric(`%CPU`)]

DecodeR_Mem <- lapply(list.files(path = "~", "DecodeR_pidstat", full.names = T), function(x) {
  dat1 <- fread(x, header = F, fill = T)[V4 != "PID" & V10 == "R"][, .SD[2, ], .(V1)]
  colnames(dat1) <- c("Date", "AM", "UID", "PID", "minflt", "majflt", "VSZ", "RSS", "%MEM", "Command")
  data.table(dat1, Time = seq_len(nrow(dat1)))
})
DecodeR_Mem <- data.table(Method = "DEMIENRS", Rep = rep(seq_along(DecodeR_Mem), mapply(nrow, DecodeR_Mem)), do.call(rbind, DecodeR_Mem))
DecodeR_Mem[, RSS := as.numeric(RSS)]


DeePlexiCon_CPU <- lapply(list.files(path = "~", "deeplexicon_pidstat", full.names = T), function(x) {
  dat1 <- fread(x, header = F, fill = T)[V4 != "PID" & V10 == "python"][, .SD[1, ], .(V1)]
  colnames(dat1) <- c("Date", "AM", "UID", "PID", "%usr", "%system", "%guest", "%CPU", "CPU", "Command")
  data.table(dat1, Time = seq_len(nrow(dat1)))
})
DeePlexiCon_CPU <- data.table(Method = "DeePlexiCon", Rep = rep(seq_along(DeePlexiCon_CPU), mapply(nrow, DeePlexiCon_CPU)), do.call(rbind, DeePlexiCon_CPU))
DeePlexiCon_CPU[, `%CPU` := as.numeric(`%CPU`)]

DeePlexiCon_Mem <- lapply(list.files(path = "~", "deeplexicon_pidstat", full.names = T), function(x) {
  dat1 <- fread(x, header = F, fill = T)[V4 != "PID" & V10 == "python"][, .SD[2, ], .(V1)]
  colnames(dat1) <- c("Date", "AM", "UID", "PID", "minflt", "majflt", "VSZ", "RSS", "%MEM", "Command")
  data.table(dat1, Time = seq_len(nrow(dat1)))
})
DeePlexiCon_Mem <- data.table(Method = "DeePlexiCon", Rep = rep(seq_along(DeePlexiCon_Mem), mapply(nrow, DeePlexiCon_Mem)), do.call(rbind, DeePlexiCon_Mem))
DeePlexiCon_Mem[, RSS := as.numeric(RSS)]


CPU <- rbind(DecodeR_CPU, DeePlexiCon_CPU)
Mem <- rbind(DecodeR_Mem, DeePlexiCon_Mem)

CPU[, Rep := paste0("fast5 file ", Rep)]
Mem[, Rep := paste0("fast5 file ", Rep)]

ggplot(CPU, aes(x = Time, y = `%CPU`, colour = Method, lty = factor(Rep))) + 
  geom_line() + 
  labs(x = "Wall clock time (s)", y = "CPU usage (%)") + 
  scale_color_manual(values = c("#F2913C", "#D5231E")) + 
  theme_classic(base_size = 15) + 
  theme(legend.position = c(0.2, 0.75), 
        legend.title = element_blank()) -> pl1

ggplot(Mem, aes(x = Time, y = (RSS / 1024) / 1024, colour = Method, lty = factor(Rep))) + 
  geom_line() + 
  labs(x = "Wall clock time (s)", y = "Memory consumption (GB)") + 
  theme_classic(base_size = 15) + 
  theme(legend.position = c(0.2, 0.75), 
        legend.background = element_blank(), 
        legend.title = element_blank()) + 
  scale_color_manual(values = c("#F2913C", "#D5231E")) -> pl2

cowplot::plot_grid(pl1, pl2, nrow = 1)



