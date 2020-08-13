setwd("E:/work/Project/wes")
library(ggplot2)
library(data.table)

#HLA fisher test####
hla <- fread("HLA_result.csv", stringsAsFactors = FALSE)
hla_freq <- fread("hla_freq.txt", stringsAsFactors = FALSE)
hla <- hla[, c("genotype", "name") := tstrsplit(`#Genotype`, "\\*")]#split by *, avoid regular expression
hla <- hla[order(`#Genotype`),]

hla[`#Genotype` %in% hla_freq$hlasim]$Case <- hla_freq$Freq
hla <- hla[!is.na(OR),]
hla <- hla[, c("ca", "co") := list(Case/66, Control/200)]

get_p <- function(y){
  # include here as.numeric to be sure that your values are numeric:
  table <-  matrix(as.numeric(c(y[2], y[3], 66, 200)), ncol = 2, byrow = TRUE)
  if(any(is.na(table))) p <- "error" else p <- fisher.test(table)$p.value
  p
} 
hla$p <- apply(hla, 1, get_p)

get_or <- function(y){
  # include here as.numeric to be sure that your values are numeric:
  table <-  matrix(as.numeric(c(y[2], y[3], 66, 200)), ncol = 2, byrow = TRUE)
  if(any(is.na(table))) p <- "error" else p <- fisher.test(table)$estimate
  p
} 
hla$or <- apply(hla, 1, get_or)

hla <- hla[p > 0.00001]
signi <- hla[p < 0.05, c(1, 2, 3, 11:14)]
signi$Case_non <- 66 - signi$Case
signi$Control_non <- 200 - signi$Control
fwrite(signi, "signi.csv")

#plot
set.seed(9527)
ggsave(filename = "HLA.svg",
       ggplot(hla) + aes(x = genotype, y = -log10(fishers), color = ifelse(fishers < 0.05, "a", "b")) + geom_jitter(size = 2) +
         geom_hline(yintercept = -log10(0.05), size = 0.9, alpha = 0.8, color = "red", linetype = "dashed") +
         geom_text(aes(label = ifelse(fishers < 0.05, `#Genotype`, '')),
                   nudge_x = 1.1,
                   color = "black", size = 3) +
         scale_colour_manual(values=c("#F44336", "#9e9e9e")) +
         labs(x = "HLA genotype", y = expression("-Log"[10]*"P Value")) +
         theme_classic() +
         theme(legend.position="none",
               axis.title.x = element_text(size = 18),
               axis.title.y = element_text(size = 18),
               axis.text.x = element_text(size = 12),
               axis.text.y = element_text(size = 16),
               axis.line.x = element_line(colour = "black", size = 1),
               axis.line.y = element_line(colour = "black", size = 1)
         ),
       width = 30, height = 10, units = "cm"
)

set.seed(9527)
ggsave(filename = "HLA.png",
       ggplot(hla) + aes(x = genotype, y = -log10(fishers), color = ifelse(fishers < 0.05, "a", "b")) + geom_jitter(size = 2) +
         geom_hline(yintercept = -log10(0.05), size = 0.9, alpha = 0.8, color = "red", linetype = "dashed") +
         geom_text(aes(label = ifelse(fishers < 0.05, `#Genotype`, '')),
                   nudge_x = 1.1,
                   color = "black", size = 3) +
         scale_colour_manual(values=c("#F44336", "#9e9e9e")) +
         labs(x = "HLA genotype", y = expression("-Log"[10]*"P Value")) +
         theme_classic() +
         theme(legend.position="none",
               axis.title.x = element_text(size = 18),
               axis.title.y = element_text(size = 18),
               axis.text.x = element_text(size = 12),
               axis.text.y = element_text(size = 16),
               axis.line.x = element_line(colour = "black", size = 1),
               axis.line.y = element_line(colour = "black", size = 1)
         ),
       width = 30, height = 10, dpi = 800, units = "cm"
)

#GWAS####
library(dplyr)
assoc <- read.table("gwas_non005.txt", header = TRUE, sep = "", stringsAsFactors = FALSE)
patient <- read.csv("patient.csv", stringsAsFactors = FALSE)


#Q-Q plot
library(qqman)
library(Cairo)

CairoSVG("qqplot.svg")
par(mar=c(5,6,4,1)+.1)
qq(assoc$P, main = NULL, cex = 0.6, cex.axis = 1, cex.lab = 1)
dev.off()

CairoPNG("qqplot.png", height = 900*4, width = 900*4, dpi = 800)
par(mar=c(5,6,4,1)+.1)
qq(assoc$P, main = NULL, cex = 0.6, cex.axis = 1, cex.lab = 1)
dev.off()


#manhattan
threshold <- 5*10^-8
don <- assoc %>% 
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(assoc, ., by=c("CHR"="CHR")) %>%
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

ggsave(filename = "manhattan.svg",
  ggplot(don) + aes(x = BPcum, y = -log10(P), color = as.factor(CHR)) +
    geom_point(alpha = 0.8, size = 2) +
    geom_hline(yintercept = -log10(threshold), size = 0.6, alpha = 0.8, color="red", linetype="dashed") +
    scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 13), breaks = seq(from = 0, to = 20, by = 2)) +
    theme_classic() +
    labs(x = "Chromosome", y = expression("-Log"[10]*"P Value")) +
    theme(
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.text = element_text(size = 12),
      axis.line.x = element_line(colour = "black", 
                                 size=0.8,
                                 lineend = "butt"),
      axis.line.y = element_line(colour = "black",
                                 size=0.8)
    ),
  width = 30, height = 10, units = "cm"
)

ggsave(filename = "manhattan.png",
       ggplot(don) + aes(x = BPcum, y = -log10(P), color = as.factor(CHR)) +
         geom_point(alpha = 0.8, size = 2) +
         geom_hline(yintercept = -log10(threshold), size = 0.6, alpha = 0.8, color="red", linetype="dashed") +
         scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
         scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
         scale_y_continuous(expand = c(0, 0), limits = c(0, 13), breaks = seq(from = 0, to = 20, by = 2)) +
         theme_classic() +
         labs(x = "Chromosome", y = expression("-Log"[10]*"P Value")) +
         theme(
           legend.position="none",
           panel.border = element_blank(),
           panel.grid.major.x = element_blank(),
           panel.grid.minor.x = element_blank(),
           axis.title.x = element_text(size = 18),
           axis.title.y = element_text(size = 18),
           axis.text = element_text(size = 12),
           axis.line.x = element_line(colour = "black", 
                                      size=0.8,
                                      lineend = "butt"),
           axis.line.y = element_line(colour = "black",
                                      size=0.8)
         ),
       width = 30, height = 10, dpi = 800, units = "cm"
)

#PCA####
mds18 <- fread("trunc18.mds", header = TRUE, sep = " ", stringsAsFactors = FALSE)
mds18$c <- factor(ifelse(nchar(mds18$IID) <= 3, "COVID-19", "Normal"), levels = c("COVID-19", "Normal"))

ggsave(filename = "population_stratification.svg",
       ggplot(mds18) + aes(C1, C2, color = c) + geom_point(size = 2) +
         scale_x_continuous(breaks = c(-0.05, 0, 0.05)) +
         scale_y_continuous(breaks = c(-0.05, 0, 0.05)) +
         labs(x = "Coordinate 1", y = "Coordinate 2") +
         theme_classic() +
         theme(axis.text.x = element_text(size=12),
               axis.text.y = element_text(size=12),
               axis.title.x = element_text(size=14),
               axis.title.y = element_text(size=14),
               legend.title = element_blank(),
               legend.text = element_text(size=14),
               panel.background = element_blank(),
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               axis.line = element_blank(),
               panel.border = element_rect(colour = "black", fill = NA, size = 1),
               plot.margin = unit(c(2.5, 1, 2.5, 1), "cm")
         ),
       width = 20, height = 20, units = "cm"
)

ggsave(filename = "population_stratification.png",
       ggplot(mds18) + aes(C1, C2, color = c) + geom_point(size = 2) +
         scale_x_continuous(breaks = c(-0.05, 0, 0.05)) +
         scale_y_continuous(breaks = c(-0.05, 0, 0.05)) +
         labs(x = "Coordinate 1", y = "Coordinate 2") +
         theme_classic() +
         theme(axis.text.x = element_text(size=12),
               axis.text.y = element_text(size=12),
               axis.title.x = element_text(size=14),
               axis.title.y = element_text(size=14),
               legend.title = element_blank(),
               legend.text = element_text(size=14),
               panel.background = element_blank(),
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               axis.line = element_blank(),
               panel.border = element_rect(colour = "black", fill = NA, size = 1),
               plot.margin = unit(c(2.5, 1, 2.5, 1), "cm")
         ),
  width = 20, height = 20, dpi = 800, units = "cm"
)
