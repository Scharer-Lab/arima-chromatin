# Setup ----

# JM 5/15/2023

library(ggplot2)
library(ggrepel)
library(ggnewscale)
library(tidyverse)
library(dplyr)
library(GenomicRanges)
source("/home/boss_lab/Apps/bitbucket/pca/PCA.theme.module.R") # Magic, plot theme
library(data.table)
library(geomtextpath) # for geom text paths, JM 05/25/2023

## Constants ----
projectName = "Parimal"
projectDir = "/home/boss_lab/Projects/Boss_Projects/HiChIP/Arima_RJ12_Raji1234_09-02-2022/" # Replace, must have '/' at the end
homeDir = paste0(projectDir, "analysis/diff/")
setwd(homeDir)

#source("HiChIPseq.diff.module.R")

manifestDir = paste0(projectDir, "pipeline/")
manifestName = "HiChIPseq.sample.manifest.txt" # Replace
manifest = read.table(paste0(manifestDir, manifestName), sep = "\t", header = T, as.is = T)
manifest = manifest[manifest$include, ]

covName = paste0(covDir, projectName, ".HiChIPseq.common.coverage.txt")
data <- fread(covName)

covDir = paste0(projectDir, "coverage/")
bin = 5000
resolution = as.character(paste0(bin/1000, "kb"))
figSize = PRESENTATION_FIG_SIZE
minCount = 12
bedpe_colnames = c("chr1", "start1", "end1", "chr2", "start2", "end2")

guide_line_alpha = 0.7

# Analysis ----

## MA Plot ----

### Get graphing parameters ----

# Get thresholds, JM 5/25/2023; count and rpm are in the same table
minRajiRpm = data$Raji_2_deep.rpm[data$Raji_2_deep==minCount][1]

limit_lines =
  list(
    geom_vline(xintercept = minRajiRpm, linetype = "dashed", color = "black", alpha = guide_line_alpha),
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = GEOM_LWD, color = "black", alpha = guide_line_alpha),
    geom_hline(yintercept = 1, linetype = "dashed", color = "blue", alpha = guide_line_alpha),
    geom_hline(yintercept =-1, linetype = "dashed", color = "blue", alpha = guide_line_alpha),
    geom_hline(yintercept = max_logFC, color = "red", alpha = guide_line_alpha),
    geom_hline(yintercept = min_logFC, color = "red", alpha = guide_line_alpha),
    geom_vline(xintercept = xmin, color = "red", alpha = guide_line_alpha))

limit_lines_labelled =
  list(
    geom_labelvline(label = "Raji = 12", xintercept = minRajiRpm, size = SMALL_TEXT / .pt, textalpha = 1, hjust = 0.85, linetype = "dashed", color = "black"),
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = GEOM_LWD, color = "black", alpha = guide_line_alpha),
    geom_hline(yintercept = 1, linetype = "dashed", color = "blue", alpha = guide_line_alpha),
    geom_hline(yintercept =-1, linetype = "dashed", color = "blue", alpha = guide_line_alpha),
    geom_labelhline(label = "RJ2.2.5 = 0", yintercept = max_logFC, size = SMALL_TEXT / .pt, textalpha = 1, hjust = 0.92, textcolor = "black", linecolor = "red", linealpha = guide_line_alpha), # line alpha doesn't work?
    geom_labelhline(label = "Raji = 0", yintercept = min_logFC, size = SMALL_TEXT / .pt, textalpha = 1, hjust = 0.92, textcolor = "black", linecolor = "red", linealpha = guide_line_alpha),
    geom_labelvline(label = "Raji = 0", xintercept = xmin, size = SMALL_TEXT / .pt, textalpha = 1, hjust = 0.8, textcolor = "black", linecolor = "red", linealpha = guide_line_alpha))

# Calc diff
data$logFC = log2(data$Raji_2_deep.rpm/data$RJ225_2_deep.rpm)
max_logFC = data[, (max(logFC[is.finite(logFC)]))]
min_logFC = data[, (min(logFC[is.finite(logFC)]))]
data[logFC > max_logFC, logFC := ceil_dec(max_logFC)]
data[logFC < min_logFC, logFC := floor_dec(min_logFC)]
max_logFC_value = max(abs(ceil_dec(max_logFC)), abs(floor_dec(min_logFC)))
ylims = c(-max_logFC_value, max_logFC_value)
xmin = min(data$Raji_2_deep.rpm[data$Raji_2_deep.rpm!=0])

### Whole Genome ----

# Save "blank" without any explicit labels
p <- make_Raji_RJ_MA(data)
if (interactive()) { print(p) }
ggsave("Raji_RJ_rep2_deep.MA.png", p, width = figSize, height = figSize)

# Save "blank" with line labels
pp <- make_Raji_RJ_MA(data, limlines = limit_lines_labelled)
ggsave("Raji_RJ_rep2_deep.MA.labelled_lines.png", pp, width = figSize, height = figSize)

# Add SE points
SE_labels =
  list(geom_point(data = data[data$SE, ], color = "black"),
       geom_label_repel(data = data[data$SE, ],
                       aes(label = paste0(hg38.tssSymbol1, "\n", hg38.tssSymbol2)),
                       max.overlaps = Inf, size = LABEL_TEXT / .pt,
                       bg.r = 0.05, bg.color = "gray", fill = rgb(1,1,1,0.5)))
ggsave("Raji_RJ_rep2_deep.MA.SE.png", p + SE_labels, width = figSize, height = figSize)

# Add CIITA points
data[, CIITA := factor((CIITA1 != "") + (CIITA2 != ""), labels = c("NOT", "XOR", "AND"))]
CIITA_labels =
  list(new_scale_color(),
       geom_point(data = data[CIITA %in% c("XOR", "AND"), ] %>% arrange(CIITA), aes(color = CIITA)),
       scale_color_manual(values = c("gray", "black")),
       theme(legend.margin = margin(c(1,1,3,3))),
       geom_label_repel(data = data[CIITA == "AND", ],
                       aes(label = paste0(hg38.tssSymbol1, "\n", hg38.tssSymbol2)),
                       max.overlaps = Inf, force = 20, min.segment.length = 0,
                       size = LABEL_TEXT / .pt, fill = rgb(1,1,1,0.5)))
ggsave("Raji_RJ_rep2_deep.MA.CIITA.png", p + CIITA_labels, width = figSize, height = figSize)

# Add MHCII loci
data[, MHCII := factor((C1 != "") + (C2 != ""), labels = c("NOT", "XOR", "AND"))]
MHCII_labels =
  list(new_scale_color(),
       geom_point(data = data[MHCII %in% c("XOR", "AND")] %>% arrange(MHCII), aes(color = MHCII)),
       scale_color_manual(values = c("gray", "black")),
       theme(legend.margin = margin(c(1,1,3,3))),
       geom_label_repel(data = data[MHCII == "AND", ],
                       aes(label = paste0(C1, "\n", C2)),
                       max.overlaps = Inf, force = 1, min.segment.length = 0,
                       size = LABEL_TEXT / .pt, fill = rgb(1,1,1,0.5)))
ggsave("Raji_RJ_rep2_deep.MA.MHCII.png", p + MHCII_labels, width = figSize, height = figSize)

# Add MHCII loci, but just the focused ones
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4381081/
focused_sites = c("C1", "XL-4", "XL-7", "XL-9", "C2")
data[, MHCII_focused := factor((C1 %in% focused_sites) + (C2 %in% focused_sites), labels = c("NOT", "XOR", "AND"))]
MHCII_focused_labels =
  list(new_scale_color(),
       geom_point(data = data[MHCII_focused %in% c("XOR", "AND")] %>% arrange(MHCII_focused), aes(color = MHCII_focused)),
       scale_color_manual(values = c("gray", "black"), name = "C1, XL-4, XL-7, \nXL-9, C2"),
       theme(legend.margin = margin(c(1,1,3,3))),
       geom_label_repel(data = data[C1 %in% focused_sites, ], aes(label = C1),
                        max.overlaps = 20, force = 1, min.segment.length = 0,
                        size = LABEL_TEXT / .pt, fill = rgb(1,1,1,0.5)),
       geom_label_repel(data = data[C2 %in% focused_sites, ], aes(label = C2),
                        max.overlaps = 20, force = 1, min.segment.length = 0,
                        size = LABEL_TEXT / .pt, fill = rgb(1,1,1,0.5)),
       geom_label_repel(data = data[MHCII_focused == "AND", ],
                        aes(label = paste0(C1, "\n", C2)),
                        max.overlaps = Inf, force = 5, min.segment.length = 0,
                        size = SMALL_TEXT / .pt))
ggsave("Raji_RJ_rep2_deep.MA.MHCII_focused.png", p + MHCII_focused_labels, width = figSize, height = figSize)

### chr6 ----

# Remove all non-chr6 data from the object named "data" and rerun the previous sections
odata = data
data = data[chr1 == "chr6", ]

p = p %+% data

# Save "blank" without any explicit labels
ggsave("Raji_RJ_rep2_deep.MA.chr6.png", p, width = figSize, height = figSize)

# Add SE points
ggsave("Raji_RJ_rep2_deep.MA.SE.chr6.png", p + SE_labels, width = figSize, height = figSize)

# Add CIITA points
CIITA_labels_chr6 =
  list(new_scale_color(),
       geom_point(data = data[CIITA %in% c("XOR", "AND"), ] %>% arrange(CIITA), aes(color = CIITA)),
       scale_color_manual(values = c("gray", "black")),
       theme(legend.margin = margin(c(1,1,3,3))),
       geom_label_repel(data = data[CIITA == "AND", ],
                        aes(label = paste0(hg38.tssSymbol1, "\n", hg38.tssSymbol2)),
                        max.overlaps = Inf, force = 20, min.segment.length = 0,
                        size = LABEL_TEXT / .pt, fill = rgb(1,1,1,0.5)))
ggsave("Raji_RJ_rep2_deep.MA.CIITA.chr6.png", p + CIITA_labels_chr6, width = figSize, height = figSize)

# Add MHCII loci
MHCII_labels_chr6 =
  list(new_scale_color(),
       geom_point(data = data[MHCII %in% c("XOR", "AND")] %>% arrange(MHCII), aes(color = MHCII)),
       scale_color_manual(values = c("gray", "black")),
       theme(legend.margin = margin(c(1,1,3,3))),
       geom_label_repel(data = data[MHCII == "AND", ],
                        aes(label = paste0(C1, "\n", C2)),
                        max.overlaps = Inf, force = 1, min.segment.length = 0,
                        size = LABEL_TEXT / .pt, fill = rgb(1,1,1,0.5)))
ggsave("Raji_RJ_rep2_deep.MA.MHCII.chr6.png", p + MHCII_labels_chr6, width = figSize, height = figSize)

# Add MHCII loci, but just the focused ones
MHCII_focused_labels_chr6 =
  list(new_scale_color(),
       geom_point(data = data[MHCII_focused %in% c("XOR", "AND")] %>% arrange(MHCII_focused), aes(color = MHCII_focused)),
       scale_color_manual(values = c("gray", "black"), name = "C1, XL-4, XL-7, \nXL-9, C2"),
       theme(legend.margin = margin(c(1,1,3,3))),
       geom_label_repel(data = data[C1 %in% focused_sites, ], aes(label = C1),
                        max.overlaps = 20, force = 1, min.segment.length = 0,
                        size = LABEL_TEXT / .pt, fill = rgb(1,1,1,0.5)),
       geom_label_repel(data = data[C2 %in% focused_sites, ], aes(label = C2),
                        max.overlaps = 20, force = 1, min.segment.length = 0,
                        size = LABEL_TEXT / .pt, fill = rgb(1,1,1,0.5)),
       geom_label_repel(data = data[MHCII_focused == "AND", ],
                        aes(label = paste0(C1, "\n", C2)),
                        max.overlaps = Inf, force = 5, min.segment.length = 0,
                        size = SMALL_TEXT / .pt))
ggsave("Raji_RJ_rep2_deep.MA.MHCII_focused.png", p + MHCII_focused_labels_chr6, width = figSize, height = figSize)

### Basemean ----
# Here we will use an average of the two samples to get a new x-axis.

data = fread(covName) # Undo the logic mess from before

data[, M := log2(data$Raji_2_deep/data$RJ225_2_deep)]
max_M = data[, (max(M[is.finite(M)]))]
min_M = data[, (min(M[is.finite(M)]))]
data[M > max_M, M := ceil_dec(max_M)]
data[M < min_M, M := floor_dec(min_M)]
max_M_value = max(abs(ceil_dec(max_M)), abs(floor_dec(min_M)))
ylims = c(-max_M_value, max_M_value)

data[, A := Raji_2_deep * RJ225_2_deep] # This is like an average
data$A[data$A == 0] = epsilon
data[, A := 0.5 * log2(A)] # As an average, this 0.5 should be 1/n

data[, CIITA := factor((CIITA1 != "") + (CIITA2 != ""), labels = c("NOT", "XOR", "AND"))]
data[, MHCII := factor((C1 != "") + (C2 != ""), labels = c("NOT", "XOR", "AND"))]
data[, MHCII_focused := factor((C1 %in% focused_sites) + (C2 %in% focused_sites), labels = c("NOT", "XOR", "AND"))]
data[, meets_coverage_threshold := factor((Raji_2_deep > 12) + (RJ225_2_deep > 12), labels = c("Neither", "One Sample", "Both"))]

limit_lines_A =
  list(
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = GEOM_LWD, color = "black", alpha = guide_line_alpha),
    geom_hline(yintercept = 1, linetype = "dashed", color = "blue", alpha = guide_line_alpha),
    geom_hline(yintercept =-1, linetype = "dashed", color = "blue", alpha = guide_line_alpha),
    geom_hline(yintercept = max_M, color = "red", alpha = guide_line_alpha),
    geom_hline(yintercept = min_M, color = "red", alpha = guide_line_alpha),
    # x = log2(Raji*RJ)/2; y = log2(Raji/RJ)
    # Solve for Raji & RJ when the other is 12
    # Then set sample_yexpression = sample_xexpression to get the slope and intercept
    geom_abline(slope = 2, intercept = -log2(144), linetype = "dashed", linewidth = GEOM_LWD, color = "black", alpha = guide_line_alpha),
    geom_abline(slope = -2, intercept = log2(144), linetype = "dashed", linewidth = GEOM_LWD, color = "black", alpha = guide_line_alpha))

limit_lines_MA =
  list(
    geom_labelabline(label = "RJ2.2.5 >= 12", hjust = 0.1, slope = 2, intercept = -log2(144), linetype = "dashed", linewidth = GEOM_LWD, color = "black", size = SMALL_TEXT / .pt),
    geom_labelabline(label = "Raji >= 12", hjust = 0.1, slope = -2, intercept = log2(144), linetype = "dashed", linewidth = GEOM_LWD, color = "black", size = SMALL_TEXT / .pt))

p <- make_Raji_RJ_MA_MA(data, limlines = limit_lines_A, colorCol = "meets_coverage_threshold") +
  labs(color = "Coverage >= 12")
ggsave("Raji_RJ_rep2_deep.MA.MA.png", p, width = figSize, height = figSize)

# Save "blank" with line labels
pp <- make_Raji_RJ_MA_MA(data, limlines = limit_lines_A, colorCol = "meets_coverage_threshold") +
  labs(color = "Coverage >= 12") + limit_lines_MA
ggsave("Raji_RJ_rep2_deep.MA.MA.labelled_lines.png", pp, width = figSize, height = figSize)

loopData = data
for (chr in c("All", "chr6")) {
  if (chr == "chr6") {
    loopData = data[chr1 == "chr6", ]
  }
  # Add SE points
  SE_labels =
    list(geom_point(data = loopData[SE == TRUE, ], color = "black"),
         geom_label_repel(data = loopData[SE == TRUE, ],
                          aes(label = paste0(hg38.tssSymbol1, "\n", hg38.tssSymbol2)),
                          max.overlaps = Inf, size = LABEL_TEXT / .pt,
                          bg.r = 0.05, bg.color = "gray", fill = rgb(1,1,1,0.5)))
  ggsave(paste0("Raji_RJ_rep2_deep.MA.MA.", chr, ".SE.png"), p + SE_labels, width = figSize, height = figSize)

  # Add CIITA points
  CIITA_labels =
    list(new_scale_color(),
         geom_point(data = loopData[CIITA %in% c("XOR", "AND"), ] %>% arrange(CIITA), aes(color = CIITA)),
         scale_color_manual(values = c("gray", "black")),
         theme(legend.margin = margin(c(1,1,3,3))),
         geom_label_repel(data = loopData[CIITA == "AND", ],
                          aes(label = paste0(hg38.tssSymbol1, "\n", hg38.tssSymbol2)),
                          max.overlaps = Inf, force = 20, min.segment.length = 0,
                          size = LABEL_TEXT / .pt, fill = rgb(1,1,1,0.5)))
  ggsave(paste0("Raji_RJ_rep2_deep.MA.MA.", chr, ".CIITA.png"), p + CIITA_labels, width = figSize, height = figSize)

  # Add MHCII loci
  MHCII_labels =
    list(new_scale_color(),
         geom_point(data = loopData[MHCII %in% c("XOR", "AND")] %>% arrange(MHCII), aes(color = MHCII)),
         scale_color_manual(values = c("gray", "black")),
         theme(legend.margin = margin(c(1,1,3,3))),
         geom_label_repel(data = loopData[MHCII == "AND", ],
                          aes(label = paste0(C1, "\n", C2)),
                          max.overlaps = Inf, force = 1, min.segment.length = 0,
                          size = LABEL_TEXT / .pt, fill = rgb(1,1,1,0.5)))
  ggsave(paste0("Raji_RJ_rep2_deep.MA.MA.", chr, ".MHCII.png"), p + MHCII_labels, width = figSize, height = figSize)

  # Add MHCII loci, but just the focused ones
  MHCII_focused_labels =
    list(new_scale_color(),
         geom_point(data = loopData[MHCII_focused %in% c("XOR", "AND")] %>% arrange(MHCII_focused), aes(color = MHCII_focused)),
         scale_color_manual(values = c("gray", "black"), name = "C1, XL-4, XL-7, \nXL-9, C2"),
         theme(legend.margin = margin(c(1,1,3,3))),
         geom_label_repel(data = loopData[C1 %in% focused_sites, ], aes(label = C1),
                          max.overlaps = 20, force = 1, min.segment.length = 0,
                          size = LABEL_TEXT / .pt, fill = rgb(1,1,1,0.5)),
         geom_label_repel(data = loopData[C2 %in% focused_sites, ], aes(label = C2),
                          max.overlaps = 20, force = 1, min.segment.length = 0,
                          size = LABEL_TEXT / .pt, fill = rgb(1,1,1,0.5)),
         geom_label_repel(data = loopData[MHCII_focused == "AND", ],
                          aes(label = paste0(C1, "\n", C2)),
                          max.overlaps = Inf, force = 5, min.segment.length = 0,
                          size = SMALL_TEXT / .pt))
  ggsave(paste0("Raji_RJ_rep2_deep.MA.MA.", chr, ".MHCII_focused.png"), p + MHCII_focused_labels, width = figSize, height = figSize)
}

# Normalize quantiles? Maybe another time

## Diff Stats ----
'
## EdgeR ----

library(edgeR)

d <- DGEList(counts=data[, .SD, .SDcols = manifest$name], group=factor(manifest$group))

# Sanity Check
dim(d)
apply(d$counts, 2, sum)

# Filter
keep <- rowSums(cpm(d)>100) >= 2
d <- d[keep,]
d$samples$lib.size <- colSums(d$counts)

# Normalize
d <- calcNormFactors(d)

png(paste0(projectName, ".MDS.BCV.png"))
plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=20)
dev.off()

d1 <- estimateCommonDisp(d, verbose=T)
d1 <- estimateTagwiseDisp(d1)

png(paste0(projectName, ".CommonDispersion.BCV.png"))
plotBCV(d1)
dev.off()

# GLM
design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d <- estimateGLMCommonDisp(d, design.mat)
d <- estimateGLMTrendedDisp(d,design.mat, method="power")
d <- estimateGLMTagwiseDisp(d, design.mat)

png(paste0(projectName, ".GlmDispersion.BCV.png"))
plotBCV(d2)
dev.off()

# Testing!
diff = exactTest(d1)
topTags(diff, n = 10)

de1 <- decideTestsDGE(diff, adjust.method="BH", p.value=0.05)
summary(de1)

de1tags12 <- rownames(d1)[as.logical(de1)]
png(paste0(projectName, ".smear.png"))
plotSmear(et12, de.tags=de1tags12)
abline(h = c(-2, 2), col = "blue")
dev.off()

# LRT???
fit <- glmFit(d2, design.mat)
lrt <- glmLRT(fit)
topTags(lrt, n = 10)

de2 <- decideTestsDGE(lrt, adjust.method="BH", p.value = 0.05)
de2tags12 <- rownames(d2)[as.logical(de2)]
png(paste0(projectName, ".smear.LRT.png"))
  plotSmear(lrt, de.tags=de2tags12)
  abline(h = c(-2, 2), col = "blue")
dev.off()


## diffloop ----

## diffHiC ----

## HICCUPS Diff ----

## multiHiCcompare ----

## selfish ----

# I need to make a merged_loops file by running hiccups on the hic file.
command = paste0("hiccupsdiff Raji.hic RJ225.hic $(hiccups Raji.hic) $(hiccups RJ225.hic) hiccups/")
'

# Scratch ----

