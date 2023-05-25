# JM 3/20/23

library(data.table)
library(tidyverse)

## Constants ----

epsilon = 0.1

options(warn = -1) # Ignore warnings
da = dup_axis(labels = NULL, name = NULL, breaks = NULL)

## Functions ----

make_Raji_RJ_MA <- function(data, limlines = limit_lines, x = "Raji_2_deep.rpm", xlabel = "RPM Raji", ylabel = "log2FC (Raji/RJ2.2.5)") {
  ggplot(data, aes(x = .data[[x]], y = logFC)) +
    geom_point(aes(color = RJ225_2_deep > minCount)) +
    labs(x = xlabel, y = ylabel,
         color = "RJ Count >= 12",
         title = "MA plot (Raji vs RJ2.2.5)") +
    limlines +
    theme_classic() + PCA_theme +
    theme(panel.grid.major = element_line(),
          plot.margin = margin(0,1,0,6),
          axis.title.y = element_text(vjust = 0)) +
    scale_y_continuous(limits = ylims, sec.axis = da, breaks = seq(from = ylims[1], to = ylims[2])) +
    scale_x_continuous(trans = "log10", sec.axis = da)
}

make_Raji_RJ_MA_MA <- function(data, limlines, colorCol = NULL,
                               xlabel = "A (log2(Raji * RJ2.2.5) / 2)",
                               ylabel = "M (log2(Raji/RJ2.2.5))") {
  p = ggplot(data, aes(x = A, y = M))
  if (is.null(colorCol)) {
    p = p + geom_point(color = "black")
  } else {
    p = p + geom_point(aes(color = .data[[colorCol]]))
  }
  p +
    labs(x = xlabel, y = ylabel, title = "MA plot (Raji vs RJ2.2.5)") +
    limlines +
    theme_classic() + PCA_theme +
    theme(panel.grid.major = element_line(),
          plot.margin = margin(0,1,0,6),
          axis.title.y = element_text(vjust = 0)) +
    scale_y_continuous(limits = ylims, sec.axis = da, breaks = seq(from = ylims[1], to = ylims[2])) +
    scale_x_continuous(sec.axis = da)
}

## Relevant Annotation Files ----

SEbed = "/home/boss_lab/Projects/Boss_Projects/HiChIP/Arima_RJ_Raji_rep2_Resequenced_10-26-22/coverage_faster/peaks_loops/SE.bed"
CIITAbed = "/home/boss_lab/Projects/Boss_Projects/HiChIP/Arima_RJ_Raji_rep2_Resequenced_10-26-22/coverage_faster/CIITA.bed"

# Scratch ----
#
# get_clusters = function(loop) {
#   loop %>%
#     group_by(ClusterLabel) %>%
#     reframe(chr1 = dplyr::first(chr1),
#             start1 = min(start1), end1 = max(end1),
#             chr2 = dplyr::first(chr2),
#             start2 = min(start2), end2 = max(end2),
#             count = sum(count), nBinPairs = n())
# }
#
# loop_list_summary = function(loop, name = "") {
#   print(paste("Number of Loops:", nrow(loop)))
#   png(paste("start", name, "png", sep = "."))
#     hist(loop$end1 - loop$start1) # Thickness of start
#   dev.off()
#   png(paste("end", name, "png", sep = "."))
#     hist(loop$end2 - loop$start2) # Thickness of end
#   dev.off()
#   dists = (loop$end2 + loop$start2) / 2 - (loop$end1 + loop$start1) / 2
#   png(paste("distance", name, "png", sep = "."))
#     hist(dists) # Distance between ends of the loop
#   dev.off()
#   png(paste("logDistance", name, "png", sep = "."))
#     hist(log10(dists)) # Distance between ends of the loop
#   dev.off()
# }
#
# check_cluster2_arrangements = function(loop, name = "name", bin = 5000) {
#   # For my own sanity regarding the "get_clusters() function
#   xadjustFactor = 5
#   yadjustFactor = 3.3
#   xadjust = bin / xadjustFactor
#   yadjust = bin / yadjustFactor
#   displayBinWidthFactor = 1.5
#
#   # Data to plot
#   paired_clusters = loop %>%
#     group_by(ClusterLabel) %>%
#     filter(n() == 2) %>% # Magic
#     reframe(distStart = dplyr::last(start1) - dplyr::first(start1), # in Clusters, replace start1 with the midpoint of the loopStart or loopEnd
#             distEnd = dplyr::last(start2) - dplyr::first(start2),
#             distStartEnd = dplyr::first(((end2 + start2) / 2 - (end1 + start1) / 2)))
#   unique_distance_pairs = paired_clusters %>% group_by(distStart, distEnd) %>% reframe(count = n())
#   minDistStart = min(paired_clusters$distStart) * displayBinWidthFactor
#   maxDistEnd = max(paired_clusters$distEnd) * displayBinWidthFactor
#
#   p = ggplot(paired_clusters, aes(x = distStart, y = distEnd)) +
#     scale_x_continuous(n.breaks = length(unique(paired_clusters$distStart))) +
#     scale_y_continuous(n.breaks = length(unique(paired_clusters$distEnd))) +
#     geom_count(aes(fill = after_stat(n)), color = "black", stroke = 1, size = 20, shape = 22) +
#     scale_fill_gradient(low = "white", high = "red", name = "Frequency") +
#     geom_text(data = unique_distance_pairs,
#               mapping = aes(label = count)) +
#     labs(x = "Distance, Starts of Loops (bp)",
#          y = "Distance, Ends of Loops (bp)",
#          title = "Configuration of the Bin Pairs\nfor Clusters with Two sets of Bin Pairs") +
#     geom_curve(data = unique_distance_pairs,
#                mapping = aes(x = distStart - xadjust, xend = distStart + xadjust,
#                              y = distEnd + yadjust, yend = distEnd + yadjust), curvature = -1) +
#     geom_curve(data = unique_distance_pairs,
#                mapping = aes(x = distStart - xadjust * (1 + (distStart / minDistStart)),
#                              xend = distStart + xadjust * (1 + (distEnd / maxDistEnd)),
#                              y = distEnd + yadjust, yend = distEnd + yadjust), curvature = -1, color = "darkgray") +
#     theme_classic() +
#     theme(plot.title = element_text(hjust = 0.5),
#           axis.text = element_text(color = "black"))
#   ggsave(paste0(name, ".cluster2.configuration.bp.png"), p)
#
#   # Data to plot # TODO: this isn't right
#   paired_clusters = paired_clusters %>%
#     mutate(distStart = distStart / bin,
#           distEnd = distEnd / bin,
#           distStartEnd = distStartEnd / bin)
#   unique_distance_pairs = paired_clusters %>% group_by(distStart, distEnd) %>% reframe(count = n())
#
#   p = p +
#     scale_x_continuous(breaks = length(unique(paired_clusters$distStart)),
#                        labels = unique(paired_clusters$distStart) / bin) +
#     scale_y_continuous(breaks = length(unique(paired_clusters$distEnd)),
#                        labels = unique(paired_clusters$distEnd) / bin) +
#     labs(x = "Distance, Starts of Loops (bin)",
#          y = "Distance, Ends of Loops (bin)",
#          title = "Configuration of the Bin Pairs\nfor Clusters with Two sets of Bin Pairs") +
#   ggsave(paste0(name, ".cluster2.configuration.bin.png"), p)
# }
#
# loop_annotation_summary = function (loop, name = "") {
#   if (!all(c("SYMBOL1", "SYMBOL2") %in% names(loop))) { print("This loop object is not annotated correctly. Quitting."); return()}
#   analyze_list = list(just_start = loop$SYMBOL1,
#                       just_end = loop$SYMBOL2,
#                       start_end_independent = c(loop$SYMBOL1, loop$SYMBOL2),
#                       by_loop = paste0(loop$SYMBOL1, "_", loop$SYMBOL2))
#   for (type in names(analyze_list)) {
#     gene_count_table = sort(table(analyze_list[[type]]), decreasing = TRUE)
#     print(type)
#     print(head(gene_count_table, 10))
#     write.table(gene_count_table, file = paste(type, name, "txt", sep = "."), row.names = FALSE, col.names = FALSE, quote = FALSE)
#   }
# }
