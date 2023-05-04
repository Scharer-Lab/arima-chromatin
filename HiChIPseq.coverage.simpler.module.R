# JM 4/24/23

library(GenomicRanges) # library(IRanges)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(latticetools) # get_repelled_labels()
library(eulerr)
library(stringr) # str_count()
library(ggrepel)
source("/home/boss_lab/Apps/bitbucket/pca/PCA.theme.module.R") # Magic, plot theme
library(data.table)

# Constants ----

projectName = "Parimal"

bed_names = c("chr", "start", "end", "name")
bedpe_colnames = c("chr1", "start1", "end1", "chr2", "start2", "end2")
loop_colnames = c("chr", "start1", "end1", "start2", "end2", "score", "name", "notes")
bedpexbed_colnames = c(bedpe_colnames, "chr_bed", "start_bed", "end_bed")

## Themes ----

barplot_aesthetics =
  list(theme_classic(),
       scale_y_continuous(expand = expansion(mult = c(0, NA))),
       theme(text = element_text("black"),
             axis.text = element_text(color = "black", size = REGULAR_TEXT),
             axis.text.x = element_text(angle = -45, hjust = 0),
             axis.title.y = element_text(size = LARGE_TEXT),
             axis.ticks.y = element_line(color = "black"),
             axis.ticks.x = element_blank(),
             panel.grid.major = element_line()))

# Functions ----

## Helper Functions ----

catn = \(x) {cat(paste0(x, '\n'))}

fwrite_default = \(...) { fwrite(sep = '\t', quote = FALSE, ...) }

ggsave_default = \(filename, plot, figSize, ...) { ggsave(filename, plot, width = figSize, height = figSize, ...)}

powerset = function (x, max, min = 0, rev = FALSE) { # inspired by rje::powerSet()
  if (missing(max)) max = length(x)
  if (max == 0) return(list(x[c()]))
  if (max > length(x)) {warning("Parameter 'max' should be less than the length of 'x'. Continuing.")}
  out = list(x[c()])
  if (length(x) == 1) return(c(out, list(x)))
  for (i in seq_along(x)) {
    if (rev) out = c(lapply(out[lengths(out) < max], \(y) c(y, x[i])), out)
    else out = c(out, lapply(out[lengths(out) < max], \(y) c(y, x[i])))
  }
  out[lengths(out)>=min]
}

bed_intersect_bed = function(bed_source, bed_query, sep = ';', idcol = "name") {

  # Define the ranges of interest
  bed_sourceGR = GRanges(seqnames = bed_source$chr, ranges = IRanges(start = bed_source$start, end = bed_source$end, names = bed_source$name))
  bed_queryGR = GRanges(seqnames = bed_query$chr, ranges = IRanges(start = bed_query$start, end = bed_query$end, names = bed_query$name))

  # Find the overlaps
  hits = findOverlaps(bed_queryGR, bed_sourceGR)

  # Handle many-to-many
  hitMap = aggregate(bed_query[[idcol]][queryHits(hits)] ~ subjectHits(hits),
                     FUN = paste, collapse = sep, drop = FALSE)

  # Make annotation column
  out = character(nrow(bed_source))
  out[hitMap[,1]] = hitMap[,2]

  return(out)
}

bedpe_intersect_bed = function(bedpe, bed, idcolpe = "ID", ...) {

  # Define the ranges of interest
  .halve_bedpe_names = \(i) { c(paste0(bed_names[1:3], i), idcolpe) }
  loop_start = bedpe[, .halve_bedpe_names(1), with = FALSE]
  loop_end = bedpe[, .halve_bedpe_names(2), with = FALSE]

  setnames(loop_start, bed_names)
  setnames(loop_end, bed_names)

  # Annotate with the overlaps
  bedpe_annot = list(bed_intersect_bed(loop_start, bed, ...),
                     bed_intersect_bed(loop_end, bed, ...))
  return(bedpe_annot)
}

## Plotting Functions----

plot_sample_loop_count_bar = \(manifest) {
  ggplot(manifest, aes(x = name, y = loop_count)) +
    geom_col() +
    labs(y = "Loops", x = "Sample") +
    barplot_aesthetics
}

plot_coverage_bar = \(manifest) {
  ggplot(manifest, aes(x = name, y = mapped_and_deduped/1e6)) +
    geom_col() +
    labs(y = "Reads (Millions)", x = "Sample") +
    barplot_aesthetics
}

draw_eulerr_repel = function(unique_loop_sample, samples, outfile, shortening = 0.03, force_pull = 0.001, point_padding_x = 0.05, ...) {
  all_comb = sapply(powerset(samples, min = 1), paste, collapse = "&") # Get rid of empty subset
  loop_origins = c(table(unique_loop_sample))
  missing_comb = setdiff(all_comb, names(loop_origins))
  loop_origins = c(loop_origins, setNames(integer(length(missing_comb)), missing_comb))
  loop_origins = loop_origins[order(stringr::str_count(names(loop_origins), "&"), names(loop_origins))]

  fit = euler(loop_origins, shape = "ellipse")

  # First plot is to see where I need to label
  p = plot(fit,
           labels = TRUE, #list(fontface = "bold", label = sort(manifest$name)), # Unnecessary, but demonstrative
           quantities = list(labels = fit$original.values),
           main = list(labels = "Sample Loops", fontface = "bold", vjust = 1.5),
           fill = rev(gg_color_hue(nrow(manifest))),
           edges = FALSE,
           legend = TRUE)
  repelled_labs = with(p$data$centers,
                       get_repelled_labels(x, y, labels = quantities,
                                           xlim = p$data$xlim, ylim = p$data$ylim,
                                           force_pull = force_pull, point_padding_x = point_padding_x, ...))
  repelled_labs = repelled_labs[, c("x_orig","y_orig","x","y", "label")]
  xcoords = repelled_labs[, c("x_orig", "x")]
  repelled_labs[, c("x_orig", "x")] = (xcoords - min(p$data$xlim)) / (max(p$data$xlim) - min(p$data$xlim))

  ycoords = repelled_labs[, c("y_orig", "y")]
  repelled_labs[, c("y_orig", "y")] = (ycoords - min(p$data$ylim)) / (max(p$data$ylim) - min(p$data$ylim))

  repelled_labs$d_orig = sqrt((repelled_labs$x-repelled_labs$x_orig)^2 +
                                (repelled_labs$y-repelled_labs$y_orig)^2)
  n = 1 - (shortening / repelled_labs$d_orig) # Ratio to reduce distance by 0.3
  repelled_labs$segment_y = with(repelled_labs, (1-n)*y_orig + n*y)
  repelled_labs$segment_x = with(repelled_labs, (1-n)*x_orig + n*x)

  pp = plot(fit,
            label = TRUE,
            main = list(labels = "Sample Loops", fontface = "bold", vjust = 1.5),
            fill = rev(gg_color_hue(nrow(manifest))),
            edges = FALSE,
            legend = TRUE)
  pdf(outfile, width = 10) # default 7, "venn.HiChIPseq.loops.repel.pdf"
    pp
    seekViewport(pp$children$canvas.grob$vp$name)
    if (any(n>0)) {
      with(repelled_labs,
           grid::grid.segments(x0 = x_orig, x1 = segment_x,
                               y0 = y_orig, y1 = segment_y))
    }
    with(repelled_labs, grid::grid.text(x, y, label = label)) # fontsize is 12
  dev.off()
}

## Workhorse Functions ----

get_unique_loops = function(manifest) {
  loop_files = setNames(as.list(paste0(manifest$dir, manifest$loop_file)), manifest$name)
  sig_bedpe = rbindlist(sapply(loop_files, \(loop_file) { fread(loop_file) },
                               simplify = FALSE, USE.NAMES = TRUE),
                        idcol = "sample")
  sig_bedpe %>%
    reframe(sample = paste(unique(sample), collapse = "&"),
          .by = all_of(bedpe_colnames)) %>%
    mutate(ID = paste0(chr1, ":", start1, "-", end1, "_",
                           chr2, ":", start2, "-", end2)) %>%
    setDT %>% setcolorder(c(bedpe_colnames, "ID", "sample"))
}

get_raw_coverage = function(cov, manifest, verbose = 0) {
  for (sample in 1:nrow(manifest)) {
    if (verbose > 0) { catn(paste("Sample:", manifest$sample[sample], Sys.time())) }
    sampCovSumPath = paste0(manifest$dir[sample], manifest$MAPS_output_dir[sample])
    sampCovSumPre = data.table()
    for (chr in 1:num_autosomes) {
      if (verbose > 1) { catn(paste("\tchr:", chr)) }
      sampCovSumPaths = paste0(sampCovSumPath, list.files(sampCovSumPath, pattern = paste0("reg_raw.chr", chr, "\\..*(xor|and)$")))
      if (verbose > 1) { catn(paste("\tfiles:", paste(sampCovSumPaths, collapse = ""))) }
      sampCovSumChr = rbindlist(lapply(sampCovSumPaths, fread))
      sampCovSumChr$chr = paste0("chr", chr)
      sampCovSumPre = rbind(sampCovSumChr, sampCovSumPre)
    }
    sampCovSum =
      sampCovSumPre %>%
      reframe(chr1 = chr, start1 = bin1_mid, end1 = bin1_mid + bin,
              chr2 = chr, start2 = bin2_mid, end2 = bin2_mid + bin,
              count)
    cov = merge.data.table(cov, sampCovSum, all.x = TRUE)
    cov$count[is.na(cov$count)] = 0
    setnames(cov, old = "count", new = manifest$name[sample])
  }

  return(cov)
}

make_genomeplot_loop_files =
  function(covDir, covName, manifest, conditions = list(all = bquote((TRUE))),
           min_genomeplot_score = 0.05, range_dist = 400*1000, verbose = 0) {

  # Load data
  rpm = fread(paste0(covDir, covName))
  rpm[, `:=`(called_samples = sample,
             score = 1,
             name = paste0("p", 1:nrow(rpm)),
             notes = "Nothing")]
  setnames(rpm, "chr1", "chr")

  # Create folder to save output
  fp = paste0(covDir, "feature_loops/")
  dir.create(fp, showWarnings = FALSE)

  # Make tracks per sample
  for (msample in manifest$name) {
    if (verbose > 0) { catn(paste("Sample:", msample, Sys.time())) }
    sampRPM = rpm[called_samples %flike% msample, ]
    dist = sampRPM$start2 - sampRPM$start1
    sampRPM$score = (dist - min(dist)) / (range_dist - min(dist))
    sampRPM$score = sampRPM$score * (1 - min_genomeplot_score) + min_genomeplot_score
    for (condition in names(conditions)) {
      if (verbose > 1) { catn(paste("\tCondition:", condition)) }
      sampLoops = sampRPM[eval(conditions[[condition]]), ..loop_colnames]
      outfile = paste0(fp, projectName, ".", msample, ".HiChIPseq.", condition, ".loop")
      fwrite_default(x = sampLoops, file = outfile)
    }
    if (verbose > 1) { catn("\tCondition: Everything Else") }
    nc = str2lang(paste0("!(", paste(sapply(conditions, deparse), collapse = " | "), ")"))
    sampLoops = sampRPM[eval(nc), ..loop_colnames]
    outfile = paste0(fp, projectName, ".", msample, ".HiChIPseq.other.loop")
    fwrite_default(x = sampLoops, file = outfile)
  }
}

# Scratch ----

# myC <- function(...) {
#   ll <- list(...)
#   ll <- lapply(ll, function(X) as.list(X))
#   ll <- do.call("c", ll)
#   as.call(c(as.symbol("and"), ll))
# }
# myC(a,b)
# fwrite_default(x = sampRPM[SE == TRUE, ..loop_colnames],
#                file = paste0(fp, projectName, ".", sample, ".HiChIPseq.SE.loop")) # Red
# fwrite_default(x = sampRPM[CIITA1 != "" | CIITA2 != "", ..loop_colnames],
#                file = paste0(fp, projectName, ".", sample, ".HiChIPseq.CIITA.loop")) # Blue
# fwrite_default(x = sampRPM[MHCII1 != "" | MHCII2 != "", ..loop_colnames],
#                file = paste0(fp, projectName, ".", sample, ".HiChIPseq.MHCII.loop")) # Green
# fwrite_default(x = sampRPM[!(SE | CIITA %in% c("XOR1", "XOR2", "AND") | MHCII1 != "" | MHCII2 != ""), ..loop_colnames],
#                file = paste0(fp, projectName, ".", sample, ".HiChIPseq.anyOther.loop")) # Black
# # library(basicPlotteR) # addText...
# library(tikzDevice) # gridToDevice()

# overlap_bedpe_bed = \(bedpe, bed, out, cut = "| cut -f 1-6") {
#   pair2bedCommand = paste("bash -c 'pairToBed -a", bedpe, "-b", bed, cut, ">", out, "'")
#   system(pair2bedCommand)
#   fread(out, header = F)
# }
#
# bedpe_intersect_bed = function(bedpe, bed) {
#
#   # Define the ranges of interest
#   loop_start = GRanges(seqnames = bedpe$chr1, ranges = IRanges(start = bedpe$start1, end = bedpe$end1, names = bedpe$name))
#   loop_end = GRanges(seqnames = bedpe$chr2, ranges = IRanges(start = bedpe$start2, end = bedpe$end2, names = bedpe$name))
#   bedGR = GRanges(seqnames = bed$chr, IRanges(start = bed$start, end = bed$end, names = bed$name))
#
#   # Find the overlaps
#   startHits = findOverlaps(loop_start, bedGR)
#   endHits = findOverlaps(loop_end, bedGR)
#
#   # Map query to subject, many-to-many
#   startMap = aggregate(bed$name[subjectHits(startHits)] ~ queryHits(startHits), FUN = paste, collapse = ";", drop = FALSE) # data.frame(q = queryHits(startHits), s = subjectHits(startHits))
#   endMap = aggregate(bed$name[subjectHits(endHits)] ~ queryHits(endHits), FUN = paste, collapse = ";", drop = FALSE)
#
#   # Annotate with the overlaps
#   bedpe$name1[startMap] = startMap$Val
#   bedpe$name2[endMap] = endMap$Val
#
#   return(bedpe)
# }

#'
#' #' Given a Set of Points and Box sizes, find locations
#' findboxes <- function(df, xcol = "x", ycol = "y", xlim, ylim, point_size = 1,
#'                       pad_point_x = 0, pad_point_y = 0, hjust = 0, vjust = 0, ...) {
#'
#'   #x and y posiitons as a dataframe
#'   posdf <- df[c(xcol,ycol)]
#'
#'   #returnd a df where columns are points
#'   boxdf <- apply(posdf,1,function(row) {
#'     xval <- row[xcol]
#'     yval <- row[ycol]
#'     return(c(xval,
#'              yval,
#'              xval + pad_point_x,
#'              yval + pad_point_y))})
#'   # columns are x1,y1,x2,y2
#'   boxmatrix = as.matrix(t(boxdf))
#'
#'   moved <- ggrepel:::repel_boxes2(data_points = as.matrix(posdf),
#'                                   point_padding_x = pad_point_x,
#'                                   point_padding_y = pad_point_y,
#'                                   point_size = rep_len(point_size, length.out = nrow(boxmatrix)),
#'                                   boxes = boxmatrix,
#'                                   xlim = xlim, ylim = ylim,
#'                                   hjust = rep_len(hjust, length.out = nrow(boxmatrix)),
#'                                   vjust = rep_len(vjust, length.out = nrow(boxmatrix)))
#'
#'   finaldf <- cbind(posdf, moved)
#'   names(finaldf) <- c("x1","y1","x2","y2")
#'   return(finaldf)
#' }
#'
#'
#' df1 <- findboxes(mtcars, xcol = 'wt', ycol='mpg',
#'                  pad_point_x = 0.5,pad_point_y = 0.5,
#'                  xlim = c(0,8), ylim=c(10,35))
#
# ggplot(df1) +
#   geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2)) +
#   geom_point(aes(x1,y1), color="black") +
#   geom_point(aes(x2,y2), color="red")

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



