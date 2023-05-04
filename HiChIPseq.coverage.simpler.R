# Description ----

# Just like with all other coverage scripts in the Boss-Scharer lab bitbucket,
# this script takes in a manifest file describing a project's samples and outputs
# several figures describing the data and outputs a tsv file of the coverage values
# for HiChIPseq loops from the ARIMA HiChIPseq MAPS pipeline.

# Setup ----

## Project Information ----

projectDir = "/home/boss_lab/Projects/Boss_Projects/HiChIP/Arima_RJ12_Raji1234_09-02-2022/" # Replace, must have '/' at the end
setwd(paste0(projectDir, "analysis/diff/"))
covDir = paste0(projectDir, "coverage/")

## Dependencies ----

source("HiChIPseq.coverage.simpler.module.R")
source("bedpe.annotate.module.R")

## Manifest ----

manifestDir = paste0(projectDir, "pipeline/")
manifestName = "HiChIPseq.sample.manifest.txt" # Replace
manifest = read.table(paste0(manifestDir, manifestName), sep = "\t", header = T, as.is = T)
manifest = manifest[manifest$include, ]

## Constants ----

# Replace, genomic binning width
bin = 5000

resolution = as.character(paste0(bin/1000, "kb"))

# Replace, genome used to map
genome = "hg38"
num_autosomes = 22 # 22 for human, 19 for mouse

# Replace, output figure size
figSize = PRESENTATION_FIG_SIZE

## Relevant Annotation Files ----

# Replace, if necessary
SEbed = "/home/boss_lab/Projects/Boss_Projects/HiChIP/Arima_RJ_Raji_rep2_Resequenced_10-26-22/coverage_faster/peaks_loops/SE.bed"
CIITAbed = "/home/boss_lab/Projects/Boss_Projects/HiChIP/Arima_RJ_Raji_rep2_Resequenced_10-26-22/coverage_faster/CIITA.bed"
Cbed =  "/home/boss_lab/Projects/Boss_Projects/HiChIP/Arima_Raji_01-22/coverage/C1C16.bed"
# H3K27AcBW = "https://webapp2.pharm.emory.edu/ucsc-downloads/ucsc-data/912.ucsc"
H3K27AcBed = "/home/boss_lab/Apps/genomes/encode_peaks/hg38.H3K27ac.encodePeaks.sorted.bed"

# Coverage ----

## Union set of bin pairs ----
unique_loops = get_unique_loops(manifest)

## Raw Coverage ----
cov = get_raw_coverage(unique_loops, manifest)

## Normalize ----
cov[, (paste0(manifest$name, ".rpm")) := lapply(.SD, \(x) {x * 1e6 / sum(x)}), .SDcols = manifest$name]

## Annotate ----

### Standard  ----

txAnnot = annotBedPEtoTxDbGene(cov, genome = genome,
                               tx = select_ref_txdb(genome), org = select_ref_org(genome),
                               get_cre = TRUE, DT = TRUE)
cov = cbind(cov, txAnnot)

### Specific ----

C = fread(Cbed, col.names = bed_names)
cov[, c("C1", "C2") := bedpe_intersect_bed(cov, C)]

cov[, SE := grepl("DR/DQ-SE", C1) | grepl("DR/DQ-SE", C2)]

CIITA = fread(CIITAbed)
CIITA[, name := paste(chr, start, end, sep = "_")]
cov[, c("CIITA1", "CIITA2") := bedpe_intersect_bed(cov, CIITA)]

H3K27Ac = fread(H3K27AcBed)
setnames(H3K27Ac, paste0("V",1:4), bed_names)
cov[, c("H3K27AcName1", "H3K27AcName2") := bedpe_intersect_bed(cov, H3K27Ac, idcolpe = "ID", idcol = "name")]
H3K27Ac[, ID := paste0(chr, '_', start, '_', end)]
cov[, c("H3K27AcID1", "H3K27AcID2") := bedpe_intersect_bed(cov, H3K27Ac, idcolpe = "ID", idcol = "ID")]

## Write Coverage ----

covBedpeName = paste0(covDir, projectName, ".common.HiChIPseq.bedpe")
fwrite_default(x = unique_loops, file = covBedpeName)

covName = paste0(covDir, projectName, ".HiChIPseq.common.coverage.txt")
fwrite(cov, file = covName, sep = '\t', quote = FALSE)

# Barplots ----

p = plot_sample_loop_count_bar(manifest)
ggsave_default("sample_loop_count.png", p, figSize = figSize)

p = plot_coverage_bar(manifest)
ggsave_default("sample_useful_reads.png", p, figSize = figSize)

# Genomeplot Track Files ----

make_genomeplot_loop_files(covDir, covName = paste0(projectName, ".HiChIPseq.common.coverage.txt"), manifest,
                           conditions = list(SE = bquote((SE == TRUE)),
                                             C = bquote((C1 != "" | C2 != "")),
                                             CIITA = bquote((CIITA1 != "" | CIITA2 != "")),
                                             H3K27Ac = bquote((H3K27AcID1 != "" | H3K27AcID2 != ""))))

# Output R session for reproducibility ----
capture.output(sessionInfo(), file = paste("Rsession.", projectName, gsub("\\D","",Sys.time()), "txt", sep="."))

# Scratch ----

quit()
