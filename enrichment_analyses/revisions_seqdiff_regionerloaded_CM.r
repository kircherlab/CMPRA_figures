library(dplyr)
library(ggplot2)
library(tidyr)
library(regioneReloaded)
library(BSgenome.Hsapiens.UCSC.hg38)

# since the test need randomizations we have to set a seed
set.seed(42)

args <- commandArgs(TRUE)
type <- args[1]
if (type == "promoter") {
  upbed <- "results/manuscript_figures/revisions/enrichment_differences_sequences/promoter_max_logFC.bed"
  up_name <- "Highest activity"
  down_name <- "Lowest activity"
  downbed <- "results/manuscript_figures/revisions/enrichment_differences_sequences/promoter_min_logFC.bed"
  outfile <- "results/manuscript_figures/revisions/enrichment_differences_sequences/CM_promoter_enrichment_minmax_activity.pdf"
  rds <- "results/manuscript_figures/revisions/enrichment_differences_sequences/CM_promoter_enrichment_minmax_activity.rds"
} else if (type == "CRE") {
  upbed <- "results/manuscript_figures/revisions/enrichment_differences_sequences/CRE_max_logFC.bed"
  up_name <- "Highest activity"
  down_name <- "Lowest activity"
  downbed <- "results/manuscript_figures/revisions/enrichment_differences_sequences/CRE_min_logFC.bed"
  outfile <- "results/manuscript_figures/revisions/enrichment_differences_sequences/CM_CRE_enrichment_minmax_activity.pdf"
  rds <- "results/manuscript_figures/revisions/enrichment_differences_sequences/CM_CRE_enrichment_minmax_activity.rds"
} else {
  stop("Invalid type argument. Use 'promoter' or 'CRE'.")
}

data_up <- read.table(upbed, header = F, sep = "\t", col.names = c("chr", "start", "end"))
data_down <- read.table(downbed, header = F, sep = "\t", col.names = c("chr", "start", "end"))

OE_GR_1 <- GenomicRanges::makeGRangesFromDataFrame(data_up, seqnames.field = "chr",start.field = "start",end.field = "end")
OE_GR_2 <- GenomicRanges::makeGRangesFromDataFrame(data_down, seqnames.field = "chr",start.field = "start",end.field = "end")

GR_OEs <- list(Highest_activity=OE_GR_1,Lowest_activity=OE_GR_2)
# Now we prepare the data of ChIP-seq
# we do the same, read the bed files and turn them into GenomicRanges
features <- c("data/HepG2_H3K4me3_ENCSR575RRX/ENCFF982DUT.bed.gz",
				"data/HepG2_H3K27ac_ENCSR000AMO/ENCFF392KDI.bed.gz",
				"data/HepG2_ATAC_ENCSR042AWH/ENCFF913MQB.bed.gz",
				"data/HepG2_H3K9me3_ENCSR000ATD/ENCFF533JQH.bed.gz",
				"data/HepG2_H3K4me1_ENCSR000APV/ENCFF428FAW.bed.gz",
				"data/HepG2_DNase_ENCSR149XIL/ENCFF209DJG.bed.gz",
				"data/HepG2_H3K27me3_ENCSR000DUE/ENCFF982CSD.bed.gz", 
				"data/HepG2_H3K9ac_ENCSR000AMD/ENCFF358CJW.bed.gz",
				"data/HepG2_H3K79me2_ENCSR000AOM/ENCFF500KVU.bed.gz",
				"data/HepG2_H3K4me2_ENCSR000AMC/ENCFF769MUQ.bed.gz",
				"data/HepG2_H4K20me1_ENCSR000AMQ/ENCFF133PSF.bed.gz",
				"data/HepG2_H2AFZ_ENCSR000AOK/ENCFF474NDN.bed.gz",
				"data/HepG2_H3K36me3_ENCSR000DUD/ENCFF524PMT.bed.gz")
names(features) <- c("H3K4me3","H3K27ac", "ATAC", "H3K9me3", "H3K4me1", "DNase", "H3K27me3", "H3K9ac", "H3K79me2", "H3K4me2", "H4K20me1", "H2AFZ", "H3K36me3")
# then we have to transform them to GenomicRanges
ftlist <- list()
for (i in 1:(length(features)))
{
  df <- read.table(features[i])
  ftlist[[i]] <- GenomicRanges::makeGRangesFromDataFrame(df,seqnames.field = "V1",start.field = "V2",end.field = "V3")
}
names(ftlist) <- names(features)


# now we performe the crosswise permutation test, using 5000 sampling of the data and 100 permutations. 
# The function that we use for randomization if "randomizeRegions" (NOTE: It can take a while, but you 
# # can do it in parallel using the argument mc.cores) and the evalutation function is "numOveralps".

lichic_chip_OE <-  crosswisePermTest(
  Alist = GR_OEs,
  Blist = ftlist,
  sampling = FALSE,
  ranFUN = "randomizeRegions", 
  evFUN = "numOverlaps", 
  ntimes = 1000,
  genome = "hg38",
  mc.cores = 20, 
  mc.set.seed = FALSE
)

# Use clusterize=F when using only two sets of regions
# Finally we generate the matrix that summarise this analysis
lichic_chip_OE <- makeCrosswiseMatrix(lichic_chip_OE, clusterize=F, transform=F) #, transform=T) #, clusterize=T)


saveRDS(lichic_chip_OE, file=rds)
#saveRDS(lichic_chip_OE, file="results/manuscript_figures/figure3/lichic_chip_noeffectBG_OE.rds")

lichic_chip_OE <- readRDS(rds)


pdf(outfile)
# we can use the function to plot implemented in the package or take the matrix and plot it ourselves
plot_enrich_prom <- plotCrosswiseMatrix(lichic_chip_OE)
plot_enrich_prom
dev.off()