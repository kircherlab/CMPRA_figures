library(dplyr)
library(ggplot2)
library(regioneReloaded)
library(BSgenome.Hsapiens.UCSC.hg38)


args <- commandArgs(TRUE)
if (length(args) == 0) {
  stop("Provide either 'chromHMM' or 'SCREEN' as an argument.")
}
experiment = args[1]
if (experiment == "chromHMM") {
  fn = "results/CMPRA5_minP/lichic_chromHMM_chip.rds"
  fig = "results/manuscript_figures/minP/minP_heatmap_chromHMM_enrichment.pdf"
} else if (experiment == "SCREEN") {
  fn = "results/CMPRA5_minP/lichic_SCREEN_chip.rds"
  fig = "results/manuscript_figures/minP/minP_heatmap_SCREEN_enrichment.pdf"
} else {
  stop("Invalid argument. Use 'chromHMM' or 'SCREEN'.")
} 

upbed <- "results/CMPRA5_minP/higher_CMPRA_CREs.bed"
downbed <- "results/CMPRA5_minP/lower_CMPRA_CREs.bed"
backgroundbed <- "results/CMPRA5_minP/similar_CREs.bed"

data_up <- read.table(upbed, header = F, sep = "\t", col.names = c("chr", "start", "end", "name", "strand", "score"))
data_down <- read.table(downbed, header = F, sep = "\t", col.names = c("chr", "start", "end", "name", "strand", "score"))
data_bg <- read.table(backgroundbed, header = F, sep="\t", col.names = c("chr", "start", "end", "name", "strand", "score"))

OE_GR_1 <- GenomicRanges::makeGRangesFromDataFrame(data_up, seqnames.field = "chr",start.field = "start",end.field = "end")
OE_GR_2 <- GenomicRanges::makeGRangesFromDataFrame(data_down, seqnames.field = "chr",start.field = "start",end.field = "end")
OE_GR_3 <- GenomicRanges::makeGRangesFromDataFrame(data_bg, seqnames.field = "chr",start.field = "start",end.field = "end")
GR_OEs <- list(high_PSI=OE_GR_1,low_PSI=OE_GR_2, mid_PSI=OE_GR_3)

# Now we prepare the data of ChIP-seq
# we do the same, read the bed files and turn them into GenomicRanges
if (experiment == "chromHMM") {
	df <- read.table("data/HepG2_chromHMM_ENCSR753TNH/E118_15_coreMarks_hg38lift_dense_4cols.bed.gz", header = F, sep = "\t", col.names = c("chr", "start", "end", "state"), skip = 1)
} else {
	df <- read.table("resources/GRCh38-cCREs_4col.bed", header = F, sep = "\t", col.names = c("chr", "start", "end", "state"))
}
states <- unique(df$state)

# then we have to transform them to GenomicRanges
ftlist <- list()
for (i in 1:(length(states)))
{
  curr <- df[df$state == states[i], ]
  ftlist[[i]] <- GenomicRanges::makeGRangesFromDataFrame(curr)
}

if (experiment == "chromHMM") {
	states_full <- recode(states, 
	"1_TssA" = "Active TSS",
	"2_TssAFlnk" = "Flanking Active TSS	Orange",
	"3_TxFlnk" = "Transcr. at gene 5' and 3'",
	"4_Tx" = "Strong transcription",
	"5_TxWk" = "Weak transcription",
	"6_EnhG" = "Genic enhancers",
	"7_Enh" = "Enhancers",
	"8_ZNF" = "Rpts	ZNF genes & repeats	Medium",
	"9_Het" = "Heterochromatin",
	"10_TssBiv" = "Bivalent/Poised TSS",
	"11_BivFlnk" = "Flanking Bivalent TSS/Enh",
	"12_EnhBiv" = "Bivalent Enhancer",
	"13_ReprPC" = "Repressed PolyComb",
	"14_ReprPCWk" = "Weak Repressed PolyComb",
	"15_Quies" = "Quiescent/Low"
	)
	names(ftlist) <- states_full
} else {
	names(ftlist) <- states
}




# since the test need randomizations we have to set a seed
set.seed(42)

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
  mc.cores = 20
)

# Use clusterize=F when using only two sets of regions
# Finally we generate the matrix that summarise this analysis
lichic_chip_OE <- makeCrosswiseMatrix(lichic_chip_OE, clusterize=T, transform=T) #, transform=T) #, clusterize=T)


saveRDS(lichic_chip_OE, file=fn)
#saveRDS(lichic_chip_OE, file="results/manuscript_figures/figure3/lichic_chip_noeffectBG_OE.rds")

lichic_chip_OE <- readRDS(fn)


pdf(fig)
# we can use the function to plot implemented in the package or take the matrix and plot it ourselves
plot_enrich_prom <- plotCrosswiseMatrix(lichic_chip_OE)
plot_enrich_prom
dev.off()
