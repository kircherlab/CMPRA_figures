library(dplyr)
library(ggplot2)
library(regioneReloaded)
library(BSgenome.Hsapiens.UCSC.hg38)


fn = "results/manuscript_figures/figure3/lichic_refseq_OE.rds"
fig = "results/manuscript_figures/figure3/heatmap_refseq_enrichment.pdf"


upbed <- "results/manuscript_figures/CRE_enrichment_plot/upregulating_distal_OEs_absolutely_nonprom.bed"
midbed <- "results/manuscript_figures/CRE_enrichment_plot/no_effect_distal_OEs_absolutely_nonprom.bed"
downbed <- "results/manuscript_figures/CRE_enrichment_plot/downregulating_distal_OEs_absolutely_nonprom.bed"

data_up <- read.table(upbed, header = F, sep = "\t", col.names = c("chr", "start", "end", "name", "strand", "score"))
data_down <- read.table(downbed, header = F, sep = "\t", col.names = c("chr", "start", "end", "name", "strand", "score"))
data_mid <- read.table(midbed, header = F, sep="\t", col.names = c("chr", "start", "end", "name", "strand", "score"))

OE_GR_1 <- GenomicRanges::makeGRangesFromDataFrame(data_up, seqnames.field = "chr",start.field = "start",end.field = "end")
OE_GR_2 <- GenomicRanges::makeGRangesFromDataFrame(data_down, seqnames.field = "chr",start.field = "start",end.field = "end")
OE_GR_3 <- GenomicRanges::makeGRangesFromDataFrame(data_mid, seqnames.field = "chr",start.field = "start",end.field = "end")
GR_OEs <- list(enhancer=OE_GR_1,silencer=OE_GR_2, neutral=OE_GR_3)


df <- read.table("/data/cephfs-2/unmirrored/groups/kircher/MPRA/CaptureCMPRA/resources/gencode.v44.introns_exons_UTRs_intergenic.bed.gz", header = F, sep = "\t", col.names = c("chr", "start", "end", "state"))
states <- unique(df$state)

# then we have to transform them to GenomicRanges
ftlist <- list()
for (i in 1:(length(states)))
{
  curr <- df[df$state == states[i], ]
  ftlist[[i]] <- GenomicRanges::makeGRangesFromDataFrame(curr)
}

names(ftlist) <- states





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
