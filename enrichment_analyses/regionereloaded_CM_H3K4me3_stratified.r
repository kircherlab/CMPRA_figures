library(dplyr)
library(ggplot2)
library(tidyr)
library(regioneReloaded)
library(BSgenome.Hsapiens.UCSC.hg38)

# since the test need randomizations we have to set a seed
set.seed(42)


upbed1 <- "results/manuscript_figures/CRE_enrichment_plot/upregulating_upto2kb_OEs_nonprom.bed"
upbed2 <- "results/manuscript_figures/CRE_enrichment_plot/upregulating_2upto5kb_OEs_nonprom.bed"
upbed3 <- "results/manuscript_figures/CRE_enrichment_plot/upregulating_5upto10kb_OEs_nonprom.bed"
upbed4 <- "results/manuscript_figures/CRE_enrichment_plot/upregulating_10upto20kb_OEs_nonprom.bed"
upbed5 <- "results/manuscript_figures/CRE_enrichment_plot/upregulating_over20kb_OEs_nonprom.bed"
up_name <- "Enhancer"
down_name <- "Silencer"
mid_name <- "No Effect"
midbed1 <- "results/manuscript_figures/CRE_enrichment_plot/no_effect_upto2kb_OEs_nonprom.bed"
midbed2 <- "results/manuscript_figures/CRE_enrichment_plot/no_effect_2upto5kb_OEs_nonprom.bed"
midbed3 <- "results/manuscript_figures/CRE_enrichment_plot/no_effect_5upto10kb_OEs_nonprom.bed"
midbed4 <- "results/manuscript_figures/CRE_enrichment_plot/no_effect_10upto20kb_OEs_nonprom.bed"
midbed5 <- "results/manuscript_figures/CRE_enrichment_plot/no_effect_over20kb_OEs_nonprom.bed"
downbed1 <- "results/manuscript_figures/CRE_enrichment_plot/downregulating_upto2kb_OEs_nonprom.bed"
downbed2 <- "results/manuscript_figures/CRE_enrichment_plot/downregulating_2upto5kb_OEs_nonprom.bed"
downbed3 <- "results/manuscript_figures/CRE_enrichment_plot/downregulating_5upto10kb_OEs_nonprom.bed"
downbed4 <- "results/manuscript_figures/CRE_enrichment_plot/downregulating_10upto20kb_OEs_nonprom.bed"
downbed5 <- "results/manuscript_figures/CRE_enrichment_plot/downregulating_over20kb_OEs_nonprom.bed"
rds <- "results/manuscript_figures/figure3/stratified_distances_H3K4me3.rds"


data_up1 <- read.table(upbed1, header = F, sep = "\t", col.names = c("chr", "start", "end", "name", "strand", "score"))
data_down1 <- read.table(downbed1, header = F, sep = "\t", col.names = c("chr", "start", "end", "name", "strand", "score"))
data_mid1 <- read.table(midbed1, header = F, sep="\t", col.names = c("chr", "start", "end", "name", "strand", "score"))
data_up2 <- read.table(upbed2, header = F, sep = "\t", col.names = c("chr", "start", "end", "name", "strand", "score"))
data_down2 <- read.table(downbed2, header = F, sep = "\t", col.names = c("chr", "start", "end", "name", "strand", "score"))
data_mid2 <- read.table(midbed2, header = F, sep="\t", col.names = c("chr", "start", "end", "name", "strand", "score"))
data_up3 <- read.table(upbed3, header = F, sep = "\t", col.names = c("chr", "start", "end", "name", "strand", "score"))
data_down3 <- read.table(downbed3, header = F, sep = "\t", col.names = c("chr", "start", "end", "name", "strand", "score"))
data_mid3 <- read.table(midbed3, header = F, sep="\t", col.names = c("chr", "start", "end", "name", "strand", "score"))
data_up4 <- read.table(upbed4, header = F, sep = "\t", col.names = c("chr", "start", "end", "name", "strand", "score"))
data_down4 <- read.table(downbed4, header = F, sep = "\t", col.names = c("chr", "start", "end", "name", "strand", "score"))
data_mid4 <- read.table(midbed4, header = F, sep="\t", col.names = c("chr", "start", "end", "name", "strand", "score"))
data_up5 <- read.table(upbed5, header = F, sep = "\t", col.names = c("chr", "start", "end", "name", "strand", "score"))
data_down5 <- read.table(downbed5, header = F, sep = "\t", col.names = c("chr", "start", "end", "name", "strand", "score"))
data_mid5 <- read.table(midbed5, header = F, sep="\t", col.names = c("chr", "start", "end", "name", "strand", "score"))  

OE_UP_1 <- GenomicRanges::makeGRangesFromDataFrame(data_up1, seqnames.field = "chr",start.field = "start",end.field = "end")
OE_DOWN_1 <- GenomicRanges::makeGRangesFromDataFrame(data_down1, seqnames.field = "chr",start.field = "start",end.field = "end")
OE_MID_1 <- GenomicRanges::makeGRangesFromDataFrame(data_mid1, seqnames.field = "chr",start.field = "start",end.field = "end")
OE_UP_2 <- GenomicRanges::makeGRangesFromDataFrame(data_up2, seqnames.field = "chr",start.field = "start",end.field = "end")
OE_DOWN_2 <- GenomicRanges::makeGRangesFromDataFrame(data_down2, seqnames.field = "chr",start.field = "start",end.field = "end")
OE_MID_2 <- GenomicRanges::makeGRangesFromDataFrame(data_mid2, seqnames.field = "chr",start.field = "start",end.field = "end")
OE_UP_3 <- GenomicRanges::makeGRangesFromDataFrame(data_up3, seqnames.field = "chr",start.field = "start",end.field = "end")
OE_DOWN_3 <- GenomicRanges::makeGRangesFromDataFrame(data_down3, seqnames.field = "chr",start.field = "start",end.field = "end")
OE_MID_3 <- GenomicRanges::makeGRangesFromDataFrame(data_mid3, seqnames.field = "chr",start.field = "start",end.field = "end")
OE_UP_4 <- GenomicRanges::makeGRangesFromDataFrame(data_up4, seqnames.field = "chr",start.field = "start",end.field = "end")
OE_DOWN_4 <- GenomicRanges::makeGRangesFromDataFrame(data_down4, seqnames.field = "chr",start.field = "start",end.field = "end")
OE_MID_4 <- GenomicRanges::makeGRangesFromDataFrame(data_mid4, seqnames.field = "chr",start.field = "start",end.field = "end")
OE_UP_5 <- GenomicRanges::makeGRangesFromDataFrame(data_up5, seqnames.field = "chr",start.field = "start",end.field = "end")
OE_DOWN_5 <- GenomicRanges::makeGRangesFromDataFrame(data_down5, seqnames.field = "chr",start.field = "start",end.field = "end")
OE_MID_5 <- GenomicRanges::makeGRangesFromDataFrame(data_mid5, seqnames.field = "chr",start.field = "start",end.field = "end")

GR_OEs <- list(Upregulating_1=OE_UP_1,Downregulating_1=OE_DOWN_1, No_Effect_1=OE_MID_1,
			   Upregulating_2=OE_UP_2,Downregulating_2=OE_DOWN_2, No_Effect_2=OE_MID_2,
			   Upregulating_3=OE_UP_3,Downregulating_3=OE_DOWN_3, No_Effect_3=OE_MID_3,
			   Upregulating_4=OE_UP_4,Downregulating_4=OE_DOWN_4, No_Effect_4=OE_MID_4,
			   Upregulating_5=OE_UP_5,Downregulating_5=OE_DOWN_5, No_Effect_5=OE_MID_5)



# Now we prepare the data of ChIP-seq
# we do the same, read the bed files and turn them into GenomicRanges
features <- c("data/HepG2_H3K4me3_ENCSR575RRX/ENCFF982DUT.bed.gz")
names(features) <- c("H3K4me3")
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

set.seed(42)
lichic_chip_OE <-  crosswisePermTest(
  Alist = GR_OEs,
  Blist = ftlist,
  sampling = TRUE,
  fraction = 0.5, # approximately 1000 regions
  min_sampling = 2000,
  ranFUN = "randomizeRegions", 
  evFUN = "numOverlaps", 
  ntimes = 1000,
  genome = "hg38",
  mc.cores = 20, 
  mc.set.seed = FALSE
)
# Use clusterize=F when using only two sets of regions
# Finally we generate the matrix that summarise this analysis
lichic_chip_OE <- makeCrosswiseMatrix(lichic_chip_OE, clusterize=F, transform=T) #, transform=T) #, clusterize=T)
saveRDS(lichic_chip_OE, file=rds)
lichic_chip_OE <- readRDS(rds)

z_scores <- as.data.frame(lichic_chip_OE@matrix$GMat)
z_scores$pvalue <- lichic_chip_OE@matrix$GMat_pv
write.table(z_scores, file="results/manuscript_figures/revisions/stratified_distances_H3K4me3_zscores.tsv", sep="\t")
