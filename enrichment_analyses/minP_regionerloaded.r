library(dplyr)
library(ggplot2)
library(regioneReloaded)
library(BSgenome.Hsapiens.UCSC.hg38)

# since the test need randomizations we have to set a seed
set.seed(42)

args <- commandArgs(TRUE)
if (length(args) == 0) {
  stop("Provide either 'TF', 'unbiasedTF', or 'accessibility' as an argument.")
}
experiment = args[1]
if (experiment == "TF") {
  fn = "results/CMPRA5_minP/lichic_TF_proms_chip.rds"
  fig = "results/manuscript_figures/minP/minP_heatmap_TF_enrichment.pdf"
} else if (experiment == "unbiasedTF") {
  fn = "results/CMPRA5_minP/lichic_unbiasedTF_proms_chip.rds"
  fig = "results/manuscript_figures/minP/minP_dotplot_unbiasedTF_enrichment.pdf"
} else if (experiment == "accessibility") {
  fn = "results/CMPRA5_minP/lichic_CM_chip.rds"
  fig = "results/manuscript_figures/minP/minP_dotplot_accessibility_enrichment.pdf"
} else {
  stop("Invalid argument. Use 'TF', 'unbiasedTF', or 'accessibility'.")
} 

# upbed <- "results/CMPRA5_minP/higher_CMPRA_CREs.bed"
# downbed <- "results/CMPRA5_minP/lower_CMPRA_CREs.bed"
# backgroundbed <- "results/CMPRA5_minP/similar_CREs.bed"

# upbed <- "results/CMPRA5_minP/high_psi_promoters.bed"
# downbed <- "results/CMPRA5_minP/low_psi_promoters.bed"
# backgroundbed <- "results/CMPRA5_minP/mid_psi_promoters.bed"

# data_up <- read.table(upbed, header = F, sep = "\t", col.names = c("chr", "start", "end", "name", "strand", "score"))
# data_down <- read.table(downbed, header = F, sep = "\t", col.names = c("chr", "start", "end", "name", "strand", "score"))
# data_bg <- read.table(backgroundbed, header = F, sep="\t", col.names = c("chr", "start", "end", "name", "strand", "score"))

# OE_GR_1 <- GenomicRanges::makeGRangesFromDataFrame(data_up, seqnames.field = "chr",start.field = "start",end.field = "end")
# OE_GR_2 <- GenomicRanges::makeGRangesFromDataFrame(data_down, seqnames.field = "chr",start.field = "start",end.field = "end")
# OE_GR_3 <- GenomicRanges::makeGRangesFromDataFrame(data_bg, seqnames.field = "chr",start.field = "start",end.field = "end")
# GR_OEs <- list(high_PSI=OE_GR_1,low_PSI=OE_GR_2, mid_PSI=OE_GR_3)

# # Now we prepare the data of ChIP-seq
# # we do the same, read the bed files and turn them into GenomicRanges
# features <- c("data/HepG2_H3K4me3_ENCSR575RRX/ENCFF982DUT.bed.gz",
# 				"data/HepG2_H3K27ac_ENCSR000AMO/ENCFF392KDI.bed.gz",
# 				"data/HepG2_ATAC_ENCSR042AWH/ENCFF913MQB.bed.gz",
# 				# "data/HepG2_CTCF_ENCSR000BIE/ENCFF254DEQ.bed.gz",
# 				# "data/HepG2_RAD21_ENCSR000EEG/ENCFF145VFI.bed.gz",
# 				"data/HepG2_H3K9me3_ENCSR000ATD/ENCFF533JQH.bed.gz",
# 				"data/HepG2_H3K4me1_ENCSR000APV/ENCFF428FAW.bed.gz",
# 				"data/HepG2_DNase_ENCSR149XIL/ENCFF209DJG.bed.gz",
# 				"data/HepG2_H3K27me3_ENCSR000DUE/ENCFF982CSD.bed.gz", 
# 				"data/HepG2_H3K9ac_ENCSR000AMD/ENCFF358CJW.bed.gz",
# 				"data/HepG2_H3K79me2_ENCSR000AOM/ENCFF500KVU.bed.gz",
# 				"data/HepG2_H3K4me2_ENCSR000AMC/ENCFF769MUQ.bed.gz",
# 				"data/HepG2_H4K20me1_ENCSR000AMQ/ENCFF133PSF.bed.gz",
# 				"data/HepG2_H2AFZ_ENCSR000AOK/ENCFF474NDN.bed.gz",
# 				"data/HepG2_H3K36me3_ENCSR000DUD/ENCFF524PMT.bed.gz",
# 				"data/HepG2_ReMAP_TF_ChIPSeq/remap2022_Hep-G2_4cols.bed.gz")
# names(features) <- c("H3K4me3","H3K27ac", "ATAC", "H3K9me3", "H3K4me1", "DNase", "H3K27me3", "H3K9ac", "H3K79me2", "H3K4me2", "H4K20me1", "H2AFZ", "H3K36me3", "TFs")

# # then we have to transform them to GenomicRanges
# ftlist <- list()
# if (experiment == "accessibility") {
# 	for (i in 1:(length(features) - 1 ))
# 	{
# 	df <- read.table(features[i])
# 	ftlist[[i]] <- GenomicRanges::makeGRangesFromDataFrame(df,seqnames.field = "V1",start.field = "V2",end.field = "V3")
# 	}
# 	names(ftlist) <- names(features[1:(length(features) - 1)])
# } else {
# 	df <- read.table(features[length(features)])
# 	TFs <- GenomicRanges::makeGRangesFromDataFrame(df,seqnames.field = "V1",start.field = "V2",end.field = "V3", keep.extra.columns = T)

# 	split_granges <- split(TFs, GenomicRanges::mcols(TFs)$V4)
# 	split_granges <- lapply(split_granges, function(gr) {
# 		GenomicRanges::mcols(gr)$V4 <- NULL
# 		return(gr)
# 		})
# 	ftlist <- c(ftlist, split_granges)
# }

# if (experiment == "TF") {
# 	specific_experiments <- c("RAD21", "CTCF", "HNF1A", "ONECUT1", "EZH2", "SUZ12", "RNF2")
# 	ftlist <- ftlist[specific_experiments]
# }


# # now we performe the crosswise permutation test, using 5000 sampling of the data and 100 permutations. 
# # The function that we use for randomization if "randomizeRegions" (NOTE: It can take a while, but you 
# # # can do it in parallel using the argument mc.cores) and the evalutation function is "numOveralps".

# lichic_chip_OE <-  crosswisePermTest(
#   Alist = GR_OEs,
#   Blist = ftlist,
#   sampling = FALSE,
#   ranFUN = "randomizeRegions", 
#   evFUN = "numOverlaps", 
#   ntimes = 1000,
#   genome = "hg38",
#   mc.cores = 20
# )

# # Use clusterize=F when using only two sets of regions
# # Finally we generate the matrix that summarise this analysis
# lichic_chip_OE <- makeCrosswiseMatrix(lichic_chip_OE, clusterize=T, transform=T) #, transform=T) #, clusterize=T)


# saveRDS(lichic_chip_OE, file=fn)
# #saveRDS(lichic_chip_OE, file="results/manuscript_figures/figure3/lichic_chip_noeffectBG_OE.rds")


lichic_chip_OE <- readRDS(fn)


#pdf(fig)
# we can use the function to plot implemented in the package or take the matrix and plot it ourselves
# plot_enrich_prom <- plotCrosswiseMatrix(lichic_chip_OE)
# plot_enrich_prom

# Extract the two datasets
up_data <- lichic_chip_OE@multiOverlaps$high_PSI %>%
  mutate(PSI = "High")  # Add a column to label the dataset
down_data <- lichic_chip_OE@multiOverlaps$low_PSI %>%
  mutate(PSI = "Low")  # Add a column to label the dataset
mid_data <- lichic_chip_OE@multiOverlaps$mid_PSI %>%
  mutate(PSI = "Mid")  # Add a column to label the dataset

if (experiment == "accessibility" || experiment == "TF") {
	# Combine all data
	combined_data <- bind_rows(up_data, mid_data, down_data)
} else {
	# Select top 10 other enriched experiments (not in specific_experiments) for each dataset
	top5_up <- up_data %>%
	#filter(adj.p_value < 0.05 & is.finite(norm_zscore)) %>%
	arrange(desc(norm_zscore)) %>%
	slice_head(n = 10) %>% 
	pull(name)

	top5_down <- down_data %>%
	#filter(adj.p_value < 0.05 & is.finite(norm_zscore)) %>%
	arrange(desc(norm_zscore)) %>%
	slice_head(n = 10) %>% 
	pull(name)

	top5_mid <- mid_data %>%
	#filter(adj.p_value < 0.05 & is.finite(norm_zscore)) %>%
	arrange(desc(norm_zscore)) %>%
	slice_head(n = 10) %>% 
	pull(name)

	# Filter for specific experiments 
	filtered_up <- up_data %>%
	filter(name %in% c(top5_up, top5_mid, top5_down))
	filtered_mid <- mid_data %>%
	filter(name %in% c(top5_up, top5_mid, top5_down))
	filtered_down <- down_data %>%
	filter(name %in% c(top5_up, top5_mid, top5_down))

	# Combine all data
	combined_data <- bind_rows(filtered_up, filtered_mid, filtered_down)
}

# Add a significance column to the data
combined_data <- combined_data %>%
  mutate(significance = ifelse(adj.p_value < 0.05, "Significant", "Not Significant"))

# wide_data <- combined_data %>%
# 	select(name, PSI, norm_zscore, adj.p_value, significance) %>%
# 	tidyr::pivot_wider(
# 		names_from = PSI,
# 		values_from = c(norm_zscore, adj.p_value, significance)
# 	) %>%
# 	# compute an ordering value (average of the two z-scores) for plotting
# 	mutate(
# 		avg_z = rowMeans(
# 			cbind(norm_zscore_Upregulating, norm_zscore_Downregulating),
# 			na.rm = TRUE
# 		),
# 	) %>%
# 	arrange(avg_z)

# Plot the data
pdf(fig, width = 6, height = 5)  # Adjust width and height
ggplot(combined_data, aes(x = norm_zscore, y = reorder(name, norm_zscore), fill = PSI, color=PSI)) +
  geom_point(
    aes(
      shape = significance  # Map significance to shape
    ),
    alpha = 0.8, stroke = 1 , size=4 # Uniform stroke for all points
  ) +
  labs(
    x = "Enrichment (normalized z-score)",
    y = "Experiment",
    title = "Enrichment of CREs"
  ) +
  theme_minimal(base_size = 12) +
	theme(
		legend.position = "top",
		legend.justification = "left",
		legend.box = "vertical",
		legend.box.just = "left",
		legend.spacing.y = unit(0, "mm"),
		legend.text = element_text(size = 12, hjust = 0),
		legend.title = element_text(size = 12),
		panel.grid.major.y = element_line(color = "gray90"),
		panel.grid.major.x = element_blank(),
		axis.text.y = element_text(size = 12),
		axis.text.x = element_text(size = 12)
	) +
	# star/asterisk for significant, filled circle for non-significant
	scale_shape_manual(
		values = c("Significant" = 8, "Not Significant" = 21),
		name = ""
	)  +
  guides(
    fill = guide_legend(override.aes = list(size = 4, shape = 21, stroke = 0)),  # No outline for fill legend
    color = guide_legend(override.aes = list(size = 4, shape = 21, stroke = 0))  # Black outline for significant points
  ) +
  coord_cartesian(clip = "off")  # Ensure no points are cut off
dev.off()
