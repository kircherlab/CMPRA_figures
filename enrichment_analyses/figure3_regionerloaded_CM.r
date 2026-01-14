library(dplyr)
library(ggplot2)
library(tidyr)
library(regioneReloaded)
library(BSgenome.Hsapiens.UCSC.hg38)

# since the test need randomizations we have to set a seed
set.seed(42)

args <- commandArgs(TRUE)
if (length(args) == 0) {
  upbed <- "results/manuscript_figures/CRE_enrichment_plot/upregulating_distal_OEs_nonprom.bed"
  up_name <- "Enhancer"
  down_name <- "Silencer"
  mid_name <- "No Effect"
  midbed <- "results/manuscript_figures/CRE_enrichment_plot/no_effect_distal_OEs_nonprom.bed"
  downbed <- "results/manuscript_figures/CRE_enrichment_plot/downregulating_distal_OEs_nonprom.bed"
  outfile <- "results/manuscript_figures/figure3/dotplot_accessibility_all_distal_enrichment.pdf"
  rds <- "results/manuscript_figures/figure3/lichic_CM_chip_including_all_distal_OE.rds"
} else {
  upbed <- args[1]
  downbed <- args[2]
  experiment <- args[3]
  up_name <- args[4]
  down_name <- args[5]
  outfile <- sprintf("results/manuscript_figures/figure3/dotplot_accessibility_enrichment_% s.pdf", experiment)
  rds <- sprintf("results/manuscript_figures/figure3/lichic_CM_chip_OE_% s.rds", experiment)
}

data_up <- read.table(upbed, header = F, sep = "\t", col.names = c("chr", "start", "end", "name", "strand", "score"))
data_down <- read.table(downbed, header = F, sep = "\t", col.names = c("chr", "start", "end", "name", "strand", "score"))

if (exists("midbed")) {
  data_mid <- read.table(midbed, header = F, sep="\t", col.names = c("chr", "start", "end", "name", "strand", "score"))
  OE_GR_3 <- GenomicRanges::makeGRangesFromDataFrame(data_mid, seqnames.field = "chr",start.field = "start",end.field = "end")
}

OE_GR_1 <- GenomicRanges::makeGRangesFromDataFrame(data_up, seqnames.field = "chr",start.field = "start",end.field = "end")
OE_GR_2 <- GenomicRanges::makeGRangesFromDataFrame(data_down, seqnames.field = "chr",start.field = "start",end.field = "end")

if (exists("OE_GR_3")) {
  GR_OEs <- list(Upregulating_OE=OE_GR_1,Downregulating_OE=OE_GR_2, No_Effect_OE=OE_GR_3)
} else {
  GR_OEs <- list(Upregulating_OE=OE_GR_1,Downregulating_OE=OE_GR_2)
}


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


# Extract the two datasets
up_data <- lichic_chip_OE@multiOverlaps$Upregulating_OE %>%
  mutate(regulation = "Upregulating")  # Add a column to label the dataset
down_data <- lichic_chip_OE@multiOverlaps$Downregulating_OE %>%
  mutate(regulation = "Downregulating")  # Add a column to label the dataset

if ("No_Effect_OE" %in% names(lichic_chip_OE@multiOverlaps)) {
  mid_data <- lichic_chip_OE@multiOverlaps$No_Effect_OE %>%
    mutate(regulation = "NoEffect")
} else {
  mid_data <- NULL
}

# Combine all data
combined_data <- bind_rows(up_data, down_data, mid_data)

# Add a significance column to the data
combined_data <- combined_data %>%
  mutate(significance = ifelse(adj.p_value < 0.05, "Significant", "Not Significant"))

# Prepare wide-format data for a dumbbell plot
wide_data <- combined_data %>%
	select(name, regulation, norm_zscore, adj.p_value, significance) %>%
	pivot_wider(
		names_from = regulation,
		values_from = c(norm_zscore, adj.p_value, significance)
	) %>%
	# compute an ordering value (average of the two z-scores) for plotting
	mutate(
		avg_z = rowMeans(
			cbind(norm_zscore_Upregulating, norm_zscore_Downregulating),
			na.rm = TRUE
		),
	) %>%
	arrange(avg_z)

# Plot as dumbbell: horizontal segment per name with two end points.
pdf(outfile, width = 6, height = 5)
ggplot(wide_data, aes(y = reorder(name, avg_z))) +
	# connecting segment
	geom_segment(
		aes(
			x = norm_zscore_Upregulating,
			xend = norm_zscore_Downregulating,
			yend = name
		),
		color = "gray70",
		size = 0.7
	) +
	# Downregulating end
	geom_point(
		aes(
			x = norm_zscore_Downregulating,
			shape = significance_Downregulating,
			fill = "Downregulating"
		),
		color = "#5D7CA8",
		stroke = 1,
		alpha = 0.8,
		size = 4
	) +
	# Upregulating end
	geom_point(
		aes(
			x = norm_zscore_Upregulating,
			shape = significance_Upregulating,
			fill = "Upregulating",
		),
		color = "#F4A261",
		stroke = 1,
		alpha = 0.8,
		size = 4
	) +
	# No Effect end
	(if ("norm_zscore_NoEffect" %in% names(wide_data)) {
		geom_point(
			aes(
				x = norm_zscore_NoEffect,
				shape = significance_NoEffect,
				fill = "No Effect"
			),
			color = "#A8DADC",
			stroke = 1,
			alpha = 0.8,
			size = 4
		)
	} else {
		NULL
	}) +
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
	scale_fill_manual(
		values = if ("norm_zscore_NoEffect" %in% names(wide_data)) {
			c("Upregulating" = "#F4A261", "Downregulating" = "#5D7CA8", "No Effect" = "#A8DADC")
		} else {
			c("Upregulating" = "#F4A261", "Downregulating" = "#5D7CA8")
		},
		labels = if ("norm_zscore_NoEffect" %in% names(wide_data)) {
			c(
				"Upregulating" = up_name,
				"Downregulating" = down_name,
				"No Effect" = if (exists("mid_name")) mid_name else "No Effect"
			)
		} else {
			c("Upregulating" = up_name, "Downregulating" = down_name)
		},
		name = ""
	) +
	# star/asterisk for significant, filled circle for non-significant
	scale_shape_manual(
		values = c("Significant" = 8, "Not Significant" = 21),
		name = ""
	) +
	guides(
		fill = guide_legend(override.aes = list(shape = 21, size = 4, stroke = 0)),
		shape = guide_legend(override.aes = list(size = 3.5, color= "black")),
		size = guide_legend(order = 3)
	) +
	coord_cartesian(clip = "off") +
	scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)))
dev.off()
