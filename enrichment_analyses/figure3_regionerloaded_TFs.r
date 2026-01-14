suppressPackageStartupMessages(library(dplyr))
library(ggplot2)
library(tidyr)
suppressPackageStartupMessages(library(regioneReloaded))
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
  outfile <- "results/manuscript_figures/figure3/dotplot_TF_enrichment.pdf"
  rds <- "results/manuscript_figures/figure3/lichic_chip_all_distal_OE.rds"
} else {
  upbed <- args[1]
  downbed <- args[2]
  experiment <- args[3]
  up_name <- args[4]
  down_name <- args[5]
  outfile <- sprintf("results/manuscript_figures/figure3/dotplot_TF_enrichment_% s.pdf", experiment)
  rds <- sprintf("results/manuscript_figures/figure3/lichic_TF_chip_OE_% s.rds", experiment)
}

data_up <- read.table(upbed, header = F, sep = "\t", col.names = c("chr", "start", "end", "name", "strand", "score"))
data_down <- read.table(downbed, header = F, sep = "\t", col.names = c("chr", "start", "end", "name", "strand", "score"))

if (exists("midbed")) {
	data_mid <- read.table(midbed, header = F, sep = "\t", col.names = c("chr", "start", "end", "name", "strand", "score"))
	OE_GR_3 <- GenomicRanges::makeGRangesFromDataFrame(data_mid, seqnames.field = "chr",start.field = "start",end.field = "end")
}

OE_GR_1 <- GenomicRanges::makeGRangesFromDataFrame(data_up, seqnames.field = "chr",start.field = "start",end.field = "end")
OE_GR_2 <- GenomicRanges::makeGRangesFromDataFrame(data_down, seqnames.field = "chr",start.field = "start",end.field = "end")

if (exists("OE_GR_3")) {
  GR_OEs <- list(Upregulating_OE=OE_GR_1,Downregulating_OE=OE_GR_2, No_Effect_OE=OE_GR_3)	
} else {
  GR_OEs <- list(Upregulating_OE=OE_GR_1,Downregulating_OE=OE_GR_2)
}

df <- read.table("data/HepG2_ReMAP_TF_ChIPSeq/remap2022_Hep-G2_4cols.bed.gz")
TFs <- GenomicRanges::makeGRangesFromDataFrame(df,seqnames.field = "V1",start.field = "V2",end.field = "V3", keep.extra.columns = T)

split_granges <- split(TFs, GenomicRanges::mcols(TFs)$V4)
split_granges <- lapply(split_granges, function(gr) {
  GenomicRanges::mcols(gr)$V4 <- NULL
  return(gr)
})
ftlist <- c(split_granges)
##specific_experiments <- c("RAD21", "CTCF", "HNF1A", "ONECUT1", "EZH2", "SUZ12", "RNF2")
# Comment these to lines for unbiased analysis using all TFs
#specific_experiments <- c("RAD21", "CTCF", "YY1", "ZBTB21", "SMC1A", "SMC3", "STAG1", "NIPBL", "FOXA1")
#ftlist <- ftlist[specific_experiments]

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

# Select top 10 other enriched experiments (not in specific_experiments) for each dataset
top5_up <- up_data %>%
  #filter(adj.p_value < 0.05 & is.finite(norm_zscore)) %>%
  arrange(desc(norm_zscore)) %>%
  slice_head(n = 15) %>% 
  pull(name)

top5_down <- down_data %>%
  #filter(adj.p_value < 0.05 & is.finite(norm_zscore)) %>%
  arrange(desc(norm_zscore)) %>%
  slice_head(n = 15) %>% 
  pull(name)

# Filter for specific experiments
filtered_up <- bind_rows(up_data, mid_data) %>%
  filter(name %in% c(top5_up))
filtered_down <- bind_rows(down_data, mid_data) %>%
  filter(name %in% c(top5_down))

# Add a significance column to the data
filtered_up <- filtered_up %>%
  mutate(significance = ifelse(adj.p_value < 0.05, "Significant", "Not Significant"))
filtered_down <- filtered_down %>%
  mutate(significance = ifelse(adj.p_value < 0.05, "Significant", "Not Significant"))

# Combine all data
combined_data <- bind_rows(filtered_up, filtered_down)


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


# # Plot the data
# pdf(outfile, width = 8, height = 6)  # Adjust width and height
# ggplot(combined_data, aes(x = norm_zscore, y = reorder(name, norm_zscore))) +
#   geom_point(
#     aes(
#       size = -log10(adj.p_value),
#       color = significance  # Map significance to color
#     ),
#     shape = 21, alpha = 0.6, stroke = 1  # Uniform stroke for all points
#   ) +
#   labs(
#     x = "Enrichment (normalized z-score)",
#     y = "Experiment",
#     title = "Enrichment of CREs"
#   ) +
#   theme_minimal(base_size = 12) +
#   theme(
#     legend.position = "top",
# 	legend.justification = "left",
# 	legend.spacing.y = unit(0, "mm"),
#     legend.box = "vertical",  # Stack legend items vertically
#     #legend.margin = margin(t = 2, b = 2),  # Add margin around the legend
#     legend.text = element_text(size = 10, hjust = 0),  # Adjust legend text size
#     legend.title = element_text(size = 12),  # Adjust legend title size
#     panel.grid.major.y = element_line(color = "gray90"),
#     panel.grid.major.x = element_blank(),
#     axis.text.y = element_text(size = 10), #, margin = margin(r = 5)),  # Reduce y-axis text size for density
#     axis.text.x = element_text(size = 10)  # Adjust x-axis text size
#   ) +
#   scale_fill_manual(
#     values = c("Upregulating" = "#33a02c", "Downregulating" = "#e31a1c"),  # Green for Enhancer, Red for Silencing
#     labels = c("Upregulating" = up_name, "Downregulating" = down_name),  # Rename labels
#     name = "Regulation"
#   ) +
#   scale_color_manual(
#     values = c("Significant" = "black", "Not Significant" = "gray"),  # Black for significant, white for non-significant
#     name = "Significance"  # Add a legend title for significance
#   ) +
#   scale_size_continuous(
#     range = c(2, 5),  # Adjust dot size range for better visibility
#     name = "-log10 p-value"
#   ) +
#   guides(
#     fill = guide_legend(override.aes = list(size = 4, shape = 21, stroke = 0)),  # No outline for fill legend
#     color = guide_legend(override.aes = list(size = 4, shape = 21, stroke = 1))  # Black outline for significant points
#   ) +
#   coord_cartesian(clip = "off")  # Ensure no points are cut off
# dev.off()

# wide_up <- filtered_up %>%
# 	select(name, regulation, norm_zscore, adj.p_value, significance) %>%
# 	tidyr::pivot_wider(
# 		names_from = regulation,
# 		values_from = c(norm_zscore, adj.p_value, significance)
# 	) 

# wide_down <- filtered_down %>%
# 	select(name, regulation, norm_zscore, adj.p_value, significance) %>%
# 	tidyr::pivot_wider(
# 		names_from = regulation,
# 		values_from = c(norm_zscore, adj.p_value, significance)
# 	) 

# x_limits <- range(c(filtered_up$norm_zscore, filtered_down$norm_zscore), na.rm = TRUE)

# # Plot the data
# pdf("results/manuscript_figures/figure3/dotplot_enhancer_TF_enrichment.pdf", width = 6, height = 5.5)  # Adjust width and height
# ggplot(wide_up, aes(y = reorder(name, norm_zscore_Upregulating))) +
# 	geom_segment(
# 		aes(
# 			x = norm_zscore_Upregulating,
# 			xend = norm_zscore_NoEffect,
# 			yend = name
# 		),
# 		color = "gray70",
# 		size = 0.7
# 	) +
# 	# Upregulating end
# 	geom_point(
# 		aes(
# 			x = norm_zscore_Upregulating,
# 			shape = significance_Upregulating,
# 			fill = "Upregulating",
# 		),
# 		color = "#F4A261",
# 		stroke = 1,
# 		alpha = 0.8,
# 		size = 4, 
# 		show.legend = TRUE
# 	) +
# 	# No Effect end
# 	(if ("norm_zscore_NoEffect" %in% names(wide_up)) {
# 		geom_point(
# 			aes(
# 				x = norm_zscore_NoEffect,
# 				shape = significance_NoEffect,
# 				fill = "No Effect"
# 			),
# 			color = "#A8DADC",
# 			stroke = 1,
# 			alpha = 0.8,
# 			size = 4, 
# 			show.legend = TRUE
# 		)
# 	} else {
# 		NULL
# 	}) +
#   labs(
#     x = "Enrichment (normalized z-score)",
#     y = "Experiment",
#     title = "Enrichment of CREs"
#   ) +
#   theme_minimal(base_size = 12) +
# 	theme(
# 		legend.position = "top",
# 		legend.justification = "left",
# 		legend.box = "vertical",
# 		legend.box.just = "left",
# 		legend.spacing.y = unit(0, "mm"),
# 		legend.text = element_text(size = 12, hjust = 0),
# 		legend.title = element_text(size = 12),
# 		panel.grid.major.y = element_line(color = "gray90"),
# 		panel.grid.major.x = element_blank(),
# 		axis.text.y = element_text(size = 12),
# 		axis.text.x = element_text(size = 12)
# 	) +
#   scale_fill_manual(
#     values = c("Upregulating" = "#F4A261", "No Effect" = "#A8DADC"), #, "Downregulating" = "#e31a1c"),  # Green for Enhancer, Red for Silencing
#     labels = c("Upregulating" = "Enhancer", "No Effect" = "No Effect"), #, "Downregulating" = "Silencer"),  # Rename labels
#     name = "regulation"
#   ) +
#   scale_shape_manual(
#     values = c("Significant" = 8, "Not Significant" = 21),
#     limits = c("Significant", "Not Significant"),
#     drop = FALSE,
#     name = "",
#     na.value = 21
#   ) +
# 	guides(
# 		fill = guide_legend(override.aes = list(shape = 21, size = 4, stroke = 0)),
# 		shape = guide_legend(override.aes = list(size = 3.5, color= "black")),
# 		size = guide_legend(order = 3)
# 	) +
#   scale_x_continuous(limits = x_limits) + 
#   coord_cartesian(clip = "off")  # Ensure no points are cut off
# dev.off()

# # Plot the data
# pdf("results/manuscript_figures/figure3/dotplot_silencer_TF_enrichment.pdf", width = 6, height = 5.5)  # Adjust width and height
# ggplot(wide_down, aes(y = reorder(name, norm_zscore_Downregulating))) +
# 	geom_segment(
# 		aes(
# 			x = norm_zscore_Downregulating,
# 			xend = norm_zscore_NoEffect,
# 			yend = name
# 		),
# 		color = "gray70",
# 		size = 0.7
# 	) +
# 	# Downregulating end
# 	geom_point(
# 		aes(
# 			x = norm_zscore_Downregulating,
# 			shape = significance_Downregulating,
# 			fill = "Downregulating",
# 		),
# 		color = "#5D7CA8",
# 		stroke = 1,
# 		alpha = 0.8,
# 		size = 4, 
# 		show.legend = TRUE
# 	) +
# 	# No Effect end
# 	(if ("norm_zscore_NoEffect" %in% names(wide_up)) {
# 		geom_point(
# 			aes(
# 				x = norm_zscore_NoEffect,
# 				shape = significance_NoEffect,
# 				fill = "No Effect"
# 			),
# 			color = "#A8DADC",
# 			stroke = 1,
# 			alpha = 0.8,
# 			size = 4, 
# 			show.legend = TRUE
# 		)
# 	} else {
# 		NULL
# 	}) +
#   labs(
#     x = "Enrichment (normalized z-score)",
#     y = "Experiment",
#     title = "Enrichment of CREs"
#   ) +
#   theme_minimal(base_size = 12) +
# 	theme(
# 		legend.position = "top",
# 		legend.justification = "left",
# 		legend.box = "vertical",
# 		legend.box.just = "left",
# 		legend.spacing.y = unit(0, "mm"),
# 		legend.text = element_text(size = 12, hjust = 0),
# 		legend.title = element_text(size = 12),
# 		panel.grid.major.y = element_line(color = "gray90"),
# 		panel.grid.major.x = element_blank(),
# 		axis.text.y = element_text(size = 12),
# 		axis.text.x = element_text(size = 12)
# 	) +
#   scale_fill_manual(
#     values = c("Downregulating" = "#5D7CA8", "No Effect" = "#A8DADC"), #, "Downregulating" = "#e31a1c"),  # Green for Enhancer, Red for Silencing
#     labels = c("Downregulating" = "Silencer", "No Effect" = "No Effect"), #, "Downregulating" = "Silencer"),  # Rename labels
#     name = "regulation"
#   ) +
#   scale_shape_manual(
#     values = c("Significant" = 8, "Not Significant" = 21),
#     limits = c("Significant", "Not Significant"),
#     drop = FALSE,
#     name = "",
#     na.value = 21
#   ) +
# 	guides(
# 		fill = guide_legend(override.aes = list(shape = 21, size = 4, stroke = 0)),
# 		shape = guide_legend(override.aes = list(size = 3.5, color= "black")),
# 		size = guide_legend(order = 3)
# 	) +
#   scale_x_continuous(limits = x_limits) + 
#   coord_cartesian(clip = "off")  # Ensure no points are cut off
# dev.off()