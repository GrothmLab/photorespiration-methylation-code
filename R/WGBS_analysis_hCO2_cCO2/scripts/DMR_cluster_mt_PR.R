### This scrips plots sample clustering over gbM- and CMT2 target-enriched DMR clusters between genotypes and CO2 levels (Hankover et al. (2026) Nature PLants).

#####
### Set working directory to directory containing the R script
#setwd("/path/to/photorespiration_methylation_code/R/WGBS_analysis_hCO2_cCO2/scripts")
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

library(methylKit)
library(genomation)
library(ggplot2)
library(pheatmap)
library(gridExtra)
library(readr)
library(RColorBrewer)
library(circlize)
library(reshape2)
library(tidyverse)
library(dendextend)
library(svglite)

options(scipen=999)

samplesheet <- read_csv("../samplesheet_mt_PR.csv")
(sample.id <- as.list(unique(sort(samplesheet$name))))
(coldata <- read_csv("../WGBS_CO2_sample_info_mt_PR.csv"))

#Line colors for metaplots
clrplt8 <- c("indianred1", "indianred2", "indianred3", "indianred4", "steelblue1", "steelblue2", "steelblue3", "steelblue4")
clrplt4 <- c("indianred", "indianred4", "steelblue", "steelblue4")
clrplt4_mt <- c("indianred2", "peru", "turquoise4", "steelblue1")


#Line colors for heatmap
#col.ramp <- colorRampPalette(brewer.pal(9, "Reds") )(255)
colors_or <- brewer.pal(9,"OrRd")
colors_r <- brewer.pal(9,"Reds")
colors_rb <- rev(brewer.pal(11,"RdBu"))
colors_rg <- rev(brewer.pal(11,"RdGy"))
pal_or <- colorRampPalette(colors_or)
pal_r <- colorRampPalette(colors_r)
pal_rb <- colorRampPalette(colors_rb)
pal_rg <- colorRampPalette(colors_rg)


#Reposition brakes at the quantiles of the data
quantile_breaks <- function(xs, n = 100) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}


####################
#Generate methylkit object
context <- "CG"

#Read DMR regions from bed file as GRanges object
DMRs_PR <- readBed("../source-data/CO2_union_CG_DMRs2wtc_40_01.bed", zero.based = FALSE)
DMRs_mt <- readBed("../../WGBS_analysis_m_t_mt_wt/source-data/merged_CG_DMRs_40_01_union.bed", zero.based = FALSE)

DMRs_mt_PR <- c(DMRs_PR,DMRs_mt)

#Merge ranges that are separated by 100 bp
DMRs.merged <- GenomicRanges::reduce(DMRs_mt_PR, drop.empty.ranges=TRUE, min.gapwidth=101)

(file.list_PR <- as.list(list.files(path = "../source-data/methylkit_in", pattern = "CG_report_sorted\\.txt$", full.names = TRUE)))
(file.list_mt <- as.list(list.files(path = "../../WGBS_analysis_m_t_mt_wt/source-data/methylkit_in", pattern = "CG_report_sorted\\.txt$", full.names = TRUE)))

(file.list <- c(file.list_PR,file.list_mt))

#Generate methylkit object for all samples
myobj.s <- methRead(file.list,
                    sample.id[c(3:6,13:16,1:2,7:12)],
                    treatment=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,0,0),
                    pipeline="bismarkCytosineReport",
                    header=FALSE,
                    assembly="TAIR10",
                    context=context,
                    mincov=4,
                    dbtype=NA
)

#Subset tiled counts to merged DMRs
myobj.s.DMRs.merged <- regionCounts(myobj.s, DMRs.merged, cov.bases=4, mc.cores=2, save.db=FALSE)
#myobj.s.DMRs.merged.filt <- filterByCoverage(myobj.s.DMRs.merged, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.95)

#Select windows covered in all samples
meth.s.DMRs.merged <- methylKit::unite(myobj.s.DMRs.merged, destrand=FALSE, mc.cores=2)

#Calculate CG methylation over DMRs
meth.s.DMRs.merged.mat <- percMethylation(meth.s.DMRs.merged, rowids=TRUE)

#Pool replicates
pooled.meth.s.DMRs.merged <- NULL
pooled.meth.s.DMRs.merged <- pool(meth.s.DMRs.merged,sample.ids = c("mc","mh","wtc","wth","m","mt","t","wt"))

#Calculate CG methylation over DMRs
pooled.meth.s.DMRs.merged.mat <- percMethylation(pooled.meth.s.DMRs.merged, rowids=TRUE)

#Transpose and scale CG methylation matrix
pooled.meth.s.DMRs.merged.mat.z <- t(pooled.meth.s.DMRs.merged.mat) %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()

#####
# Add DMR annotation

# Import bed file as GRanges object
Ac2_DMRs <- readBed("../../WGBS_analysis_m_t_mt_wt/source-data/CG_DMRs_cluster2.bed", zero.based = FALSE) %>%
  GenomicRanges::reduce(drop.empty.ranges=TRUE, min.gapwidth=101)
myobj.s.Ac2_DMRs.merged.merged <- regionCounts(myobj.s, Ac2_DMRs, cov.bases=4, mc.cores=2, save.db=FALSE)
#myobj.s.DMRs.merged.filt <- filterByCoverage(myobj.s.DMRs.merged, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.95)

meth.s.Ac2_DMRs.merged.merged <- methylKit::unite(myobj.s.Ac2_DMRs.merged.merged, destrand=FALSE, mc.cores=2)
#meth.s.mh_mc_hyperDMRs.merged.mat <- percMethylation(meth.s.mh_mc_hyperDMRs.merged, rowids=TRUE)
pooled.meth.s.Ac2_DMRs.merged.merged <- pool(meth.s.Ac2_DMRs.merged.merged,sample.ids = c("mc","mh","wtc","wth","m","mt","t","wt"))
pooled.meth.s.Ac2_DMRs.merged.merged.mat <- percMethylation(pooled.meth.s.Ac2_DMRs.merged.merged, rowids=TRUE)

Bc3_DMRs <- readBed("../source-data/CO2_CG_cluster_3_DMRs.bed", zero.based = FALSE) %>%
  GenomicRanges::reduce(drop.empty.ranges=TRUE, min.gapwidth=101)
myobj.s.Bc3_DMRs.merged.merged <- regionCounts(myobj.s, Bc3_DMRs, cov.bases=4, mc.cores=2, save.db=FALSE)
#myobj.s.DMRs.merged.filt <- filterByCoverage(myobj.s.DMRs.merged, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.95)

meth.s.Bc3_DMRs.merged.merged <- methylKit::unite(myobj.s.Bc3_DMRs.merged.merged, destrand=FALSE, mc.cores=2)
#meth.s.Bc3_DMRs.merged.merged.mat <- percMethylation(meth.s.Bc3_DMRs.merged.merged, rowids=TRUE)
pooled.meth.s.Bc3_DMRs.merged.merged <- pool(meth.s.Bc3_DMRs.merged.merged,sample.ids = c("mc","mh","wtc","wth","m","mt","t","wt"))
pooled.meth.s.Bc3_DMRs.merged.merged.mat <- percMethylation(pooled.meth.s.Bc3_DMRs.merged.merged, rowids=TRUE)

mat_list <- list(pooled.meth.s.Ac2_DMRs.merged.merged.mat, pooled.meth.s.Bc3_DMRs.merged.merged.mat)
names(mat_list) <- c("Ac2","Bc3")

#####
# Plots fpr Ac2 DMRs

plot_title <- "Ac2 DMRs"

#Remove m and wt (redundant)
pooled.meth.s.Ac2_DMRs.merged.merged.mat.2 <- pooled.meth.s.Ac2_DMRs.merged.merged.mat[,-c(5,8)]

#Euclidean distance between rows (DMRs)
meth_dist <- dist(pooled.meth.s.Ac2_DMRs.merged.merged.mat.2)

# Pairwise correlation between rows (DMRs)
#meth_cor <- cor(t(pooled.meth.s.Ac2_DMRs.merged.merged.mat.2), use = "pairwise.complete.obs", method = "pearson")
#meth_dist <- as.dist(1-meth_cor)

meth_hclust <- hclust(meth_dist, method="ward.D2")
plot(meth_hclust, labels = FALSE)
abline(h = 2, col = "red", lwd = 1)
#abline(h = 8.2, col = "orange", lwd = 1)

sample_dist <- dist( t( pooled.meth.s.Ac2_DMRs.merged.merged.mat.2 ) )
sample_hclust <- hclust(sample_dist, method="complete")
plot(sample_hclust)
sample_hclust <- rotate(sample_hclust, order = c(5,6,4,1,3,2))
plot(sample_hclust)

dend <- as.dendrogram(sample_hclust)

# Customize labels
dend <- dend %>%
  set("labels_colors", c("indianred","indianred4","peru","turquoise4","steelblue4","steelblue")) %>%  # Set label colors
  set("labels_cex", 1.2)  # Set font size  # Set bold font

# Plot the dendrogram

#svglite(paste("../out/Complete_linkage_cluster_", context, "m_Ac2_DMRs_2.svg", sep = "_"), width=2.7, height=2.5)
plot(dend)
#dev.off()

# Summarise ntd counts 
meth_mean <- as.data.frame(pooled.meth.s.Ac2_DMRs.merged.merged.mat.2)
meth_mean$pos <- rownames(meth_mean)

# convert to long format
meth_mean <- pivot_longer(meth_mean, cols = 1:6, names_to = "sample_name", values_to = "meth")  %>% 
  # for each gene and condition
  group_by(pos, sample_name) %>%
  # calculate the mean (scaled) meth
  summarise(mean_meth = mean(meth),
            nrep = n()) %>% 
  ungroup()

meth_mean$sample_name <- factor(meth_mean$sample_name, levels = c("mc", "mh", "mt", "t", "wth", "wtc"))

m.cl <- meth_mean %>% 
  ggplot(aes(sample_name, mean_meth, fill = sample_name)) +
  geom_violin(width = 1, trim = TRUE) +
  scale_fill_manual(values = c("indianred","indianred4","peru","turquoise4","steelblue4","steelblue")) +
  scale_y_continuous(breaks = seq(0, 100, by = 50)) +
  theme_bw() +
  ggtitle(plot_title) +
  ylab(paste0(context,"m in DMRs")) +
  theme(axis.title.x=element_blank(), axis.title.y=element_text(size = 13), axis.text = element_text(size = 12), plot.title = element_text(size = 14, hjust = 0.5), legend.position="none") +
  geom_boxplot(width = 0.2, notch = FALSE, outlier.shape = 19, outlier.size = 0.02, outlier.alpha = 0.8, fill = "snow1") +
  geom_line(stat = "summary", fun = "median", colour = "brown", linewidth = 1, 
            aes(group = 1), alpha = 0.3)


#svglite(paste("../out/Cluster_CO2", context, "m_Ac2_DMRs_2.svg", sep = "_"), width=2.5, height=2) #Fig. 4b
plot(m.cl)
#dev.off()

#####
#Repeat for Bc3 cluster

plot_title <- "Bc3 DMRs"

#Remove m and wt (redundant)
pooled.meth.s.Bc3_DMRs.merged.merged.mat.2 <- pooled.meth.s.Bc3_DMRs.merged.merged.mat[,-c(5,8)]

#Euclidean distance between rows (DMRs)
meth_dist <- dist(pooled.meth.s.Bc3_DMRs.merged.merged.mat.2)

# Pairwise correlation between rows (DMRs)
#meth_cor <- cor(t(pooled.meth.s.Bc3_DMRs.merged.merged.mat.2), use = "pairwise.complete.obs", method = "pearson")
#meth_dist <- as.dist(1-meth_cor)

meth_hclust <- hclust(meth_dist, method="ward.D2")
plot(meth_hclust, labels = FALSE)
abline(h = 2, col = "red", lwd = 1)
#abline(h = 8.2, col = "orange", lwd = 1)

sample_dist <- dist( t( pooled.meth.s.Bc3_DMRs.merged.merged.mat.2 ) )
sample_hclust <- hclust(sample_dist, method="complete")
plot(sample_hclust)
sample_hclust <- rotate(sample_hclust, order = c(5,6,3,4,2,1))
plot(sample_hclust)

dend <- as.dendrogram(sample_hclust)

# Customize labels
dend <- dend %>%
  set("labels_colors", c("indianred","indianred4","peru","turquoise4","steelblue4","steelblue")) %>%  # Set label colors
  set("labels_cex", 1.2)  # Set font size  # Set bold font

# Plot the dendrogram

#svglite(paste("../out/Complete_linkage_cluster_", context, "m_Bc3_DMRs_2.svg", sep = "_"), width=2.7, height=2.5)
plot(dend)
#dev.off()

# Summarise ntd counts 
meth_mean <- as.data.frame(pooled.meth.s.Bc3_DMRs.merged.merged.mat.2)
meth_mean$pos <- rownames(meth_mean)

# convert to long format
meth_mean <- pivot_longer(meth_mean, cols = 1:6, names_to = "sample_name", values_to = "meth")  %>% 
  # for each gene and condition
  group_by(pos, sample_name) %>%
  # calculate the mean (scaled) meth
  summarise(mean_meth = mean(meth),
            nrep = n()) %>% 
  ungroup()

meth_mean$sample_name <- factor(meth_mean$sample_name, levels = c("mc", "mh", "mt", "t", "wth", "wtc"))

m.cl <- meth_mean %>% 
  ggplot(aes(sample_name, mean_meth, fill = sample_name)) +
  geom_violin(width = 1, trim = TRUE) +
  scale_fill_manual(values = c("indianred","indianred4","peru","turquoise4","steelblue4","steelblue")) +
  scale_y_continuous(breaks = seq(0, 100, by = 50)) +
  theme_bw() +
  ggtitle(plot_title) +
  ylab(paste0(context,"m in DMRs")) +
  theme(axis.title.x=element_blank(), axis.title.y=element_text(size = 13), axis.text = element_text(size = 12), plot.title = element_text(size = 14, hjust = 0.5), legend.position="none") +
  geom_boxplot(width = 0.2, notch = FALSE, outlier.shape = 19, outlier.size = 0.02, outlier.alpha = 0.8, fill = "snow1") +
  geom_line(stat = "summary", fun = "median", colour = "brown", linewidth = 1, 
            aes(group = 1), alpha = 0.3)


#svglite(paste("../out/Cluster_CO2", context, "m_Bc3_DMRs_2.svg", sep = "_"), width=2.5, height=2) #Fig. 4b
plot(m.cl)
#dev.off()

#####
#Clean up
rm(list=ls(pattern="^meth.s"))
rm(list=ls(pattern="^myobj.s"))
rm(list=ls(pattern="^pooled.meth.s"))
gc()
