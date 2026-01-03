### This scrips computes C/T counts over union sets of DMRs from all pairwise comparisons with mutants and wild type 
### and generates clustered heatmaps of DNA methylation ratios in DMRs, visualizes distributions of DNA methylation levels 
### per cluster, and DMR overlaps with genomic regions per cluster. The DMRs were computed using the standard methylkit
### procedure, as described in the Methods of Hankover et al. (2026) Nature Plants, and are provided as source-data.

#####
### Set working directory to directory containing the R script
#setwd("/path/to/photorespiration_methylation_code/R/WGBS_analysis_m_t_mt_wt/scripts")
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

samplesheet <- read_csv("../samplesheet.csv")
(sample.id <- as.list(sort(samplesheet$sample)))
(coldata <- read_csv("../WGBS_mt_sample_info.csv"))

#Line colors for metaplots
#clrplt <- c("indianred3", "indianred2", "peru", "navajowhite3", "turquoise4", "turquoise", "steelblue4", "steelblue")
clrplt4 <- c("indianred3", "peru", "turquoise4", "steelblue4")

#Line colors for heatmap
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

#Script fpr counting region overlaps 
source("./count_chromatinstate_overlaps_GR.R")

#####
###Sites for methylation stats
#GbM sites
GbM_sites <- readRDS("../source-data/GbM_sites.RDS") %>%
  GenomicRanges::reduce(drop.empty.ranges=TRUE, min.gapwidth=101)

#CMT2 sites
cmt2_sites <- readBed("../source-data/cmt2_forddcc_hypoCHH_DMR_vs3reps_min4filter_mg200_sort.bed", zero.based = FALSE) %>%
  GenomicRanges::reduce(drop.empty.ranges=TRUE, min.gapwidth=101)

#PolV sites
PolV_sites <- readBed("../source-data/Pol_V_Col_031715_nochr.bed", zero.based = FALSE) %>%
  GenomicRanges::reduce(drop.empty.ranges=TRUE, min.gapwidth=101)

#Hypervariable sites
HV_sites <- readBed("../source-data/merged_bins_sd0_2_CDMR_72_Ecker_Weigel_union_spontaneous_DMRs.bed") %>%
  GenomicRanges::reduce(drop.empty.ranges=TRUE, min.gapwidth=101)

#MTHFD1 sites
mthfd1_CG_sites <- readBed("../source-data/m_vs_wt_CG_hypoDMRs_40_01.bed", zero.based = FALSE) %>%
  GenomicRanges::reduce(drop.empty.ranges=TRUE, min.gapwidth=101)
mthfd1_CHG_sites <- readBed("../source-data/m_vs_wt_CHG_hypoDMRs_20_01.bed", zero.based = FALSE) %>%
  GenomicRanges::reduce(drop.empty.ranges=TRUE, min.gapwidth=101)
mthfd1_CHH_sites <- readBed("../source-data/m_vs_wt_CHH_hypoDMRs_10_01.bed", zero.based = FALSE) %>%
  GenomicRanges::reduce(drop.empty.ranges=TRUE, min.gapwidth=101)

#Whole genome bed for rest
other_sites <- readBed("../source-data/TAIR10_1to5.bed")

####################
context <- "CG"

#Read DMR regions from bed file as GRanges object
DMRs <- readBed("../source-data/merged_CG_DMRs_40_01_union.bed", zero.based = FALSE)

#Merge ranges that are separated by 100 bp
DMRs.merged <- GenomicRanges::reduce(DMRs, drop.empty.ranges=TRUE, min.gapwidth=101)

(file.list <- as.list(list.files(path = "../source-data/methylkit_in", pattern = "CG_report_sorted\\.txt$", full.names = TRUE)))

#Generate methylkit object for m, mt, t, and wt
myobj.s <- methRead(file.list[1:8],
                    as.list(sample.id[1:8]),
                    treatment=c(1,1,2,2,3,3,0,0),
                    pipeline="bismarkCytosineReport",
                    header=FALSE,
                    assembly="TAIR10",
                    context=context,
                    mincov=4,
                    dbtype=NA
)

#####
#Make chromosome plots
myobj.s.10kb <- NULL
myobj.s.10kb <- tileMethylCounts(myobj.s,win.size=10000,step.size=10000,cov.bases=10,mc.cores=2)
meth.s.10kb <- methylKit::unite(myobj.s.10kb, destrand=FALSE, mc.cores=2)
pooled.meth.s.10kb <- NULL
pooled.meth.s.10kb <- pool(meth.s.10kb, sample.ids = c("m","mt","t","wt"))
pooled.gr.10kb <- as(pooled.meth.s.10kb, "GRanges")
myobj.grl.10kb <- c()
for (i1 in seq(1,12,3)) {
  mcols(pooled.gr.10kb, level="within")[,"coverage"] <- mcols(pooled.gr.10kb, level="within")[,i1]
  mcols(pooled.gr.10kb, level="within")[,"numCs"] <- mcols(pooled.gr.10kb, level="within")[,i1+1]
  mcols(pooled.gr.10kb, level="within")[,"numTs"] <- mcols(pooled.gr.10kb, level="within")[,i1+2]
  myobj.grl.10kb <- c(myobj.grl.10kb, pooled.gr.10kb[,c("coverage","numCs","numTs")])
}
names(myobj.grl.10kb) <- c("m","mt","t","wt")

pc <- list()
pd <- list()
for (c in c("1","2","3","4","5")) {
  chr <- list()
  plot_title <- paste0("Chr ", c)
  dfcb <- tibble()
  for (i1 in 1:length(myobj.grl.10kb)){
    grl <- myobj.grl.10kb[[i1]]
    seqlevels(grl, pruning.mode="coarse") <- c
    grl$meth.lvl <- grl$numCs / grl$coverage
    methcol <- names(myobj.grl.10kb)[i1]
    chr[[i1]] <- data.frame(pos = start(grl)/1000000, meth.lvl = grl$meth.lvl)
    colnames(chr[[i1]])[[2]] <- methcol
    dfc <- melt(chr[[i1]], "pos", methcol, variable.name="sample")
    dfcb <- rbind(dfcb, dfc)
  }
  dfcb <- dfcb %>%
    group_by(pos) %>%
    #divide by wt value
    mutate(meth_diff = value - value[sample == "wt"]) %>%
    ungroup()
  plotc <- ggplot(dfcb)
  plotc <- plotc + theme_classic()
  plotc <- plotc + geom_smooth(aes_string(x='pos', y='value', col='sample'), method="loess", span=0.1, se=FALSE, linewidth=0.5)
  plotc <- plotc + scale_colour_manual(values=clrplt4)
  plotc <- plotc + ggtitle(plot_title)
  plotc <- plotc + theme(axis.line = element_blank(), axis.ticks =element_blank(), axis.title.x = element_blank(), axis.text = element_blank(), plot.title = element_text(size = 14, hjust = 0.5), legend.position = "none", 
                         axis.title.y = element_blank(), plot.margin = unit(c(2,0,0,0),"mm"))
  plotc <- plotc + coord_cartesian(ylim = c(0, 0.9), expand = TRUE, default = TRUE)
  plotc <- plotc + scale_y_continuous(breaks = seq(0, 0.8, by = 0.8))
  plotc <- plotc + scale_x_continuous(breaks = seq(0, max(dfcb$pos), by = 8))
  pc[[c]] <- plotc
  
  plotd <- ggplot(dfcb)
  plotd <- plotd + theme_classic()
  plotd <- plotd + geom_smooth(aes_string(x='pos', y='meth_diff', col='sample'), method="loess", span=0.1, se=FALSE, linewidth=0.5)
  plotd <- plotd + scale_colour_manual(values=clrplt4)
  plotd <- plotd + ggtitle(plot_title)
  plotd <- plotd + theme(axis.line.y = element_blank(), axis.ticks.y =element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_blank(), plot.title = element_blank(), legend.position = "none", 
                         axis.title.y = element_blank(), plot.margin = unit(c(0,0,2,0),"mm"))
  plotd <- plotd + coord_cartesian(ylim = c(-0.4, 0.04), expand = TRUE, default = TRUE)
  plotd <- plotd + scale_y_continuous(breaks = seq(-0.3, 0.0, by = 0.3))
  plotd <- plotd + scale_x_continuous(breaks = seq(0, max(dfcb$pos), by = 8))
  pd[[c]] <- plotd
}
pc[[1]] <- pc[[1]] + ylab(paste0(context,'m'))
pc[[1]] <- pc[[1]] + theme(axis.title.y = element_text(size = 13, angle = 90),
                           axis.text.y = element_text(size = 12, angle = 90, hjust = 0.2, margin = margin(r=2, unit="pt")),
                           axis.line.y = element_line(linewidth = 0.5), axis.ticks.y = element_line(linewidth = 0.5))
pd[[1]] <- pd[[1]] + ylab(~ Delta)
pd[[1]] <- pd[[1]] + theme(axis.title.y = element_text(size = 13, angle = 90, hjust = 0.3), 
                           axis.text.y = element_text(size = 12, angle = 90, hjust = 0.9, margin = margin(r=2, unit="pt")),
                           axis.line.y = element_line(linewidth = 0.5), axis.ticks.y = element_line(linewidth = 0.5))
pl <- pc[[5]]
pl <- pl + theme(legend.position = "top", legend.text = element_text(size = 12), legend.title = element_text(size = 16))

pl_grid <- grid.arrange(pc[[1]],pc[[2]],pc[[3]],pc[[4]],pc[[5]],pd[[1]],pd[[2]],pd[[3]],pd[[4]],pd[[5]], nrow=2, heights = c(1, 0.8), widths = c(4.5,2.3,2.5,2.2,2.8)) #adjusted plot size to 4.45 x 1.58 inches

#svglite(paste0("../out/", context,"m_chrom_all.svg"), width = 5, height = 1.7)
plot(pl_grid)
#dev.off()

#####
#Subset tiled counts to merged DMRs
myobj.s.DMRs.merged <- regionCounts(myobj.s, DMRs.merged, cov.bases=4, mc.cores=2, save.db=FALSE)
#myobj.s.DMRs.merged.filt <- filterByCoverage(myobj.s.DMRs.merged, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.95)

#Select windows covered in all samples
meth.s.DMRs.merged <- methylKit::unite(myobj.s.DMRs.merged, destrand=FALSE, mc.cores=2)

#####
#Calculate CG methylation over DMRs
meth.s.DMRs.merged.mat <- percMethylation(meth.s.DMRs.merged, rowids=TRUE)

#Pool replicates
pooled.meth.s.DMRs.merged <- NULL
pooled.meth.s.DMRs.merged <- pool(meth.s.DMRs.merged,sample.ids = c("m","mt","t","wt"))

#Calculate CG methylation over DMRs
pooled.meth.s.DMRs.merged.mat <- percMethylation(pooled.meth.s.DMRs.merged, rowids=TRUE)

#Transpose and scale CG methylation matrix
pooled.meth.s.DMRs.merged.mat.z <- t(pooled.meth.s.DMRs.merged.mat) %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()

#####
#Euclidean distance between rows (DMRs)
#meth_dist <- dist(pooled.meth.s.DMRs.merged.mat.z)

# Pairwise correlation between rows (DMRs)
meth_cor <- cor(t(pooled.meth.s.DMRs.merged.mat.z), use = "pairwise.complete.obs", method = "pearson")
meth_dist <- as.dist(1-meth_cor)

meth_hclust <- hclust(meth_dist, method="ward.D2")
plot(meth_hclust, labels = FALSE)
abline(h = 2, col = "red", lwd = 1)
#abline(h = 8.2, col = "orange", lwd = 1)

sample_dist <- dist( t( pooled.meth.s.DMRs.merged.mat.z ) )
sample_hclust <- hclust(sample_dist, method="complete")
plot(sample_hclust)

#####
#Make clustered heatmap and DNA methylation lots for DMRs

kmeans <- NA  
#if (length(DMRs) > 1000) {
#  kmeans <- 1000
#}

#Add clusters to row annotation
meth_cluster <- cutree(meth_hclust, k = 4) %>% 
  # turn the named vector into a tibble
  enframe(name = "pos", value = "cluster")


#Find features and genes next to DMRs
#####

sites_list <- list(HV_sites, GbM_sites, PolV_sites, cmt2_sites, mthfd1_CG_sites, other_sites)
names(sites_list) <- c("Hyper", "GbM", "RdDM", "CMT2", "Other m", "Rest")

clusters.features.ann <- annotateWithFeatures(DMRs.merged,GRangesList(sites_list))

#svglite(paste0("../out/", context, "_DMRs_features_overlaps_main.svg"), width=7.3, height=3.7)
genomation::plotTargetAnnotation(clusters.features.ann,precedence=TRUE,
                     col=brewer.pal(6,"Dark2"),
                     main="DMR features overlap")
#dev.off()

boxvp_list <- list()

for (i in 1:length(sites_list)) {
  myobj.s.sites <- regionCounts(myobj.s.DMRs.merged, sites_list[[i]], mc.cores=2, save.db=FALSE) %>%
    filterByCoverage(lo.count=4, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)
  if (sum(unlist(lapply(myobj.s.sites, nrow))) > 0) {
    meth.s.sites <- methylKit::unite(myobj.s.sites, min.per.group=1L, mc.cores=2)
    pooled.meth.s.sites <- pool(meth.s.sites,sample.ids = c("m","mt","t","wt"))
    sites_meth <- as.data.frame(percMethylation(pooled.meth.s.sites)) %>%
      pivot_longer(cols = 1:4, names_to = "samples", values_to = "meth") %>%
      add_column(group = rep(c("m","mt","t","wt"), nrow(meth.s.sites))) 
    
    #Boxplots for methylation over regions
    plot_title <- paste0(context, "m over \n", names(sites_list)[i], " sites")
    boxvp_list[[i]] <- ggplot(sites_meth, aes(x=samples, y=meth, fill=as.factor(samples))) + 
      geom_violin(width = 1, trim = TRUE) +
      scale_fill_manual(values = clrplt4) +
      coord_cartesian(ylim=c(0,100), expand = TRUE, default = TRUE) +
      theme_bw() +
      ggtitle(plot_title) +
      ylab(paste0(context,"m (%)")) +
      theme(axis.title.x=element_blank(), axis.title.y=element_text(size = 13), axis.text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(size = 14, hjust = 0.5), legend.position="none") +
      geom_boxplot(width = 0.2, notch = FALSE, outlier.shape = 19, outlier.size = 0.02, outlier.alpha = 0.8, fill = "snow1")
    
    # Unpaired Wilcoxon rank test (two-sided)
    kruskal_result <- kruskal.test(meth ~ group, data = sites_meth)
    wilcoxon_result <- pairwise.wilcox.test(sites_meth$meth, sites_meth$group, p.adjust.method = "BH")
    #outfile <- paste0("../out/", context, "m_over_", names(sites_list)[i], "_sites.csv")
    #sink(outfile)
    
    # Write Kruskal-Wallis results
    cat("Kruskal-Wallis Test Results\n")
    cat("---------------------------\n")
    cat(paste("Statistic:", kruskal_result$statistic, "\n"))
    cat(paste("P-Value:", kruskal_result$p.value, "\n\n"))
    
    # Write Wilcoxon results
    cat("Pairwise Wilcoxon Test Results\n")
    cat("------------------------------\n")
    print(wilcoxon_result)
    
    # Close the connection
    #sink()
  }
  else {
    boxvp_list[[i]] <- plot.new()
  }
}  
#svglite(paste0("../out/context, "m_over_sites_main.svg"), width=7, height=6)
grid.arrange(boxvp_list[[1]],boxvp_list[[2]],boxvp_list[[3]],boxvp_list[[4]], boxvp_list[[5]], boxvp_list[[6]], nrow=2)
#dev.off()

#Data frame list for overlaps of cluster DMRs with chromatin states
ol_df_list <- list()

for (c in 1:max(meth_cluster$cluster)){
  cluster.DMRs.ann <- dplyr::filter(meth_cluster, cluster == c) %>%
    separate_wider_delim(cols = pos, delim = ".", names = c("chr", "start", "end")) %>%
    as.data.frame() %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                             ignore.strand=TRUE)
  ol_df_list[[c]] <- countCSoverlapsGR(cluster.DMRs.ann)
  clusters.features.ann <- annotateWithFeatures(cluster.DMRs.ann,GRangesList(sites_list))
  #svglite(paste0("../out/context, "_DMRs_cluster_", c, "_features_overlaps_main.svg"), width=7.3, height=3.7)
  genomation::plotTargetAnnotation(clusters.features.ann,precedence=TRUE,
                       col=brewer.pal(8,"Dark2"),
                       main=paste0("Cluster ", c," features overlap"))
  #dev.off()
  
  boxvp_list <- list()
  myobj.s.DMRs.merged.cluster <- regionCounts(myobj.s.DMRs.merged, cluster.DMRs.ann, mc.cores=2, save.db=FALSE)
  for (i in 1:length(sites_list)) {
    myobj.s.sites <- regionCounts(myobj.s.DMRs.merged.cluster, sites_list[[i]], mc.cores=2, save.db=FALSE) %>%
      filterByCoverage(lo.count=4, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)
    if (sum(unlist(lapply(myobj.s.sites, nrow))) > 0) {
      meth.s.sites <- methylKit::unite(myobj.s.sites, min.per.group=1L, mc.cores=2)
      pooled.meth.s.sites <- pool(meth.s.sites,sample.ids = c("m","mt","t","wt"))
      sites_meth <- as.data.frame(percMethylation(pooled.meth.s.sites)) %>%
        pivot_longer(cols = 1:4, names_to = "samples", values_to = "meth") %>%
        add_column(group = rep(c("m","mt","t","wt"), nrow(meth.s.sites))) 
    
      plot_title <- paste0("Cluster ", c, " ", context, "m over \n", names(sites_list)[i], " sites")
      boxvp_list[[i]] <- ggplot(sites_meth, aes(x=samples, y=meth, fill=as.factor(samples))) + 
        geom_violin(width = 1, trim = TRUE) +
        scale_fill_manual(values = clrplt4) +
        coord_cartesian(ylim=c(0,100), expand = TRUE, default = TRUE) +
        theme_bw() +
        ggtitle(plot_title) +
        ylab(paste0(context,"m (%)")) +
        theme(axis.title.x=element_blank(), axis.title.y=element_text(size = 13), axis.text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(size = 14, hjust = 0.5), legend.position="none") +
        geom_boxplot(width = 0.2, notch = FALSE, outlier.shape = 19, outlier.size = 0.02, outlier.alpha = 0.8, fill = "snow1")
    
    }
    else {
      boxvp_list[[i]] <- plot.new()
    }
  }  
  #svglite(paste0("../out/Cluster_", c, "_", context, "m_over_sites_main.svg"), width=7, height=6)
  grid.arrange(boxvp_list[[1]],boxvp_list[[2]],boxvp_list[[3]],boxvp_list[[4]], boxvp_list[[5]], boxvp_list[[6]], nrow=2)
  #dev.off()
} 

#####

#Row annotation df
meth_cluster.df <- data.frame(cluster = factor(meth_cluster$cluster))
rownames(meth_cluster.df) <- meth_cluster$pos

#Col annotation df
sample_anno.df <- data.frame(sample = factor(colnames(pooled.meth.s.DMRs.merged.mat.z)))
rownames(sample_anno.df) <- colnames(pooled.meth.s.DMRs.merged.mat.z)

#Define annotation colors
ann_colors = list(
  cluster = brewer.pal(length(levels(meth_cluster.df$cluster)),("Set3")),
  sample = c(m="indianred", mt="peru", t="turquoise4", wt="steelblue4")
)
names(ann_colors$cluster) <- levels(meth_cluster.df$cluster)

#Plot chromatin states overlaps
#Add CS overlaps for all DMRs
ol_df_list[[c+1]] <- countCSoverlapsGR(DMRs.merged)
cs_ol_df <- Reduce(function(df1, df2) merge(df1, df2, by = "state", all = TRUE), ol_df_list)
colnames(cs_ol_df)[2:6] <- c("c1","c2","c3","c4","all")
#Calculate column sums
column_sums <- colSums(cs_ol_df[, -1])  # Exclude the 'state' column
#Convert values to percentages
cs_mat <- mapply('/', cs_ol_df[, -1], column_sums)
#Transpose and scale
cs_mat_z <- t(cs_mat[,1:4]) %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()
# Replace NA values with 0
cs_mat_z[is.na(cs_mat_z)] <- 0
#cs_mat_long <- pivot_longer(cs_mat, c1:c4, names_to = "column", values_to = "percentage")
rownames(cs_mat) <- cs_ol_df$state

plot_title <- "DMR-CS overlap\nper cluster"

#Reposition breaks
mat_breaks <- quantile_breaks(cs_mat_z, n = 101)

col_anno <- data.frame(cluster = 1:4)
rownames(col_anno) <- c("c1","c2","c3","c4")

pm_csz <- pheatmap::pheatmap(cs_mat_z,
                           color = pal_rg(length(mat_breaks) - 1),
                           show_rownames = FALSE,
                           show_colnames = TRUE,
                           fontsize = 16,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           treeheight_col = 16,
                           border_color = NA,
                           breaks=mat_breaks,
                           main=plot_title,
                           #annotation_row=meth_cluster.df,
                           #annotation_names_row=FALSE,
                           annotation_col=col_anno,
                           annotation_names_col=FALSE,
                           annotation_colors=ann_colors,
                           annotation_legend=TRUE
)

#All DMRs only
mat_d <- as.data.frame(cs_mat[,5])
colnames(mat_d) <- "all"
#Reposition breaks
mat_breaks <- quantile_breaks(mat_d$all, n = 101)

col_anno <- data.frame(cluster = "all")
rownames(col_anno) <- "all"
ann_colors = list(
  cluster = "white"
)

names(ann_colors$cluster) <- "all"

pm_csa <- pheatmap::pheatmap(mat_d,
                             color = pal_r(length(mat_breaks) - 1),
                             show_rownames = TRUE,
                             show_colnames = TRUE,
                             fontsize = 16,
                             cluster_rows = FALSE,
                             cluster_cols = FALSE,
                             treeheight_col = 16,
                             border_color = NA,
                             breaks=mat_breaks,
                             main=plot_title,
                             #annotation_row=meth_cluster.df,
                             #annotation_names_row=FALSE,
                             annotation_col=col_anno,
                             annotation_names_col=FALSE,
                             annotation_colors=ann_colors,
                             annotation_legend=FALSE
)

#svglite(paste("../out/Heatmap_", context, "DMRs_CS_overlaps.svg", sep = "_"), width=5, height=8)
grid.arrange(pm_csa[[4]],pm_csz[[4]],nrow=1, widths=c(1,2.1))
#dev.off()

#####
#Heatmaps of DMR methylation

plot_title <- paste0(context, "m over ",context," DMRs")

#Row annotation df
meth_cluster.df <- data.frame(cluster = factor(meth_cluster$cluster))
rownames(meth_cluster.df) <- meth_cluster$pos

#Define annotation colors
ann_colors = list(
  cluster = brewer.pal(length(levels(meth_cluster.df$cluster)),("Set3")),
  sample = c(m="indianred", mt="peru", t="turquoise4", wt="steelblue4")
)
names(ann_colors$cluster) <- levels(meth_cluster.df$cluster)

#Reposition breaks
mat_breaks <- quantile_breaks(pooled.meth.s.DMRs.merged.mat.z, n = 101)

ph.z <- pheatmap::pheatmap(pooled.meth.s.DMRs.merged.mat.z,
                           kmeans_k = kmeans,
                           color = pal_rb(length(mat_breaks) - 1),
                           show_rownames = FALSE,
                           show_colnames = FALSE,
                           fontsize = 24,
                           cutree_rows = 4,
                           cluster_rows = meth_hclust,
                           cluster_cols = sample_hclust,
                           treeheight_row = 50,
                           treeheight_col = 24,
                           border_color = NA,
                           breaks=mat_breaks,
                           main=plot_title,
                           annotation_row=meth_cluster.df,
                           annotation_names_row=FALSE,
                           annotation_col=sample_anno.df,
                           annotation_names_col=FALSE,
                           annotation_colors=ann_colors,
                           annotation_legend=FALSE
)
dev.off()

#svglite(paste("../out/Heatmap_", context, "union_DMRs_zscores.svg", sep = "_"), width=3.4, height=10.2)
print(ph.z)
#dev.off()

#####
# Summarise ntd counts 
meth_mean <- as.data.frame(meth.s.DMRs.merged.mat)
meth_mean$pos <- rownames(meth_mean)

# convert to long format
meth_mean <- pivot_longer(meth_mean, cols = 1:8, names_to = "sample_name", values_to = "meth")  %>% 
  # join with sample info table
  full_join(coldata, by = ("sample_name")) %>%
  # for each gene and condition
  group_by(pos, genotype) %>%
  # calculate the mean (scaled) meth
  summarise(mean_meth = mean(meth),
            nrep = n()) %>% 
  ungroup()

#Plot counts in clusters 
#Sort according to pheatmap rows
meth_cluster <- meth_cluster[ph.z$tree_row$order,]

meth_cluster <- right_join(meth_cluster, meth_mean, by = "pos")

head(meth_cluster)

plot_title <- "Cluster"

#Change order of clusters in plot to match heatmap
c.order <- unique(meth_cluster$cluster)
meth_cluster$cluster <- factor(meth_cluster$cluster, levels=c.order)

m.cl <- meth_cluster %>% 
  ggplot(aes(genotype, mean_meth, fill = genotype)) +
  geom_violin(width = 1, trim = TRUE) +
  scale_fill_manual(values = clrplt4) +
  scale_y_continuous(breaks = seq(0, 100, by = 50)) +
  theme_bw() +
  ggtitle(plot_title) +
  ylab(paste0(context,"m in DMRs")) +
  theme(axis.title.x=element_blank(), axis.title.y=element_text(size = 13), axis.text = element_text(size = 12), plot.title = element_text(size = 14, hjust = 0.5), legend.position="none") +
  geom_boxplot(width = 0.2, notch = FALSE, outlier.shape = 19, outlier.size = 0.02, outlier.alpha = 0.8, fill = "snow1") +
  geom_line(stat = "summary", fun = "median", colour = "brown", linewidth = 1, 
            aes(group = 1), alpha = 0.3) +
  facet_wrap(vars(cluster), ncol = 1)

#svglite(paste("../out/Cluster", context, "m_union_DMRs_2.svg", sep = "_"), width=1.8, height=5.4)
plot(m.cl)
#dev.off()

#Number of DMRs per cluster
dmr.c <- count(meth_cluster, cluster) %>% mutate(n = n/4)
#write_tsv(dmr.c,paste("../out/Number_of",context,"DMRs_per_cluster.tsv", sep="_"))

#Export cluster positions
clusters_pos_mean_meth <- pivot_wider(meth_cluster,names_from = genotype, values_from = mean_meth) %>%
  separate(pos, c("chr","start","end"), convert = TRUE)
#write_tsv(clusters_pos_mean_meth,paste0("../out/context,"_cluster_pos_mean_meth.tsv"))

#Clean up
rm(meth_cluster)
rm(meth_cluster.df)
rm(meth_mean)
rm(list=ls(pattern="myobj"))
rm(list=ls(pattern="meth"))
rm(list=ls(pattern="CG"))
rm(list=ls(pattern="cluster"))
rm(list=ls(pattern="pooled"))
gc()
