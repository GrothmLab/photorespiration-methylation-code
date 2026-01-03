### This scripts plots DNA methylation over gbM and CMT2-target regions and sample correlation for Arabidopsis plants 
### grown under long-day (LD) and short-day (SD) conditions (Hankover et al. (2026) Nature Plants).

#####
### Set working directory to directory containing the R script
#setwd("/path/to/photorespiration_methylation_code/R/WGBS_analysis_SD_LD/scripts")
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

library(methylKit)
library(genomation)
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)
library(readr)
library(RColorBrewer)
library(svglite)

########### WGBS_LD_SD_project ############

samplesheet <- read_delim("../samplesheet.csv", delim = ";")

colnames(samplesheet)[3] <- "condition"
samplesheet$condition <- gsub("_[1-2]","",samplesheet$condition)

(sample.id <- as.list(samplesheet$sample))

id_wt_ld <- sample.id[1:2]

id_wt_sd <- sample.id[3:4]

id_m_ld <- sample.id[5:6]

id_m_sd <- sample.id[7:8]


ids_treatment <- list(id_wt_sd,id_m_sd,id_m_ld)
ids_control <- list(id_wt_ld, id_m_ld)

idtreatments <- list("wt_sd","m_sd","m_ld")
idcontrols <- list("wt_ld","m_ld")

#Line colors for plots
clrplt4 <- c("indianred", "brown4", "steelblue", "dodgerblue4")

######

context <- "CG"

#read methylkit format CX reports

(file.list <- as.list(list.files(path = "../source-data/methylkit_in", pattern = "CG_report_sorted\\.txt$", full.names = TRUE)))

f_wt_ld <- file.list[1:2]
f_wt_sd <- file.list[3:4]
f_m_ld <- file.list[5:6]
f_m_sd <- file.list[7:8]

files_treatment <- list(f_wt_sd, f_m_sd, f_m_ld)
files_control <- list(f_wt_ld, f_m_ld)

#Generate methylkit object for all samples
myobj.s <- methRead(file.list,
                    sample.id,
                    treatment=c(rep(0,2), rep(1,2), rep(2,2), rep(3,2)),
                    pipeline="bismarkCytosineReport",
                    header=FALSE,
                    assembly="TAIR10",
                    context=context,
                    mincov=0,
                    dbtype=NA
)

#####
#Plot methylation over GbM (intersected with met1 DMRs), and CMT2 sites

#GbM sites
GbM_met1_sites <- readRDS("../source-data/GbM_met1_sites.RDS")

#CMT2 sites
cmt2_sites <- readBed("../../WGBS_analysis_m_t_mt_wt/source-data/cmt2_forddcc_hypoCHH_DMR_vs3reps_min4filter_mg200_sort.bed", zero.based = FALSE)

# Tile due to low coverage
myobj.s.100bp <- tileMethylCounts(myobj.s,win.size=100,step.size=100,cov.bases=0,mc.cores=2) %>%
  filterByCoverage(lo.count=NULL, lo.perc=NULL, hi.count=NULL, hi.perc=99.95)
meth.s.100bp <- methylKit::unite(myobj.s.100bp, destrand=FALSE, min.per.group=1L, mc.cores=2)  

pooled.meth.s.100bp <- pool(meth.s.100bp,sample.ids = c("WT_LD", "WT_SD", "m_LD", "m_SD"))

#####
# Plot over sites
#Count over specific sites
sites_list <- list(GbM_met1_sites, cmt2_sites)
names(sites_list) <- c("GbM-MET1 sites", "CMT2 sites")

boxvp_list <- list()

for (i in 1:length(sites_list)) {
  meth.s.sites <- regionCounts(pooled.meth.s.100bp, sites_list[[i]], cov.bases=10, mc.cores=2, save.db=FALSE)
  sites_meth <- as.data.frame(percMethylation(meth.s.sites)) %>%
    pivot_longer(cols = 1:4, names_to = "samples", values_to = "meth")
  
  plot_title <- paste0("m", context, " over ", names(sites_list)[i])
  boxvp_list[[i]] <- ggplot(sites_meth, aes(x=samples, y=meth, fill=as.factor(samples))) + 
    geom_violin(width = 1, trim = TRUE) +
    scale_fill_manual(values = clrplt4) +
    coord_cartesian(ylim=c(0,100), expand = TRUE, default = TRUE) +
    scale_y_continuous(breaks = seq(0, 100, by = 25)) +
    theme_bw() +
    ggtitle(plot_title) +
    ylab(paste0("m",context,"(%)")) +
    theme(axis.title.x=element_blank(), axis.title.y=element_text(size = 13), axis.text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(size = 14, hjust = 0.5), legend.position="none") +
    geom_boxplot(width = 0.3, notch = FALSE, outlier.shape = 19, outlier.size = 0.02, outlier.alpha = 0.8, fill = "snow1")
}

#svglite(paste0("../out/", context, "m_over_sites_m_wt_pooled_GbM_CMT2.svg"), width=4, height=2.5)
grid.arrange(boxvp_list[[1]],boxvp_list[[2]],
             nrow=1)
#dev.off()

#####
#Tile and unite C for correlation analysis
myobj.s.1kb <- tileMethylCounts(myobj.s,win.size=10000,step.size=10000,cov.bases=4,mc.cores=2) %>%
  filterByCoverage(lo.count=4, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)
meth.s.1kb <- methylKit::unite(myobj.s.1kb, destrand=FALSE, mc.cores=2)

#svglite(paste0("../out/", context, "m_1kb_dist_m_wt.svg"), width = 8, height = 6)
clusterSamples(meth.s.1kb, dist="euclidean", method="complete", plot=TRUE)
#dev.off()

######
context <- "CHG"

#read methylkit format CX reports

(file.list <- as.list(list.files(path = "../source-data/methylkit_in", pattern = "CHG_report_sorted\\.txt$", full.names = TRUE)))

f_wt_ld <- file.list[1:2]
f_wt_sd <- file.list[3:4]
f_m_ld <- file.list[5:6]
f_m_sd <- file.list[7:8]

files_treatment <- list(f_wt_sd, f_m_sd, f_m_ld)
files_control <- list(f_wt_ld, f_m_ld)

#Generate methylkit object for all samples
myobj.s <- methRead(file.list,
                    sample.id,
                    treatment=c(rep(0,2), rep(1,2), rep(2,2), rep(3,2)),
                    pipeline="bismarkCytosineReport",
                    header=FALSE,
                    assembly="TAIR10",
                    context=context,
                    mincov=0,
                    dbtype=NA
)

# Tile due to low coverage
myobj.s.100bp <- tileMethylCounts(myobj.s,win.size=100,step.size=100,cov.bases=0,mc.cores=2) %>%
  filterByCoverage(lo.count=NULL, lo.perc=NULL, hi.count=NULL, hi.perc=99.95)
meth.s.100bp <- methylKit::unite(myobj.s.100bp, destrand=FALSE, min.per.group=1L, mc.cores=2)  

pooled.meth.s.100bp <- pool(meth.s.100bp,sample.ids = c("WT_LD", "WT_SD", "m_LD", "m_SD"))

#####
# Plot over sites
#Count over specific sites

sites_list <- list(cmt2_sites)
names(sites_list) <- c("CMT2 sites")

boxvp_list <- list()

for (i in 1:length(sites_list)) {
  meth.s.sites <- regionCounts(pooled.meth.s.100bp, sites_list[[i]], cov.bases=10, mc.cores=2, save.db=FALSE)
  sites_meth <- as.data.frame(percMethylation(meth.s.sites)) %>%
    pivot_longer(cols = 1:4, names_to = "samples", values_to = "meth")
  
  plot_title <- paste0("m", context, " over ", names(sites_list)[i])
  boxvp_list[[i]] <- ggplot(sites_meth, aes(x=samples, y=meth, fill=as.factor(samples))) + 
    geom_violin(width = 1, trim = TRUE) +
    scale_fill_manual(values = clrplt4) +
    coord_cartesian(ylim=c(0,100), expand = TRUE, default = TRUE) +
    scale_y_continuous(breaks = seq(0, 100, by = 25)) +
    theme_bw() +
    ggtitle(plot_title) +
    ylab(paste0("m",context,"(%)")) +
    theme(axis.title.x=element_blank(), axis.title.y=element_text(size = 13), axis.text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(size = 14, hjust = 0.5), legend.position="none") +
    geom_boxplot(width = 0.3, notch = FALSE, outlier.shape = 19, outlier.size = 0.02, outlier.alpha = 0.8, fill = "snow1")
}

#svglite(paste0("../out/", context, "m_over_sites_m_wt_pooled_CMT2_2.svg"), width=2.3, height=3)
grid.arrange(boxvp_list[[1]],
             nrow=1)
#dev.off()

#####
#Tile and unite C for correlation analysis
myobj.s.1kb <- tileMethylCounts(myobj.s,win.size=10000,step.size=10000,cov.bases=4,mc.cores=2) %>%
  filterByCoverage(lo.count=4, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)
meth.s.1kb <- methylKit::unite(myobj.s.1kb, destrand=FALSE, mc.cores=2)

#svglite(paste0("../out/", context, "m_1kb_dist_m_wt.svg"), width = 8, height = 6)
clusterSamples(meth.s.1kb, dist="euclidean", method="complete", plot=TRUE)
#dev.off()

######
context <- "CHH"

#read methylkit format CX reports

(file.list <- as.list(list.files(path = "../source-data/methylkit_in", pattern = "CHH_report_sorted\\.txt$", full.names = TRUE)))

f_wt_ld <- file.list[1:2]
f_wt_sd <- file.list[3:4]
f_m_ld <- file.list[5:6]
f_m_sd <- file.list[7:8]

files_treatment <- list(f_wt_sd, f_m_sd, f_m_ld)
files_control <- list(f_wt_ld, f_m_ld)

#Generate methylkit object for all samples
myobj.s <- methRead(file.list,
                    sample.id,
                    treatment=c(rep(0,2), rep(1,2), rep(2,2), rep(3,2)),
                    pipeline="bismarkCytosineReport",
                    header=FALSE,
                    assembly="TAIR10",
                    context=context,
                    mincov=0,
                    dbtype=NA
)

#####
# Tile due to low coverage
myobj.s.100bp <- tileMethylCounts(myobj.s,win.size=100,step.size=100,cov.bases=0,mc.cores=2) %>%
  filterByCoverage(lo.count=4, lo.perc=NULL, hi.count=NULL, hi.perc=99.99)
meth.s.100bp <- methylKit::unite(myobj.s.100bp, destrand=FALSE, min.per.group=1L, mc.cores=2)

pooled.meth.s.100bp <- pool(meth.s.100bp,sample.ids = c("WT_LD", "WT_SD", "m_LD", "m_SD"))

#####
# Plot over sites
#Count over specific sites

sites_list <- list(cmt2_sites)
names(sites_list) <- c("CMT2 sites")

boxvp_list <- list()

for (i in 1:length(sites_list)) {
  meth.s.sites <- regionCounts(pooled.meth.s.100bp, sites_list[[i]], cov.bases=4, mc.cores=2, save.db=FALSE)
  sites_meth <- as.data.frame(percMethylation(meth.s.sites)) %>%
    pivot_longer(cols = 1:4, names_to = "samples", values_to = "meth")
  
  plot_title <- paste0(context, "m over ", names(sites_list)[i])
  boxvp_list[[i]] <- ggplot(sites_meth, aes(x=samples, y=meth, fill=as.factor(samples))) + 
    geom_violin(width = 1, trim = TRUE) +
    scale_fill_manual(values = clrplt4) +
    coord_cartesian(ylim=c(0,40), expand = TRUE, default = TRUE) +
    scale_y_continuous(breaks = seq(0, 40, by = 10)) +
    theme_bw() +
    ggtitle(plot_title) +
    ylab(paste0(context,"m (%)")) +
    theme(axis.title.x=element_blank(), axis.title.y=element_text(size = 13), axis.text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(size = 14, hjust = 0.5), legend.position="none") +
    geom_boxplot(width = 0.3, notch = FALSE, outlier.shape = 19, outlier.size = 0.02, outlier.alpha = 0.8, fill = "snow1")
  
}

#svglite(paste0("../out/", context, "m_over_sites_m_wt_pooled_CMT2sites.svg"), width=2.3, height=3)
grid.arrange(boxvp_list[[1]],
             nrow=1)
#dev.off()

#####
#Tile and unite C for correlation analysis
myobj.s.1kb <- tileMethylCounts(myobj.s,win.size=10000,step.size=10000,cov.bases=4,mc.cores=2) %>%
  filterByCoverage(lo.count=4, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)
meth.s.1kb <- methylKit::unite(myobj.s.1kb, destrand=FALSE, mc.cores=2)

#svglite(paste0("../out/", context, "m_1kb_dist_m_wt.svg"), width = 8, height = 6)
clusterSamples(meth.s.1kb, dist="euclidean", method="complete", plot=TRUE)
#dev.off()

#Clean up
rm(list=ls(pattern="myobj"))
rm(list=ls(pattern="meth"))
rm(list=ls(pattern="pooled"))
gc()

