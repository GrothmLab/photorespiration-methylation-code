### Correlation analysis between DETE expression and DNA methylation in m and wt grown under cCO2 and hCO2 conditions (Hankover et al. (2026) Nature PLants). ###

### Set working directory to directory containing the R script
#setwd("/path/to/photorespiration_methylation_code/R/WGBS_analysis_hCO2_cCO2/scripts")
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

#####
library(methylKit)
library(genomation)
library(rtracklayer)
library(tidyverse)
library(dplyr)
library(GenomicFeatures)
library(svglite)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(FSA)

#####
# Load gene expression data
#Load tpm data
#filename <- file.choose()
tpm_df <- read_tsv("../source-data/RNAseq2_salmon.merged.TE_tpm.tsv") 
keep <- rowSums(tpm_df[,3:14]) != 0
tpm_df <- tpm_df[keep,]

#Average and normalize gene expression data
tpm_df_longer <- dplyr::select(tpm_df, -2) %>%
  pivot_longer(cols = m_c_1:wt_h_3, names_to = "sample", values_to = "tpm", values_drop_na=FALSE) %>%
  mutate(block = rep(1:ceiling(n()/3), each = 3)[1:n()]) %>%
  group_by(block) %>%
  mutate(mean_tpm = mean(tpm),
         sample = str_replace(sample,"\\_[1-3]","")) %>% 
  ungroup() %>%
  dplyr::select(-block, -tpm) %>%
  slice(seq(1, n(), by = 3)) %>%
  mutate(mean_tpm = log2(mean_tpm + 1)) #log2 normalization to reduce variation

#####
#Create GRanges object for all differentially expressed TEs/genes
###Here make GRanges object for TEs
TE_gr <- import("../../WGBS_analysis_m_t_mt_wt/source-data/Panda_AtTEs_transcripts_v1_ensembl_fixed.gtf")
TE_txdb <- makeTxDbFromGRanges(TE_gr)
TE_gr <- genes(TE_txdb)

# Subset TEs to differentially expressed TEs
CO2_DETEs <- readRDS("../source-data/CO2_DETEs.rds")
TE_gr <- subset(TE_gr, names(TE_gr) %in% CO2_DETEs)

# Remove Pt and Mt
TE_gr <- subset(TE_gr, !seqnames(TE_gr) %in% c("Pt","Mt"))

#####

### Make GRanges objects with methyl data for correlation analysis
#####
options(scipen=999)

samplesheet <- read_csv("../samplesheet.csv")
(sample.id <- as.list(unique(sort(samplesheet$name))))
(coldata <- read_csv("../WGBS_CO2_sample_info.csv"))

id1 <- as.list(sample.id[1:2])
id2 <- as.list(sample.id[3:4])
id3 <- as.list(sample.id[5:6])
id4 <- as.list(sample.id[7:8])

#ids <- list(id1,id2,id3)
sample.id.t <- list(id1,id2,id4)
sample.id.c <- list(id3)

idtreatments <- list("m_c","m_h","wt_h")
idcontrols <- list("wt_c", "m_h")

#Line colors for metaplots
clrplt4 <- c("indianred1", "indianred4", "steelblue1", "steelblue4")
clrplt8 <- c("indianred1", "indianred1", "indianred4", "indianred4", "steelblue1", "steelblue1","steelblue4", "steelblue4")

#Line colors for heatmap
col.ramp <- colorRampPalette(brewer.pal(9, "Reds") )(255)

#####
context <- "CG"

#read methylkit format CX reports

(file.list <- as.list(list.files(path = "../source-data/methylkit_in", pattern = "CG_report_sorted\\.txt$", full.names = TRUE)))
f1 <- file.list[1:2]
f2 <- file.list[3:4]
f3 <- file.list[5:6]
f4 <- file.list[7:8]

treatments <- list(f1,f2,f4)
controls <- list(f3)

#Generate methylkit object for all samples
myobj.s <- methRead(file.list[c(1:4,7:8,5:6)],
                    sample.id[c(1:4,7:8,5:6)],
                    treatment=c(1,1,2,2,3,3,0,0),
                    pipeline="bismarkCytosineReport",
                    header=FALSE,
                    assembly="TAIR10",
                    context=context,
                    mincov=4,
                    dbtype=NA
)

#Pool replicates per genotype
myobj.s.m.c <- reorganize(myobj.s,
                        c("mR_c_1","mR_c_2"),
                        treatment=c(0,0))
meth.s.m.c <- methylKit::unite(myobj.s.m.c, destrand=FALSE, min.per.group=1L, mc.cores=2)
pooled.meth.s.m.c <- pool(meth.s.m.c,sample.ids = "m_c")

myobj.s.m.h <- reorganize(myobj.s,
                         c("mR_h_1","mR_h_2"),
                         treatment=c(0,0))
meth.s.m.h <- methylKit::unite(myobj.s.m.h, destrand=FALSE, min.per.group=1L, mc.cores=2)
pooled.meth.s.m.h <- pool(meth.s.m.h,sample.ids = "m_h")

myobj.s.wt.c <- reorganize(myobj.s,
                        c("wtR_c_1","wtR_c_2"),
                        treatment=c(0,0))
meth.s.wt.c <- methylKit::unite(myobj.s.wt.c, destrand=FALSE, min.per.group=1L, mc.cores=2)
pooled.meth.s.wt.c <- pool(meth.s.wt.c,sample.ids = "wt_c")

myobj.s.wt.h <- reorganize(myobj.s,
                         c("wtR_h_1","wtR_h_2"),
                         treatment=c(0,0))
meth.s.wt.h <- methylKit::unite(myobj.s.wt.h, destrand=FALSE, min.per.group=1L, mc.cores=2)
pooled.meth.s.wt.h <- pool(meth.s.wt.h,sample.ids = "wt_h")

pooled.meth.s <- c()
pooled.meth.s <- list(pooled.meth.s.m.c, pooled.meth.s.m.h, pooled.meth.s.wt.c, pooled.meth.s.wt.h)

#Convert to GRanges objects for metaplotting
myobj.grl <- c()
for (i1 in 1:length(pooled.meth.s)) {
  myobj.grl <- c(myobj.grl, as(pooled.meth.s[[i1]],"GRanges"))
  colnames(mcols(myobj.grl[[i1]])) <- c("coverage","numCs","numTs")
}
names(myobj.grl) <- c("m_c","m_h","wt_c","wt_h")

#####
### Now using myobj.grl and TE_gr for correlation analysis

#Set variables
annotation <- TE_gr
context <- context
data <- myobj.grl
windowsize=250; insidewindows=4; range=1000

#Subset data to annotated regions +/- range
#Firts convert data from list to GRangesList
data <- GRangesList(data)

subsetByOverlaps_list <- function(query_gr, subject_gr) {
  subsetByOverlaps(query_gr, subject_gr, ignore.strand=TRUE)
}

#Remove annotations that do not overlap with data
for (i in 1:length(data)){
  annotation <- subsetByOverlaps(annotation, data[[i]], ignore.strand=TRUE, minoverlap = 1L) #minoverlap = length(data) makes sure that all sapmles overlap with data
}
#Resize annotation
data_range <- resize(x = annotation, width = width(annotation) + range, fix = "start")
data_range <- resize(x = data_range, width = width(data_range) + range, fix = "end")
data_range <- sort(data_range)

#Subset data
data <- lapply(data, subsetByOverlaps_list, data_range)
#for (i in 1:length(data)){
#  data[[i]] <- subsetByOverlaps(data[[i]], data_range, ignore.strand=TRUE)
#}
message("Generating metaplots with")
message("data: ")
print(data)
message("context: ", context)
message("up-/downstream window size (bp): ", windowsize)
message("up-/downstream range (bp): ", range)
message("inside bin number: ", insidewindows)

#breaks <- c(c(-1, -0.5, 0, (1-1/insidewindows), (1-1/insidewindows+0.5), (1-1/insidewindows+1)) * range)
breaks <- c(c(-1, 0, (1-1/insidewindows), (1-1/insidewindows+1)) * range)
#labels <- c(paste0(-range/1000,"kb"), paste0(-range/2000,"kb"), 'Start', 'Stop', paste0(range/2000,"kb"), paste0(range/1000,"kb"))
labels <- c(paste0(-range/1000,"kb"), 'start', 'stop', paste0(range/1000,"kb"))

category <- context

overlaps.meth.lvl <- array(NA, dim=c(length(data), 1, range%/%windowsize, length(annotation), 2), dimnames=list(sample=names(data), category=context, distance=1:(range%/%windowsize),  gene=names(annotation), where=c('upstream', 'downstream')))
ioverlaps.meth.lvl <- array(NA, dim=c(length(data), 1, insidewindows, length(annotation), 1), dimnames=list(sample=names(data), category=context, distance=1:insidewindows, gene=names(annotation), where=c('inside')))

for (i in 1:length(data)){
  
  if (class(data[[i]]) == 'GRanges') {
    message("Calculating upstream, downstream and inside counts for ", names(data)[[i]])
  } else {
    stop("Argument 'data' needs to be of class 'GRanges'.")
  }
  
  
  #Loop through each gene
  for (j in 1:(length(annotation))) {
    ## Subset annotation to data range
    annotation.s <- subsetByOverlaps(annotation[j], data[[i]], ignore.strand=TRUE)
    
    ## Subset data to annotation +/- range
    data_range <- resize(x = annotation.s, width = width(annotation.s) + range, fix = "start")
    data_range <- resize(x = data_range, width = width(data_range) + range, fix = "end")
    data.s <- subsetByOverlaps(data[[i]], data_range)
    
    ## Upstream and downstream annotation
    annotation.up <- resize(x = annotation.s, width = 1, fix = 'start')
    annotation.down <- resize(x = annotation.s, width = 1, fix = 'end')
    
    for (i1 in 1:(range %/% windowsize)) {
      anno.up <- suppressWarnings( promoters(x = annotation.up, upstream = i1*windowsize, downstream = 0) )
      anno.up <- suppressWarnings( resize(x = anno.up, width = windowsize, fix = 'start') )
      ind <- findOverlaps(anno.up, data.s, ignore.strand=TRUE)
      overlaps.meth.lvl[names(data)[[i]], category, as.character(i1), names(annotation.s), 'upstream'] <- (sum(data.s$numCs[ind@to], na.rm=TRUE)) / (sum(data.s$coverage[ind@to], na.rm=TRUE))
    }
    for (i1 in 1:(range %/% windowsize)) {
      anno.down <- suppressWarnings( promoters(x = annotation.down, upstream = 0, downstream = i1*windowsize) )
      anno.down <- suppressWarnings( resize(x = anno.down, width = windowsize, fix = 'end') )
      ind <- findOverlaps(anno.down, data.s, ignore.strand=TRUE)
      overlaps.meth.lvl[names(data)[[i]], category, as.character(i1), names(annotation.s), 'downstream'] <- (sum(data.s$numCs[ind@to], na.rm=TRUE)) / (sum(data.s$coverage[ind@to], na.rm=TRUE))
    }
    
    ## Inside annotation
    widths <- width(annotation.s)
    
    for (i1 in 1:insidewindows) {
      anno.in <- annotation.s
      if (as.logical(strand(annotation.s) == '+')) {
        end(anno.in) <- start(annotation.s) + round(widths * i1 / insidewindows)
        start(anno.in) <- start(annotation.s) + round(widths * (i1-1) / insidewindows)
      } else {  
        start(anno.in) <- end(annotation.s) - round(widths * i1 / insidewindows)
        end(anno.in) <- end(annotation.s) - round(widths * (i1-1) / insidewindows)
      }
      ind <- findOverlaps(anno.in, data.s, ignore.strand=TRUE)
      
      ioverlaps.meth.lvl[names(data)[[i]], category, as.character(i1), names(annotation.s), 'inside'] <- (sum(data.s$numCs[ind@to], na.rm=TRUE)) / (sum(data.s$coverage[ind@to], na.rm=TRUE))
    }
  }      
}

## Prepare data.frames for plotting
dfo <- reshape2::melt(overlaps.meth.lvl)
dfo$distance <- as.numeric(dfo$distance)

# Adjust distances for plot
dfo$distance[dfo$where == 'upstream'] <- -dfo$distance[dfo$where == 'upstream'] * windowsize + windowsize/2
dfo$distance[dfo$where == 'downstream'] <- dfo$distance[dfo$where == 'downstream'] * windowsize + range - range/insidewindows - windowsize/2

dfi <- reshape2::melt(ioverlaps.meth.lvl)
dfi$distance <- as.numeric(dfi$distance)

# Adjust distances for plot
dfi$distance <- (dfi$distance -1) * range / insidewindows

dfm <- rbind(dfo, dfi)
colnames(dfm)[4] <- "gene_id"

#####
# Calculate correlation between DNA methylation and expression for each window

df_all <- merge(dfm, tpm_df_longer, by = c("gene_id", "sample"), all.x = TRUE, sort = FALSE)

#Normalize methylation
df_all <- df_all %>%
  group_by(gene_id) %>%
  arrange(distance, sample, .by_group = TRUE) %>%
  dplyr::select(-category, -where) %>%
  ungroup() %>%
  mutate(nt_value = log2(value + 0.01))

#####
#Show value densities
ggplot(df_all, aes(x=value)) + geom_density(aes(group=sample,color=sample,fill=sample), alpha=0.3, bw="nrd")

p <- ggplot(df_all, aes(x=nt_value)) + geom_density(aes(group=sample,color=sample,fill=sample), alpha=0.3, bw="nrd")
p + geom_vline(xintercept = -4, color = "red", linetype = "dashed", size = 1) #+ coord_cartesian(xlim = c(0, 0.2),) # ylim = c(-1, 1)

df_wt <- subset(df_all, df_all$sample %in% c("wt_c","wt_h"))

ggplot(df_wt, aes(x=value)) + geom_density(aes(group=sample,color=sample,fill=sample), alpha=0.3, bw="nrd")

p <- ggplot(df_wt, aes(x=nt_value)) + geom_density(aes(group=sample,color=sample,fill=sample), alpha=0.3, bw="nrd")
p + geom_vline(xintercept = -4, color = "red", linetype = "dashed", size = 1) #+ coord_cartesian(xlim = c(0, 0.2),) # ylim = c(-1, 1)

summary(df_wt$nt_value)

#Subset genes to WT-methylated TEs (mCG ration > 0.25) as these are usually silenced
df_wt_m <- subset(df_wt, df_wt$value >= 0.25)
ggplot(df_wt_m, aes(x=nt_value)) + geom_density(aes(group=sample,color=sample,fill=sample), alpha=0.3, bw="nrd") 
meth_genes <- unique(df_wt_m$gene_id)
length(meth_genes)

#####
# Remove genes below methylation threshold
df_all_m <- subset(df_all, df_all$gene_id %in% meth_genes)

#Reshape the data frame for plotting or export
df_all_m_list <- split(df_all_m, df_all_m$sample)

# Define a function to pivot wider for each data frame
pivot_wider_df <- function(df) {
  pivoted_df <- pivot_wider(df[,-6], names_from = distance, values_from = value)
  return(pivoted_df)
}

# Apply the function to each data frame in df_list
df_all_m_list_w <- lapply(df_all_m_list, pivot_wider_df)

# Only check m correlation
df_all_m_list_w_m <- list()
df_all_m_list_w_m <- df_all_m_list_w[1:2]

# Define a function to make matrix for each data frame
make_meth_matrix <- function(df) {
  meth_matrix <- as.matrix(df[,-(1:3)])
  rownames(meth_matrix) <- df$gene_id
  return(meth_matrix)
}

make_meth_matrix_inside <- function(df) {
  as.matrix(df[8:11])
  #rownames(meth_matrix) <- df$gene_id
  #return(meth_matrix)
}

# Apply the function to each data frame in df_list
meth_mat_list <- lapply(df_all_m_list_w, make_meth_matrix)

meth_mat_list_inside <- lapply(df_all_m_list_w, make_meth_matrix_inside)

# Get mean meth levels for each row
meth_mat_list_mean <- lapply(meth_mat_list_inside, rowMeans, na.rm = TRUE)

meth_mat_mean <- do.call(cbind, meth_mat_list_mean)
colnames(meth_mat_mean) <- c("m_c_meth", "m_h_meth", "wt_c_meth", "wt_h_meth") #names(meth_mat_list_mean)
rownames(meth_mat_mean) <- df_all_m_list_w[[1]]$gene_id

#Get count data from df_all, as this has the right ordering and reshape for plotting
mean_tpm_tbl <- df_all_m[which(df_all_m$distance == 0),c(1,2,5)] %>% 
  pivot_wider(names_from = c(sample), values_from = mean_tpm) 

#Calculate Spearman correlations for m and mt only
cr_m_sp <- df_all_m %>%
  filter(sample %in% c("m_c","m_h")) %>%
  group_by(gene_id, distance) %>%
  summarise(spearman_correlation = cor(as.numeric(value), as.numeric(mean_tpm), method = "spearman")) %>%
  ungroup()

cr_m_sp_w <- as.data.frame(pivot_wider(cr_m_sp, names_from = c(distance), values_from = spearman_correlation))
rownames(cr_m_sp_w) <- cr_m_sp_w$gene_id

#save(meth_mat_list, c_mat, cr_m_sp_w, df_all_m, file="../out/DETEs_CG_corr_data_new.RData")
#saveRDS(meth_genes,"../out/CGmeth_DETEs.rds")

#####
#Combine data into one table
master_table_mean <- as.data.frame(mean_tpm_tbl)
rownames(master_table_mean) <- mean_tpm_tbl$gene_id

# Annotate rows in the heatmap (here I want to see the DETEs for each comparison)
# Import gene names for each category
mc_up <- readRDS("../source-data/DETEsUp_mc_wt.rds") 
mh_up <- readRDS("../source-data/DETEsUp_mh_wt.rds")
mh_mc_down <- readRDS("../source-data/DETEsDown_mh_mc.rds")

DETEs_list <- list()
DETEs_list <- list(mc_up, mh_up, mh_mc_down)
names(DETEs_list) <- c("mc_wtc_up","mh_wtc_up","mh_mc_down")

DETE_names_list <- list(rownames(mc_up), rownames(mh_up), rownames(mh_mc_down))
names(DETE_names_list) <- c("mc_wtc_up","mh_wtc_up","mh_mc_down")

master_table_mean <- master_table_mean[,-1]
for (m in names(DETE_names_list)){
  anno_rows_2 <- as.data.frame(DETE_names_list[[m]])
  rownames(anno_rows_2) <- anno_rows_2$`DETE_names_list[[m]]`
  label <- m
  anno_rows_2$label <- label
  anno_rows_2$`DETE_names_list[[m]]` <- NULL
  master_table_mean <- merge(master_table_mean, anno_rows_2, by="row.names", all.x=TRUE)
  rownames(master_table_mean) <- master_table_mean$Row.names
  master_table_mean <- master_table_mean[,-1]
  master_table_mean[is.na(master_table_mean)] <- "no" 
}
#master_table_mean$label[is.na(master_table_mean$label)] <- "other"

colnames(master_table_mean)[5:7] <- names(DETE_names_list)

master_table_mean$mc_wtc_up <- as.factor(master_table_mean$mc_wtc_up)
master_table_mean$mh_wtc_up <- as.factor(master_table_mean$mh_wtc_up)
master_table_mean$mh_mc_down <- as.factor(master_table_mean$mh_mc_down)

master_table_all <- merge(master_table_mean, meth_mat_mean, by="row.names")
rownames(master_table_all) <- master_table_all$Row.names
master_table_all <- master_table_all[,-1]
master_table_all <- merge(master_table_all, cr_m_sp_w, by="row.names")
rownames(master_table_all) <- master_table_all$Row.names
master_table_all <- master_table_all[,-1]
master_table_all <- master_table_all[,c(12,5:7,1:4,8:11,13:24)]
colnames(master_table_all)[c(1,5:12)] <- c("GeneID", "m_c_ntpm","m_h_ntpm","wt_c_ntpm","wt_h_ntpm","m_c_meth","m_h_meth","wt_c_meth","wt_h_meth")

#####
# Adding log2FC values to master_table_all

#load log2FC data
filename <- "../source-data/res_lfcShrink_ddsC_PR.RData"

load(filename)

#See what was loaded
#Create a new environment
new_env <- new.env()
load(filename, envir = new_env)
variables <- ls(envir = new_env)
print(variables)

#Change into data frames with log2LFC data
res_m_c_wt_c_lfc <- as.data.frame(res_m_c_wt_c[,1:2])
colnames(res_m_c_wt_c_lfc) <- c("baseMean_mc","log2FC_mc")

res_m_h_wt_c_lfc <- as.data.frame(res_m_h_wt_c[,1:2])
colnames(res_m_h_wt_c_lfc) <- c("baseMean_mh","log2FC_mh")

res_m_h_m_c_lfc <- as.data.frame(res_m_h_m_c[,1:2])
colnames(res_m_h_m_c_lfc) <- c("baseMean_mhVSmc","log2FC_mhVSmc")

res_wt_h_wt_c_lfc <- as.data.frame(res_wt_h_wt_c[,1:2])
colnames(res_wt_h_wt_c_lfc) <- c("baseMean_wth","log2FC_wth")

#Merge lfc data with hclust matrix
master_table_all <- merge(master_table_all, res_m_c_wt_c_lfc, by="row.names", all.x = TRUE)
rownames(master_table_all) <- master_table_all$Row.names
master_table_all <- master_table_all[,-1]

master_table_all <- merge(master_table_all, res_m_h_wt_c_lfc, by="row.names", all.x = TRUE)
rownames(master_table_all) <- master_table_all$Row.names
master_table_all <- master_table_all[,-1]

master_table_all <- merge(master_table_all, res_m_h_m_c_lfc, by="row.names", all.x = TRUE)
rownames(master_table_all) <- master_table_all$Row.names
master_table_all <- master_table_all[,-1]

master_table_all <- merge(master_table_all, res_wt_h_wt_c_lfc, by="row.names", all.x = TRUE)
rownames(master_table_all) <- master_table_all$Row.names
master_table_all <- master_table_all[,-1]

# Add family info
filename <- "../source-data/salmon_tx2gene.tsv"
tx2gene = filename
info = file.info(tx2gene)
if (info$size == 0) {
  tx2gene = NULL
} else {
  rowdata = read.csv(tx2gene, sep="\t", header = FALSE)
  colnames(rowdata) = c("tx", "gene_id", "gene_name")
  tx2gene = rowdata[,1:2]
}

TEanno <- read_tsv("../source-data/TAIR_IDs_Panda_AtTEs_annotation_v1_ensembl.tsv") %>% rename(tx = Transcript_ID)
TEanno <- left_join(TEanno, tx2gene, by = "tx")
#Problem: Some transcript ids match to multiple TAIR IDs and vice versa
#Not a problem: Count data is based on gene IDs, TE annotation file contains transcript IDs, but if there are multiple transcripts for a TE annotation they all belong to the same gene ID. So for annotation we can simply ignore duplicated entries per gene ID.
TEanno_uniq <- TEanno %>% distinct(gene_id, .keep_all = TRUE) %>% dplyr::select(SubFamily, Family, gene_id) %>% rename(GeneID = gene_id)

master_table_all <- left_join(master_table_all, TEanno_uniq, by="GeneID")
rownames(master_table_all) <- master_table_all$GeneID
#master_table_all <- master_table_all[,-1]

# Add additional columns to subset data into different intersections
master_table_all <- master_table_all %>%
  mutate(section = paste(mc_wtc_up, mh_wtc_up, mh_mc_down, sep = "_")) %>%
  mutate(context = "CG")

#####
### Plot log2FC, DNA methylation, and correlation data for different intersects of PR DETEs
# Reduce table to only two sections: "mc_wtc_up_mh_wtc_up_no" (CO2is) and the others combined (CO2dy)
master_table_all_2sec <- master_table_all
master_table_all_2sec$section <- gsub("mc_wtc_up_mh_wtc_up_no", "C02is", master_table_all_2sec$section)
master_table_all_2sec$section <- gsub("mc_wtc_up_no_no|mc_wtc_up_mh_wtc_up_mh_mc_down|mc_wtc_up_no_mh_mc_down", "C02dy", master_table_all_2sec$section)

# Remove sections that are irrelevant and plot LFCs
lfc_section <- master_table_all_2sec %>%
  filter(Family != ".") %>%
  dplyr::select(c(26,28,34)) %>%
  pivot_longer(cols = log2FC_mc:log2FC_mh, names_to = "comparison", values_to = "log2FC")

cl_counts <- lfc_section %>% 
  ggplot(aes(Family, log2FC, fill = Family)) +
  geom_boxplot(aes(color = comparison), position = position_dodge(width = 0.8), width = 0.7, notch = FALSE,  outlier.shape = NA, outlier.size = 0.5, alpha = 0.9, outlier.alpha =0.9, linewidth = 0.4) + #outlier.shape = NA,
  coord_cartesian(ylim = c(-3,16)) +
  ylab(paste0("log2 fold change")) +
  scale_color_manual(values = c("indianred", "indianred4")) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() +
  theme(text = element_text(size = 12), axis.title.x=element_blank(), axis.text.x = element_blank(), plot.title = element_text(size = 14, hjust = 0.5), legend.position = "right")
#geom_line(stat = "summary", fun = "median", colour = "brown", linewidth = 1, 
#          aes(group = 1), alpha = 0.3) +
#facet_wrap(vars(section), ncol = 1, scales = "fixed")

#svg("../out/DETEs_log2fc_mc_mh_Famililes.svg", width=4, height=3) #ED Fig. 7e
plot(cl_counts)
#dev.off()

#for (i in c("Copia", "En-Spm", "Gypsy", "Harbinger", "HAT", "Helitron", "L1", "MuDR")) {
#  family_count <- length(lfc_section$Family[lfc_section$Family == i])
#  print(family_count/2)
#}

# Kruskal-Wallis test for each class
kruskal_results <- lfc_section %>%
  group_by(Family) %>%
  summarize(
    kruskal_p_value = kruskal.test(log2FC ~ comparison)$p.value
  )

# Run Dunn's test for each class if Kruskal-Wallis is significant
dunn_results <- lfc_section %>%
  group_by(Family) %>%
  do(
    dunnTest(log2FC ~ comparison, data = ., method = "bonferroni")$res
  ) %>%
  ungroup()

# View results
print(kruskal_results)
print(dunn_results, n = Inf)
#write_tsv(kruskal_results, paste0("../out/Kruskal_DETEs_log2fc_mc_mh_Famililes.tsv"))
#write_tsv(dunn_results, paste0("../out/Dunn_DETEs_log2fc_mc_mh_Famililes.tsv"))

# Remove sections that are irrelevant and plot LFCs
lfc_section <- master_table_all_2sec %>%
  dplyr::select(c(26,28,30,32,35)) %>%
  filter(section == "C02is" | section == "C02dy") %>%
  pivot_longer(cols = log2FC_mc:log2FC_wth, names_to = "comparison", values_to = "log2FC")

cl_counts <- lfc_section %>%
  ggplot(aes(comparison, log2FC, fill = comparison)) +
  geom_line(stat = "summary", fun = "median", colour = "brown", linewidth = 1, 
            aes(group = 1), alpha = 0.3) +
  geom_boxplot(aes(color = section), position = position_dodge(width = 0.8), width = 0.7, notch = TRUE,  outlier.shape = NA,linewidth = 0.4) + #outlier.size = 0.5, alpha = 0.9, outlier.alpha =0.9,
  scale_color_manual(values = c("gray20", "gray60")) +
  scale_fill_manual(values = c("indianred", "indianred4", "papayawhip", "steelblue4")) +
  coord_cartesian(ylim = c(-5,15)) +
  ylab(paste0("log2 fold change")) +
  theme_bw() +
  theme(text = element_text(size = 12), axis.title.x=element_blank(), axis.text.x = element_blank(), plot.title = element_text(size = 14, hjust = 0.5), legend.position = "right")

#svg("../out/DETEs_log2fc_2sections_dodged.svg", width=4, height=1.8) #ED Fig. 7d
plot(cl_counts)
#dev.off()

#for (i in c("C02dy", "C02is")) {
#  section_count <- length(lfc_section$section[lfc_section$section == i])
#  print(section_count/4)
#}

# Kruskal-Wallis test for each class
kruskal_results <- lfc_section %>%
  group_by(comparison) %>%
  summarize(
    kruskal_p_value = kruskal.test(log2FC ~ section)$p.value
  )

# Run Dunn's test for each class if Kruskal-Wallis is significant
dunn_results <- lfc_section %>%
  group_by(comparison) %>%
  do(
    dunnTest(log2FC ~ section, data = ., method = "bonferroni")$res
  ) %>%
  ungroup()

# View results
print(kruskal_results)
print(dunn_results, n = Inf)
#write_tsv(kruskal_results, paste0("../out/Kruskal_DETEs_log2fc_2sections.tsv"))
#write_tsv(dunn_results, paste0("../out/Dunn_DETEs_log2fc_2sections.tsv"))

# Remove sections that are irrelevant and plot correlation
scc_section <- master_table_all_2sec %>%
  dplyr::select(c(35,13:24)) %>%
  filter(section == "C02is" | section == "C02dy") %>%
  pivot_longer(cols = 2:13, names_to = "position", values_to = "pcc")

cor_summary <- scc_section %>%
  group_by(section, position) %>% #category,
  summarize(mean_pcc = mean(pcc, na.rm = TRUE), .groups = "drop") %>%
  mutate(position = as.integer(position)) %>%
  arrange(section,position)

cl_counts <- cor_summary %>% 
  ggplot(aes(position, mean_pcc, color = section, group = section)) + #interaction(Family, class)
  geom_smooth(linewidth = 1.3, method = "loess", se = FALSE, alpha = 0.3) + #interaction(class, category)
  coord_cartesian(ylim = c(-1,1)) +
  ggtitle("Spearman correlation") +
  scale_color_manual(values = brewer.pal(n=9, name = "Set3")[1:4]) +
  theme_bw() +
  theme(text = element_text(size = 12), axis.title=element_blank(), plot.title = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

#svg("../out/DETEs_mCG_scc_monly_2sections.svg", width=2.5, height=1.5) #ED Fig. 7d
plot(cl_counts)
#dev.off()

# Plot TE families in each subgroup
fam_section <- master_table_all_2sec %>%
  dplyr::select(c(34,35)) %>%
  filter(Family != ".") %>%
  filter(section == "C02is" | section == "C02dy")

fam_ratios <- fam_section %>%
  group_by(section, Family) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(section) %>%
  mutate(Ratio = Count / sum(Count)) %>%
  unnest(c(Ratio))

cl_counts <- fam_ratios %>%
  ggplot(aes(x= 1, y = Ratio, fill = Family)) +  
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = brewer.pal(length(unique(fam_ratios$Family)),("Set3"))) +
  labs(y = "Ratio") +
  theme_minimal() +                                                    
  theme(text = element_text(size = 12), axis.title.x=element_blank(), plot.title = element_text(size = 4), axis.text.x = element_blank(), legend.position = "right") +
  facet_wrap(vars(section), ncol = 4)

#svglite::svglite("../out/Barchart_DETEs_families_2sections.svg", width=2.6, height=2.4)
plot(cl_counts)
#dev.off()

# Plot barchart of familiy sizes for all mthfd1 DETEs
fam_ratios_2 <- fam_section %>%
  group_by(Family) %>%
  summarise(Count = n(), .groups = "drop") %>%
  #group_by(section) %>%
  mutate(Ratio = Count / sum(Count)) %>%
  unnest(c(Ratio))

cl_counts_2 <- fam_ratios_2 %>%
  ggplot(aes(x= 1, y = Ratio, fill = Family)) +  
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = brewer.pal(length(unique(fam_ratios_2$Family)),("Set3"))) +
  labs(y = "Ratio") +
  theme_minimal() +                                                    
  theme(text = element_text(size = 12), axis.title.x=element_blank(), plot.title = element_text(size = 4), axis.text.x = element_blank(), legend.position = "right")
#facet_wrap(vars(section), ncol = 4)

#svglite::svglite("../out/Barchart_DETEs_families_all_mc.svg", width=2.2, height=2.4)
plot(cl_counts_2)
#dev.off()

# Plot family sizes in all expressable TEs
# Read full annotation
TE_fullanno <- read_tsv("../source-data/Panda_AtTEs_annotation_v1_ensembl.tsv")

fam_all_expr <- TE_fullanno %>%
  dplyr::select(c(4,8,19)) %>%
  filter(Family != ".") %>%
  filter(Expression == "ExpressedAndAnnotated")

colnames(fam_all_expr)[c(1)] <- c("TAIR_ID")
#fam_all_expr <- merge(TEanno_uniq[,-5], fam_all_expr, by="TAIR_ID")

fam_ratios_3 <- fam_all_expr %>%
  group_by(Family) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Ratio = Count / sum(Count)) %>%
  unnest(c(Ratio))

cl_counts_3 <- fam_ratios_3 %>%
  ggplot(aes(x= 1, y = Ratio, fill = Family)) +  
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = brewer.pal(length(unique(fam_ratios_3$Family)),("Set3"))) +
  labs(y = "Ratio") +
  theme_minimal() +                                                    
  theme(text = element_text(size = 12), axis.title.x=element_blank(), plot.title = element_text(size = 4), axis.text.x = element_blank(), legend.position = "right")

#svglite::svglite("../out/Barchart_DETEs_families_all_exprANDanno.svg", width=2.2, height=2.4)
plot(cl_counts_3)
#dev.off()

# Test for significant differences in ratios
# Use lapply to add the dataset name to each dataset, then combine with bind_rows

# Add data set column
fam_ratios_3 <- fam_ratios_3 %>%
  mutate(section = "all")

fam_ratios_2 <- fam_ratios_2 %>%
  mutate(section = "mc_DETEs")

fam_ratios <- fam_ratios[,c(2:4,1)]

fam_ratios_list <- list(fam_ratios,fam_ratios_2,fam_ratios_3)
names(fam_ratios_list) <- c("secDETES","mcDETES","allDETEs")

combined_data <- bind_rows(
  lapply(names(fam_ratios_list), function(name) {
    fam_ratios_list[[name]] %>% mutate(Dataset = name)
  })
)

# Create a contingency table
contingency_table <- xtabs(Count ~ Family + section, data = combined_data)

# Run Chi-squared test for homogeneity
chi_test_result <- chisq.test(contingency_table)
print(chi_test_result)

# Example of pairwise comparison for each Category
pairwise_tests <- combined_data %>%
  group_by(section) %>%
  mutate(Totals = sum(Count)) %>%
  group_by(Family) %>%
  summarize(
    Counts = list(Count),               # Store the counts for each group
    Sums = list(Totals),  # Store n values for each group
    test_result = list(
      pairwise.prop.test(
        x = Count,
        n = Totals,  # Set n as the total counts for each dataset within the Value group
        p.adjust.method = "bonferroni"
      )
    ),
    .groups = "drop"
  )

# Inspecting each `Value` with counts and totals
pairwise_tests %>%
  mutate(
    Counts = map_chr(Counts, ~ paste(.x, collapse = ", ")),
    Sums = map_chr(Sums, ~ paste(.x, collapse = ", "))
  ) %>%
  dplyr::select(Family, Counts, Sums, test_result)

# Print results
names(pairwise_tests$test_result) <- pairwise_tests$Family

#sink("../out/Pairwise_test_fam_ratios_new.txt")
print(chi_test_result)
# Inspecting each `Value` with counts and totals
pairwise_tests %>%
  mutate(
    Counts = map_chr(Counts, ~ paste(.x, collapse = ", ")),
    Sums = map_chr(Sums, ~ paste(.x, collapse = ", "))
  ) %>%
  dplyr::select(Family, Counts, Sums, test_result)
pairwise_tests$test_result
#sink()

#Plot all sections in one plot
#Set order of sections for plotting
combined_data$section <- factor(combined_data$section, levels = c("all", "mc_DETEs", "C02is", "C02dy"))

cl_counts <- combined_data %>%
  ggplot(aes(x= 1, y = Ratio, fill = Family)) +  
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = brewer.pal(length(unique(fam_ratios_3$Family)),("Set3"))) +
  labs(y = "Ratio") +
  theme_minimal() +                                                    
  theme(text = element_text(size = 12), axis.title.x=element_blank(), plot.title = element_text(size = 4), axis.text.x = element_blank(), legend.position = "right") +
  facet_wrap(vars(section), ncol = 4)

#svglite::svglite("../out/Barchart_DETEs_families_new_final_2.svg", width=3.2, height=2.4) #ED Fig. 7e
plot(cl_counts)
#dev.off()

meth_section <- master_table_all_2sec %>%
  dplyr::select(c(9:12,35)) %>%
  filter(section == "C02is" | section == "C02dy") %>%
  pivot_longer(cols = m_c_meth:wt_h_meth, names_to = "comparison", values_to = "meth")

cl_counts <- meth_section %>%
  ggplot(aes(comparison, meth, fill = comparison)) +
  geom_line(stat = "summary", fun = "median", colour = "brown", linewidth = 1, 
            aes(group = 1), alpha = 0.3) +
  scale_fill_manual(values = c("indianred", "indianred4", "papayawhip", "steelblue4")) +
  ylab(paste0("mCG")) +
  theme_bw() +
  geom_boxplot(aes(color = section), position = position_dodge(width = 0.8), width = 0.6, notch = FALSE, outlier.shape = NA, outlier.size = 0.02, outlier.alpha = 0.8 ) + #fill = "snow1"
  scale_color_manual(values = c("gray20", "gray60")) +
  theme(text = element_text(size = 12), axis.title.x=element_blank(), axis.text.x = element_blank(), plot.title = element_text(size = 14, hjust = 0.5), legend.position = "right")

#svg("../out/DETEs_mCG_2sections_dodged.svg", width=3.6, height=1.8) #ED Fig. 7d
plot(cl_counts)
#dev.off()

# Perform t-test for each combination of `comparison` and `section`
t_test_results <- meth_section %>%
  group_by(comparison) %>%
  summarise(
    t_test_p_value = t.test(meth ~ section, data = pick(everything()))$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adjusted = p.adjust(t_test_p_value, method = "BH") # Adjust for multiple testing
  )

# Display results
print(t_test_results)

# Wilcoxon rank test
wilcoxon_results <- meth_section %>%
  group_by(comparison) %>%
  summarise(
    wilc_p_value = wilcox.test(meth ~ section, data = pick(everything()))$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adjusted = p.adjust(wilc_p_value, method = "BH") # Adjust for multiple testing
  )

print(wilcoxon_results)
#write_tsv(wilcoxon_results, paste0("../out/Wilcoxon_DETEs_m", context,"_2sections_section.tsv"))

#Clean up
rm(list=ls(pattern="^meth.s"))
rm(list=ls(pattern="^myobj.s"))
rm(list=ls(pattern="^pooled.meth.s"))


