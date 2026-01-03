#Set working directory
#setwd("/path/to/photorespiration_methylation_code/R/WGBS_analysis_m_t_mt_wt/scripts")
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

library(readr)
library(genomation)
library(GenomicFeatures)
library(rtracklayer)
library(methylKit)
library(pheatmap)
library(dendextend)
library(ComplexHeatmap)
library(svglite)
library(tidyverse)
library(dplyr)
library(RColorBrewer)
library(reshape2)

samplesheet <- read_csv("../samplesheet.csv")
(sample.id <- as.list(sort(samplesheet$sample)))
(coldata <- read_csv("../WGBS_mt_sample_info.csv"))

#Line colors for metaplots
clrplt4 <- c("indianred3", "peru", "turquoise4", "steelblue4")

#Line colors for heatmap
#col.ramp <- colorRampPalette(brewer.pal(9, "Reds") )(255)
colors_or <- brewer.pal(9,"OrRd")
colors_r <- brewer.pal(9,"Reds")
colors_rb <- rev(brewer.pal(11,"RdBu"))
colors_rg <- rev(brewer.pal(11,"RdGy"))
colors_ygb <- rev(brewer.pal(9,"YlGnBu"))
pal_or <- colorRampPalette(colors_or)
pal_r <- colorRampPalette(colors_r)
pal_rb <- colorRampPalette(colors_rb)
pal_rg <- colorRampPalette(colors_rg)
pal_ygb <- colorRampPalette(colors_ygb)

#Function to reposition brakes at the quantiles for heatmaps
quantile_breaks <- function(xs, n = 100) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n), na.rm = TRUE)
  breaks[!duplicated(breaks)]
}

#Read RNAseq count data, i.e.normalized counts (log2(counts+1)) from DESeq2 
counts <- read_rds("../source-data/m_t_mt_DETEs_ntd_counts.rds")
DETEs <- rownames(counts)

###Make GRanges object for TEs
TE_gr <- import("../source-data/Panda_AtTEs_transcripts_v1_ensembl_fixed.gtf")
TE_txdb <- makeTxDbFromGRanges(TE_gr)
TE_gr <- genes(TE_txdb)

#Subset TE_gr to DETEs and merge with tx2gene to have metadata with tx ids
DETE_gr <- TE_gr[TE_gr$gene_id %in% DETEs,]

#Make sure DETE_gr is sorted by coordinates
#DETE_gr_sort <- GenomicRanges::sort(DETE_gr) #I don't knkow why it's not sorting properly! Therefore doing this:
DETE_gr_df <- as.data.frame(DETE_gr)
DETE_gr_sort_df <- DETE_gr_df[order(DETE_gr_df$seqnames, DETE_gr_df$start, DETE_gr_df$end), ]
DETE_gr_sorted <- GenomicRanges::makeGRangesFromDataFrame(DETE_gr_sort_df, keep.extra.columns=TRUE)

#####
#Calculated pooled methylation ratios for DETEs

context <- "CG"

(file.list <- as.list(list.files(path = "../source-data/methylkit_in", pattern = "CG_report_sorted\\.txt$", full.names = TRUE)))

#Generate methylkit object for m, mt, t and wt
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
#Pool replicates per genotype
myobj.s.m <- reorganize(myobj.s,
                        c("m_1","m_2"),
                        treatment=c(0,0))
meth.s.m <- methylKit::unite(myobj.s.m, destrand=FALSE, min.per.group=1L, mc.cores=2)
pooled.meth.s.m <- pool(meth.s.m,sample.ids = "m")

myobj.s.mt <- reorganize(myobj.s,
                         c("mt_1","mt_2"),
                         treatment=c(0,0))
meth.s.mt <- methylKit::unite(myobj.s.mt, destrand=FALSE, min.per.group=1L, mc.cores=2)
pooled.meth.s.mt <- pool(meth.s.mt,sample.ids = "mt")

myobj.s.t <- reorganize(myobj.s,
                        c("t_1","t_2"),
                        treatment=c(0,0))
meth.s.t <- methylKit::unite(myobj.s.t, destrand=FALSE, min.per.group=1L, mc.cores=2)
pooled.meth.s.t <- pool(meth.s.t,sample.ids = "t")

myobj.s.wt <- reorganize(myobj.s,
                         c("WT_1","WT_2"),
                         treatment=c(0,0))
meth.s.wt <- methylKit::unite(myobj.s.wt, destrand=FALSE, min.per.group=1L, mc.cores=2)
pooled.meth.s.wt <- pool(meth.s.wt,sample.ids = "wt")

pooled.meth.s <- c()
pooled.meth.s <- list(pooled.meth.s.m, pooled.meth.s.mt, pooled.meth.s.t, pooled.meth.s.wt)

#Convert to GRanges objects for metaplotting
myobj.grl <- c()
for (i1 in 1:length(pooled.meth.s)) {
  myobj.grl <- c(myobj.grl, as(pooled.meth.s[[i1]],"GRanges"))
  colnames(mcols(myobj.grl[[i1]])) <- c("coverage","numCs","numTs")
}
names(myobj.grl) <- c("m","mt","t","wt")

#####
### Calculate methylation in windows over DETEs
annotation <- DETE_gr_sorted
data <- myobj.grl
windowsize=500; insidewindows=4; range=1000

#Subset data to annotated regions +/- range
#Firts convert data from list to GRangesList
data <- GRangesList(data)

#Remove annotations that do not overlap with data
for (i in 1:length(data)){
  annotation <- subsetByOverlaps(annotation, data[[i]], ignore.strand=TRUE, minoverlap = 1L) #minoverlap = length(data) makes sure that all sapmles overlap with data
}

#Resize annotation
data_range <- resize(x = annotation, width = width(annotation) + range, fix = "start")
data_range <- resize(x = data_range, width = width(data_range) + range, fix = "end")
data_range <- sort(data_range)

#Subset data
for (i in 1:length(data)){
  data[[i]] <- subsetByOverlaps(data[[i]], data_range, ignore.strand=TRUE)
}

breaks <- c(c(-1, 0, (1-1/insidewindows), (1-1/insidewindows+1)) * range)

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

## Prepare data.frames for combining with count data
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
### Compute methylation over DETEs
myobj.s.DETEs <- regionCounts(myobj.s, DETE_gr_sorted, cov.bases=4, mc.cores=2, save.db=FALSE)
#myobj.s.DETEs.filt <- filterByCoverage(myobj.s.DETEs, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.95)

#Select windows covered in all samples
meth.s.DETEs <- methylKit::unite(myobj.s.DETEs, destrand=FALSE, min.per.group=1L, mc.cores=2)

#Get methylated DETEs names
DETE_gr_sorted_sub <- subsetByOverlaps(DETE_gr_sorted, as(meth.s.DETEs, "GenomicRanges"))
DETEs_sub <- names(DETE_gr_sorted_sub)

#Remove DETEs  that are not covered by methylation reads from count tables
counts_sub <- subset(counts, rownames(counts) %in% DETEs_sub)

#Add DETE IDs as rownames
rownames(meth.s.DETEs) <- names(DETE_gr_sorted_sub)

#Pool replicates
pooled.meth.s.DETEs <- pool(meth.s.DETEs, sample.ids = c("m","mt","t","wt"))

#Order pooled.meth.s.DETEs based on expression mat
ordered_ind <- match(rownames(counts_sub), rownames(pooled.meth.s.DETEs))
pooled.meth.s.DETEs <- pooled.meth.s.DETEs[ordered_ind, , drop = FALSE]

#Calculate pooled CG methylation over DETEs
pooled.meth.s.DETEs.mat <- percMethylation(pooled.meth.s.DETEs, rowids=TRUE)
rownames(pooled.meth.s.DETEs.mat) <- rownames(pooled.meth.s.DETEs)

#####
### Load RNAseq sample info 
coldata_rnaseq <- read_tsv("../mtRNAseq_sampleinfo.tsv")

### Load expression data for plotting together with DNA methylation
### Used TPM data for plotting
#Load tpm data
tpm_df <- read_tsv("../source-data/salmon.merged.gene_tpm.tsv") %>% 
  filter(gene_id %in% DETEs_sub) #limit to DETEs

#Subset and transform data
tpm_df_mean_longer <- pivot_longer(tpm_df[,-2], cols = m_2_LD:wt_4_LD, names_to = "sample_name", values_to = "tpm", values_drop_na=FALSE) %>%
  # join with sample info table
  full_join(coldata_rnaseq, by = ("sample_name")) %>%
  # for each gene and genotype
  group_by(gene_id, condition) %>%
  # mean and normalized mean
  summarise(mean_tpm = mean(tpm),
            nt_mean_tpm = log2(mean_tpm + 1),
            nrep = n()) %>%
  rename(sample = condition) %>%
  mutate(sample = str_replace(sample,"\\_.*","")) %>%
  ungroup()

tpm_df_mean <- dplyr::select(tpm_df_mean_longer, c(gene_id, sample, nt_mean_tpm)) %>%
  pivot_wider(names_from = sample, values_from = nt_mean_tpm)
tpm_matrix_mean <- as.matrix(tpm_df_mean[,2:5])
rownames(tpm_matrix_mean) <- tpm_df_mean$gene_id

#Order rows according to counts, so that row clustering fits for all heatmaps
tpm_matrix_mean <- tpm_matrix_mean[match(rownames(counts_sub), rownames(tpm_matrix_mean)),]

#####
### Draw heatmaps
### Use count matrix for clustering the heatmaps
hclust_matrix <- counts_sub
hclust_matrix_mean <- as.data.frame(hclust_matrix)
hclust_matrix_mean$gene <- rownames(hclust_matrix_mean)

#Calculate mean counts
hclust_matrix_mean <- pivot_longer(hclust_matrix_mean, cols = m_2_LD:wt_4_LD, names_to = "sample_name", values_to = "cts")  %>% 
  # join with sample info table
  full_join(coldata_rnaseq, by = ("sample_name")) %>%
  # for each gene and genotype
  group_by(gene, condition) %>%
  # calculate the mean cts
  summarise(mean_cts = mean(cts),
            nrep = n()) %>% 
  ungroup()

hclust_matrix_mean <- pivot_wider(hclust_matrix_mean, names_from = condition, values_from = mean_cts)
hclust_matrix_mean_mat <- as.matrix(hclust_matrix_mean[,3:6])
rownames(hclust_matrix_mean_mat) <- hclust_matrix_mean$gene
hclust_matrix_mean_mat <- as.matrix(hclust_matrix_mean_mat)

#Order rows according to counts_sub, so that row clustering fits for all heatmaps
hclust_matrix_mean_mat <- hclust_matrix_mean_mat[match(rownames(counts_sub), rownames(hclust_matrix_mean_mat)),]

#Scale mean counts
hclust_matrix_mean_mat <- hclust_matrix_mean_mat %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()

# Pairwise pearson correlation between genes
#gene_cor <- cor(t(hclust_matrix_mean_mat), use = "pairwise.complete.obs", method = "pearson")
#gene_dist <- as.dist(1-gene_cor)
#gene_hclust <- hclust(gene_dist, method="complete")

#Euclidean distance between genes
gene_dist <- dist(hclust_matrix_mean_mat, method = "euclidean") 
gene_hclust <- hclust(gene_dist, method="ward.D2") #"ward.D2"

#plot(gene_hclust, labels = FALSE)
#abline(h = 8.3, col = "red", lwd = 1)
#abline(h = 8.2, col = "orange", lwd = 1)

sample_dist <- dist( t( tpm_matrix_mean ) )
sample_hclust <- hclust(sample_dist, method="complete")
plot(sample_hclust)

df <- data.frame(genotype = c("m","mt","t","wt"))
rownames(df) <- c("m","mt","t","wt")

#Set legend and colors
genotype_colors <- clrplt4
df_colors <- list(genotype = genotype_colors)
names(df_colors$genotype) <- c("m","mt","t","wt")

mat_breaks <- quantile_breaks(tpm_matrix_mean, n = 101)

pm <- pheatmap::pheatmap(tpm_matrix_mean,
               cutree_rows = 3,
               cluster_rows=gene_hclust,
               cluster_cols=sample_hclust,
               show_rownames=FALSE,
               annotation_col=df,
               annotation_colors=df_colors,
               treeheight_row=16,
               treeheight_col=8,
               color=pal_rb(length(mat_breaks) - 1),
               border_color=NA,
               breaks=mat_breaks,
               fontsize=7,
               fontsize_col=4
)

row_order <- pm$tree_row$order

#Set legend and colors
genotype_colors <- clrplt4
df_colors <- list(genotype = genotype_colors)
names(df_colors$genotype) <- c("m","mt","t","wt")

cpm <- ComplexHeatmap::pheatmap(tpm_matrix_mean,
                                name = "log2(TPM+1)",
                                cutree_rows = 3,
                                cluster_rows=gene_hclust,
                                cluster_cols=sample_hclust,
                                show_rownames=FALSE,
                                annotation_col=df,
                                annotation_names_col = FALSE,
                                annotation_colors=df_colors,
                                treeheight_row=16,
                                treeheight_col=16,
                                color=pal_ygb(length(mat_breaks) - 1),
                                border_color=NA,
                                breaks=mat_breaks,
                                fontsize=7,
                                fontsize_col=6,
                                cellwidth=8,
                                heatmap_legend_param=list(legend_direction = "horizontal", legend_width = unit(3, "cm"), annotation_legend_side = "bottom")
)            
draw(cpm, heatmap_legend_side="bottom", annotation_legend_side = "bottom")

ht_cpm <- draw(cpm, heatmap_legend_side="bottom", annotation_legend_side = "bottom")
row_dend_cpm <- row_dend(ht_cpm)

###Get the cluster order for splitting the heatmaps
ht <- draw(cpm)
ro.list <- row_order(ht)

library(magrittr)
sum(unlist(lapply(ro.list, function(x) length(x))))
length(gene_hclust$labels)

###Get the cluster order for splitting the heatmaps
clu.df <- lapply(1:length(ro.list), function(i) {
  if (i < 10) {
    out <- data.frame(GeneID = gene_hclust$labels[ro.list[[i]]],
                      Cluster = paste0("c", i),
                      stringsAsFactors = FALSE)
    return(out)
  } else {
    out <- data.frame(GeneID = gene_hclust$labels[ro.list[[i]]],
                      Cluster = paste0("c", i),
                      stringsAsFactors = FALSE)
    return(out)
  }
}) %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu.df'
  do.call(rbind, .)

rownames(clu.df) <- clu.df$GeneID
clu.df <- clu.df[gene_hclust$labels,]
clu.df$Order <- NA
clu.df$Order <- gene_hclust$order
#clu.df <- clu.df[names(DD_genes),]
split <-  clu.df$Cluster
row_order <- clu.df$Order
lgd = Legend(at = c("c1", "c2", "c3"), title = "Clusters", 
             type = "lines", legend_gp = gpar(col = brewer.pal(n=3, name="Dark2"), lwd = c(1.5,1.5,1.5)))

### Draw methylation heatmap
#####

### Draw mean meth per DETE
# meth mat is pooled.meth.s.DETEs.mat
sample_dist <- dist( t( pooled.meth.s.DETEs.mat ) )
sample_hclust <- hclust(sample_dist, method="complete")
plot(sample_hclust)

df <- data.frame(genotype = c("m","mt","t","wt"))
rownames(df) <- c("m","mt","t","wt")

#Set legend and colors
df_colors <- list(genotype = genotype_colors)
names(df_colors$genotype) <- c("m","mt","t","wt")

colors_or <- brewer.pal(9,"OrRd")
pal_or <- colorRampPalette(colors_or)

mat_breaks <- quantile_breaks(pooled.meth.s.DETEs.mat, n = 101)

cpm2 <- ComplexHeatmap::pheatmap(pooled.meth.s.DETEs.mat,
                                 name = "% mCG",
                                 cutree_rows = 3,
                                 cluster_rows=gene_hclust,
                                 treeheight_row=0,
                                 treeheight_col=16,
                                 cluster_cols=sample_hclust,
                                 color = pal_or(length(mat_breaks) - 1),
                                 na_col = "black",
                                 show_rownames = FALSE,
                                 show_colnames = TRUE,
                                 annotation_col=df,
                                 annotation_names_col = FALSE,
                                 annotation_colors=df_colors,
                                 fontsize_col=6,
                                 fontsize = 7,
                                 border_color=NA,
                                 breaks=mat_breaks,
                                 cellwidth=8,
                                 heatmap_legend_param=list(legend_direction = "horizontal", legend_width = unit(2, "cm"), annotation_legend_side = "bottom")
)            
draw(cpm2, heatmap_legend_side="bottom", annotation_legend_side = "bottom")

### Compute correlation between normalized tpm values and DNA methylation
#####
#Calculate correlation between DNA methylation and expression for each window
#Join DNA methylation with expression dataframe
df_all <- merge(dfm, tpm_df_mean_longer, by = c("gene_id", "sample"), all.y = TRUE, sort = FALSE)
df_all <- df_all %>%
  group_by(gene_id) %>%
  arrange(distance, sample, .by_group = TRUE) %>%
  dplyr::select(-category, -where, -mean_tpm, -nrep) %>%
  ungroup()

#Calculate Pearson correlations
cr <- df_all %>%
  group_by(gene_id, distance) %>%
  summarise(spearman_correlation = cor(as.numeric(value), as.numeric(nt_mean_tpm), method = "spearman")) %>% #use = "pairwise.complete.obs"
  ungroup()
cr_w <- pivot_wider(cr, names_from = c(distance), values_from = spearman_correlation)
cr_w <- as.data.frame(cr_w)
rownames(cr_w) <- cr_w$gene_id
cr_w <- as.matrix(cr_w[,-1])

#Order rows according to counts, so that row clustering fits for all heatmaps
cr_w <- cr_w[match(rownames(counts_sub), rownames(cr_w)),]

### Draw correlation heatmap
#####
quantile_breaks <- function(xs, n = 100) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n), na.rm = TRUE)
  breaks[!duplicated(breaks)]
}
mat_breaks <- quantile_breaks(cr_w, n = 101)

colnames(cr_w) <- c("-1kb","","TSS","","","TTS","","1kb")

colors_sp <- brewer.pal(11,"Spectral")
pal_sp <- colorRampPalette(colors_sp)
colors_rb <- rev(brewer.pal(9,"YlGnBu"))
pal_rb <- colorRampPalette(colors_rb)

### Calculate mean correlation per columns and cluster
clusters <- clu.df$Cluster

# Calculate mean values for each cluster
mean_cor_values <- t(apply(cr_w, 2, function(col) tapply(col, clusters, mean, na.rm = TRUE)))

# Create the line plot annotation
line_plot_anno <- anno_lines(mean_cor_values,
                             which = "column",
                             gp = gpar(col = brewer.pal(n=3, name="Dark2"), lwd = c(1.5,1.5,1.5)),
                             ylim=c(-1,1),
                             add_points = FALSE)

cpm3 <- ComplexHeatmap::pheatmap(cr_w,
                                 name = "PCC",
                                 cutree_rows = 3,
                                 cluster_rows=gene_hclust,
                                 treeheight_row=0,
                                 cluster_cols=FALSE,
                                 color = pal_sp(length(mat_breaks)),
                                 na_col = "black",
                                 show_rownames = FALSE,
                                 show_colnames = TRUE,
                                 annotation_names_col = FALSE,
                                 fontsize_col=7,
                                 fontsize = 7,
                                 border_color=NA,
                                 breaks=mat_breaks,
                                 cellwidth=4,
                                 heatmap_legend_param=list(legend_direction = "horizontal", legend_width = unit(3, "cm"), annotation_legend_side = "bottom"),
                                 top_annotation = HeatmapAnnotation(clusters = line_plot_anno)
)
draw(cpm3, heatmap_legend_side="bottom", annotation_legend_side = "bottom")

# Make heatmaps -----------------------------------------------------------
ht_list <- Heatmap(split, col = structure(brewer.pal(n=3, name="Dark2")), name = "cluster", 
                   show_row_names = FALSE, show_heatmap_legend = FALSE, width = unit(2, "mm")) + 
           cpm + cpm2 + cpm3

ht_list <- draw(ht_list, main_heatmap = "log2(TPM+1)", show_row_dend = TRUE, annotation_legend_list = list(lgd), heatmap_legend_side = "bottom", ht_gap = unit(c(2,2,6), "mm"))

#svglite("../out/CHM_DETEs_1_05_ntTPM_CG_meth_corr.svg", width=3, height=5) #Fig. 1e
draw(ht_list, row_split = split, row_order = row_order, annotation_legend_list = list(lgd), heatmap_legend_side = "bottom", ht_gap = unit(c(2,2,6), "mm"))
#dev.off()

#Clean up
rm(list=ls(pattern="myobj"))
rm(list=ls(pattern="meth"))
rm(list=ls(pattern="CG"))
rm(list=ls(pattern="cluster"))
rm(list=ls(pattern="pooled"))
gc()

