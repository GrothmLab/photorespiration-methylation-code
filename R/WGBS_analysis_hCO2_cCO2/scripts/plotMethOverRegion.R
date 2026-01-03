library(methylKit)
library(genomation)
library(ggplot2)
library(reshape2)
library(TxDb.Athaliana.BioMart.plantsmart28)

plotMethOverRegion <- function(data, annotation="genes", context="CG", windowsize=100, insidewindows=20, range=1000, plot=TRUE, plot_title, clrs=clrplt) {
  
  y_limits <- c()
  if (class(annotation) == 'GRanges') {
    message("Using annotation:") 
    print(annotation)
  } else if (annotation == "genes") {
    message("Using annotation: genes")
    genes <- genes(TxDb.Athaliana.BioMart.plantsmart28)
    seqlevels(genes) <- c("1","2","3","4","5","Pt","Mt")
    annotation <- genes
    #Adjust y_limits to region and context 
    plot_title <- "Genes"
    if (context == "CG") {
      y_limits <- c(0,0.4)
    } else if (context == "CHG") {
      y_limits <- c(0,0.02)
    } else if (context == "CHH") {
      y_limits <- c(0,0.01)
    }
  } else if (annotation == "TEs") {
    message("Using annotation: TEs")
    load("../source-data/TEs.RData")
    annotation <- TEs
    plot_title <- "TEs"
    if (context == "CG") {
      y_limits <- c(0,1)
    } else if (context == "CHG") {
      y_limits <- c(0,0.4)
    } else if (context == "CHH") {
      y_limits <- c(0,0.1)
    }
  } else {
    stop("Argument 'annotation' needs to be of class 'GRanges'.")
  }  
  
  if(missing(data)) {
    stop("No data provided!\nUsage: plotMethOverRegion(data, annotation='genes', context='CG', windowsize=100, insidewindows=20, range=1000, plot=TRUE, clrs=c('black','grey','blue','green','red','orange'))")
  }
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
  
  overlaps.meth.lvl <- array(NA, dim=c(length(data), 1, range%/%windowsize, 2, 2), dimnames=list(sample=names(data), category=context, distance=1:(range%/%windowsize), what=c('mean', 'weight'), where=c('upstream', 'downstream')))
  ioverlaps.meth.lvl <- array(NA, dim=c(length(data), 1, insidewindows, 2, 1), dimnames=list(sample=names(data), category=context, distance=1:insidewindows, what=c('mean', 'weight'), where=c('inside')))
  
  for (i in 1:length(data)){
    
    if (class(data[[i]]) == 'GRanges') {
      message("Calculating upstream, downstream and inside counts for ", names(data)[[i]])
    } else {
      stop("Argument 'data' needs to be of class 'GRanges'.")
    }
    
    ## Subset annotation to data range
    annotation <- subsetByOverlaps(annotation, data[[i]], ignore.strand=TRUE)
    
    ## Upstream and downstream annotation
    annotation.up <- resize(x = annotation, width = 1, fix = 'start')
    annotation.down <- resize(x = annotation, width = 1, fix = 'end')
    
    for (i1 in 1:(range %/% windowsize)) {
      anno.up <- suppressWarnings( promoters(x = annotation.up, upstream = i1*windowsize, downstream = 0) )
      anno.up <- suppressWarnings( resize(x = anno.up, width = windowsize, fix = 'start') )
      ind <- findOverlaps(anno.up, data[[i]], ignore.strand=TRUE)
      overlaps.meth.lvl[names(data)[[i]], category, as.character(i1), 'mean', 'upstream'] <- (sum(data[[i]]$numCs[ind@to], na.rm=TRUE)) / (sum(data[[i]]$coverage[ind@to], na.rm=TRUE))
    }
    for (i1 in 1:(range %/% windowsize)) {
      anno.down <- suppressWarnings( promoters(x = annotation.down, upstream = 0, downstream = i1*windowsize) )
      anno.down <- suppressWarnings( resize(x = anno.down, width = windowsize, fix = 'end') )
      ind <- findOverlaps(anno.down, data[[i]], ignore.strand=TRUE)
      overlaps.meth.lvl[names(data)[[i]], category, as.character(i1), 'mean', 'downstream'] <- (sum(data[[i]]$numCs[ind@to], na.rm=TRUE)) / (sum(data[[i]]$coverage[ind@to], na.rm=TRUE))
    }
    
    ## Inside annotation
    
    mask.plus <- as.logical(strand(annotation) == '+' | strand(annotation) == '*')
    mask.minus <- !mask.plus
    widths <- width(annotation)
    
    
    for (i1 in 1:insidewindows) {
      anno.in <- annotation
      end(anno.in[mask.plus]) <- start(annotation[mask.plus]) + round(widths[mask.plus] * i1 / insidewindows)
      start(anno.in[mask.plus]) <- start(annotation[mask.plus]) + round(widths[mask.plus] * (i1-1) / insidewindows)
      start(anno.in[mask.minus]) <- end(annotation[mask.minus]) - round(widths[mask.minus] * i1 / insidewindows)
      end(anno.in[mask.minus]) <- end(annotation[mask.minus]) - round(widths[mask.minus] * (i1-1) / insidewindows)
      ind <- findOverlaps(anno.in, data[[i]], ignore.strand=TRUE)
      
      ioverlaps.meth.lvl[names(data)[[i]], category, as.character(i1), 'mean', 'inside'] <- (sum(data[[i]]$numCs[ind@to], na.rm=TRUE)) / (sum(data[[i]]$coverage[ind@to], na.rm=TRUE))
    }
    
  }
  
  ## Prepare data.frames for plotting
  dfo <- reshape2::melt(overlaps.meth.lvl)
  df <- dfo[dfo$what == 'mean',]
  names(df)[6] <- 'mean'
  #df$weight <- dfo$value[dfo$what == 'weight']
  df$distance <- as.numeric(df$distance)
  
  # Adjust distances for plot
  df$distance[df$where == 'upstream'] <- -df$distance[df$where == 'upstream'] * windowsize + windowsize/2
  df$distance[df$where == 'downstream'] <- df$distance[df$where == 'downstream'] * windowsize + range - range/insidewindows - windowsize/2
  
  dfo <- reshape2::melt(ioverlaps.meth.lvl)
  idf <- dfo[dfo$what == 'mean',]
  names(idf)[6] <- 'mean'
  #idf$weight <- dfo$value[dfo$what == 'weight']
  idf$distance <- as.numeric(idf$distance)
  
  # Adjust distances for plot
  idf$distance <- (idf$distance -1) * range / insidewindows
  
  df.plot <- rbind(df, idf)
  
  if (!plot) {
    return(df.plot)
  }
  
  # Enrichment profile
  ggplt <- ggplot(df.plot) + geom_vline(xintercept = 0, color = "grey80")
  ggplt <- ggplt + geom_vline(xintercept = ((1-1/insidewindows)*range), color = "grey80")
  ggplt <- ggplt + geom_line(aes_string(x='distance', y='mean', col='sample'), linewidth=0.3)
  ggplt <- ggplt + geom_line(aes_string(x='distance', y='mean', col='sample'), stat="smooth", method=NULL, span=0.2, se=FALSE, linewidth=1, alpha=0.5)
  ggplt <- ggplt + theme_classic()
  ggplt <- ggplt + scale_colour_manual(values=clrs)
  ggplt <- ggplt + scale_x_continuous(breaks=breaks, labels=labels)
  ggplt <- ggplt + coord_cartesian(ylim=y_limits, expand = TRUE, default = TRUE)
  ggplt <- ggplt + ggtitle(plot_title)
  ggplt <- ggplt + theme(axis.title = element_text(size = 13), axis.text = element_text(size = 12), axis.text.x = element_text(angle=45, hjust = 1), plot.title = element_text(size = 14, hjust = 0.5), legend.text = element_text(size = 12), legend.title = element_text(size = 13))
  ggplt <- ggplt + xlab('')
  ggplt <- ggplt + ylab(paste0('Mean ', context, 'm'))
  
  # Cytosine density
  #ggpltw <- ggplot(df.plot) + geom_line(aes_string(x='distance', y='weight', col='category'))
  #ggpltw <- ggpltw + scale_x_continuous(breaks=breaks, labels=labels)
  #ggpltw <- ggpltw + xlab('Distance from annotation') + ylab('Normalized C frequency\nper bp of annotation')
  
  #cowplt <- cowplot::plot_grid(c(ggplt, ggpltw), align = 'v', nrow=2, rel_heights=c(2,1))
  
  if (plot) {
    return(ggplt)
  }
}

