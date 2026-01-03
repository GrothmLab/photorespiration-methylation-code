library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(stringr)

#Function to count overlaps and return a data frame
count_overlaps <- function(grl, gr) {
  overlaps <- lapply(grl, function(x) sum(countOverlaps(x, gr, ignore.strand=TRUE)))
  data.frame(state = as.integer(names(grl)), overlaps = unlist(overlaps)) %>%
    arrange(state)
}

#Function to name states_grl
statenames <- function(n) {
  basename(n) %>%
    str_remove("\\..*") %>%
    str_remove("State")
}

countCSoverlapsGR <- function(q_gr) {

  #Path to states bed files
  (bed_files <- list.files(path = "/Users/martin.groth/GNMX/Refdata/PCSD/Arabidopsis", pattern = "\\.bed$", full.names = TRUE))
  
  #Read the states bed files into a GRangesList
  states_grl <- lapply(bed_files, import.bed)
  names <- lapply(bed_files, statenames)
  names(states_grl) <- names
  
  #Count overlaps and store results in a data frame
  ol_df <- count_overlaps(states_grl, q_gr)
  return(ol_df)
}  