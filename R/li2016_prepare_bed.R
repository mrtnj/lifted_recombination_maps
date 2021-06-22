## Prepare Dryad table from Li et al. (2016) for lifting

library(dplyr)
library(purrr)
library(readr)


li2016 <- read_delim("data/dryad_data_cattle_rmap_Li.txt",
                     delim = " ")

li2016_chr <- split(li2016, li2016$Chr)

haldane_cM <- function(r) -50 * log(1 - 2 * r)

get_marker_itervals <- function(chr) {
    
    window_start <- chr$Location[-nrow(chr)]
    window_end <- chr$Location[-1]
    window_length_bp <- window_end - window_start
    average_r <- (chr$map_f[-nrow(chr)] + chr$map_m[-nrow(chr)])/2
    window_length_cM <- haldane_cM(average_r)
    
    data.frame(chr = unique(chr$Chr),
               start = window_start,
               end = window_end,
               window_length_bp = window_length_bp,
               window_length_cM = window_length_cM,
               rec_rate = window_length_cM/window_length_bp * 1e6)
}


windows <- map_dfr(li2016_chr, get_marker_itervals)


