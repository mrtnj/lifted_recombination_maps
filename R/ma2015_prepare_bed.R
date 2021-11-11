## Prepare Dryad table from Ma et al. (2015) for lifting

library(dplyr)
library(purrr)
library(readr)

source("R/helper_functions.R")


ma2015 <- read_delim("data/dryad_data_cattle_rmap_Ma2015.txt",
                     delim = " ")




## Marker positions

bed <- data.frame(chr = paste("chr", ma2015$Chr, sep = ""),
                  start = ma2015$Location - 1,
                  end = ma2015$Location,
                  name = ma2015$Name,
                  stringsAsFactors = FALSE)


write.table(bed,
            file = "outputs/ma2015_windows.bed",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)





## Marker intervals

ma2015_chr <- split(ma2015, ma2015$Chr)

get_marker_itervals <- function(chr) {
    
    window_start <- chr$Location[-nrow(chr)]
    window_end <- chr$Location[-1]
    window_end_marker <- chr$Name[-1]
    window_length_bp <- window_end - window_start
    average_r <- (chr$map_f[-nrow(chr)] + chr$map_m[-nrow(chr)])/2
    window_length_cM <- haldane_cM(average_r)
    window_end_position_cM <- cumsum(window_length_cM)
    
    data.frame(chr = paste("chr", unique(chr$Chr), sep = ""),
               start = window_start,
               end = window_end,
               window_length_bp = window_length_bp,
               window_length_cM = window_length_cM,
               window_end_marker,
               window_end_position_cM,
               rec_rate = window_length_cM/window_length_bp * 1e6)
}


windows <- map_dfr(ma2015_chr, get_marker_itervals)

windows$window_id <- make_window_id(windows)

write.table(windows,
            file = "outputs/ma2015_windows.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)


windows_bed <- data.frame(chr = windows$chr,
                          start = windows$start - 1,
                          end = windows$end,
                          name = windows$window_id)


write.table(windows_bed,
            file = "outputs/ma2015_windows.bed",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)