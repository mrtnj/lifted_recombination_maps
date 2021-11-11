## Combine lifted positions with cM positions from Ma et al. (2015)

library(assertthat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(purrr)
library(readr)

source("R/helper_functions.R")

lifted <- read_bed("outputs/ma2015_windows_hglift_ars-ucd1.2_0.6remap.bed")

windows <- read_tsv("outputs/ma2015_windows.txt")

colnames(windows)[1:3] <- paste(colnames(windows)[1:3], "_old", sep = "")



combined <- inner_join(lifted,
                       windows,
                       by = c("name" = "window_id"))



plot_lifted_positions <- qplot(x = end, y = window_end_position_cM,
                               data = combined) +
    facet_wrap(~ chr, scale = "free")



## Check the order and distance on the lifted map

lifted_map_positions <- data.frame(chr = combined$chr,
                                   position_bp = combined$end,
                                   position_bp_old = combined$end_old,
                                   name = combined$window_end_marker,
                                   position_cM = combined$window_end_position_cM,
                                   stringsAsFactors = FALSE)

lifted_map_positions <- filter(lifted_map_positions,
                               chr %in% paste("chr", 1:29, sep = ""))

lifted_map_positions <- lifted_map_positions[order(lifted_map_positions$chr,
                                                   lifted_map_positions$position_bp),]

lifted_map_chr <- split(lifted_map_positions, lifted_map_positions$chr)



## Filter map by removing markers that have flipped cM order relative to many
## other markers

is_cM_ordered <- function(map_chr) {
    
    genetic_orders <- map(map_chr,
                          function(chr) order(chr$position_cM))
    
    is_ordered <- map_lgl(genetic_orders,
                          function(o) identical(o, 1:length(o)))
    
    all(is_ordered)
}


## Count the number of other markers that are on the "wrong side" of each marker
## based on their genetic positions; this gives a score of how many other markers
## disagree with a particular marker.

get_marker_flip_count_chr <- function(chr) {
    
    n_markers <- nrow(chr)
    n_flipped <- numeric(n_markers)
    
    for (marker_ix in 1:n_markers) {
        distance <- chr$position_cM - chr$position_cM[marker_ix]
        n_flipped_before <- sum(distance[1:marker_ix] > 0)
        n_flipped_after <- sum(distance[(marker_ix + 1):n_markers] < 0)
        if (is.na(n_flipped_after)) {
            n_flipped_after <- 0
        }
        n_flipped[marker_ix] <- n_flipped_before + n_flipped_after
    }
    n_flipped
}


## Filter a map by removing, for each chromosome, the markers that have the
## highest number of disagreements with other markers. If several have the same
## number all of those are removed.

filter_marker_flip <- function(map_chr) {
    
    flip_counts <- map(map_chr, get_marker_flip_count_chr)
    
    n_chr <- length(map_chr)
    
    for (chr_ix in 1:n_chr) {
        if (any(flip_counts[[chr_ix]] > 0)) {
            max_flip <- max(flip_counts[[chr_ix]])
            map_chr[[chr_ix]] <- map_chr[[chr_ix]][flip_counts[[chr_ix]] < max_flip,]
        }
    }
    map_chr
}


## Iteratively remove the markers that disagree the most until no further
## disagreements about marker order remain.

iterative_filter_marker_flip <- function(map_chr) {

    k <- 1
    
    while (!is_cM_ordered(map_chr)) {
        print(k)
        map_chr <- filter_marker_flip(map_chr)   
        k <- k + 1
    }
    
    map_chr
}

lifted_map_chr_flip_filter <- iterative_filter_marker_flip(lifted_map_chr)


## Check how far each interval has been lifted, compared to the median distance

get_distance_lifted <- function(chr) {
    chr$distance_lifted <- abs(chr$position_bp - chr$position_bp_old)
    chr$median_distance_lifted <- median(chr$distance_lifted)
    chr
}

lifted_map_chr <- map(lifted_map_chr, get_distance_lifted)

lifted_map_chr_filtered <- map(lifted_map_chr,
                               function(chr) {
                                   filter(chr, distance_lifted < median_distance_lifted + 2e6)
                               })

sum(unlist(lapply(lifted_map_chr, nrow)))

sum(unlist(lapply(lifted_map_chr_filtered, nrow)))



## Create windows from retained lifted markers

windows_lifted <- map_dfr(lifted_map_chr_filtered,
                          function(on_chr) {
                              get_windows_chr(chr = on_chr$chr,
                                              position_bp = on_chr$position_bp,
                                              position_cM = on_chr$position_cM,
                                              marker_name = on_chr$name)
                          })




plot_windows_lifted_diagnostic <- qplot(x = window_length_bp,
                                        y = window_length_cM,
                                        data = windows_lifted)

## Filtering

windows_lifted_filtered <- filter(windows_lifted, 
                                  window_length_bp > 0 &
                                      window_length_bp < 1e6 &
                                      window_length_cM >= 0 & 
                                      rec_rate < 100)


## Diagnostic plots

plot_windows_lifted_filtered_diagnostic <- qplot(x = window_length_bp,
                                                 y = window_length_cM,
                                                 data = windows_lifted_filtered)

plot_windows_lifted_diagnostic + plot_windows_lifted_filtered_diagnostic


plot_lifted_filtered_positions <- qplot(x = end, y = window_end_position_cM,
                                        data = windows_lifted_filtered) +
    facet_wrap(~ chr, scale = "free")


plot_lifted_positions + plot_lifted_filtered_positions


## Create map file (SHAPEIT2 format)

windows_lifted_filtered_chr <- split(windows_lifted_filtered,
                                     windows_lifted_filtered$chr)


make_shapeit_map <- function(windows) {
    tibble(pos = windows$end,
           rate = windows$rec_rate,
           cM = windows$window_end_position_cM)
}


for (chr_ix in 1:length(windows_lifted_filtered_chr)) {
    write.table(make_shapeit_map(windows_lifted_filtered_chr[[chr_ix]]),
                file = paste("outputs/ma2015_shapeit/",
                             names(windows_lifted_filtered_chr)[chr_ix],
                             ".txt",
                             sep = ""),
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE,
                sep = "\t")
}
