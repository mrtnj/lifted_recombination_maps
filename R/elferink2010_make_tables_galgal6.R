## Combine lifted positions with cM positions from Elferink et al. (2010)

library(assertthat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(purrr)
library(readr)

source("R/helper_functions.R")

lifted <- read_bed("outputs/elferink2010_hglift_galgal6.bed")

positions <- read_tsv("outputs/elferink2010_marker_positions.txt",
                      col_types = "ccncnc")


combined <- inner_join(lifted, positions[, c("marker", "position_cM", "position_bp", "chr_linkage_map")],
                       by = c("name" = "marker"))


## Remove markers that have switched chromsomes

combined <- filter(combined,
                   chr == sub(chr_linkage_map, pattern = "GGA|GG", replacement = "chr"))


plot_lifted_positions <- qplot(x = end, y = position_cM,
                               data = combined) +
    facet_wrap(~ chr, scale = "free")



## Check the order and distance on the lifted map

lifted_map_positions <- data.frame(chr = combined$chr,
                                   position_bp = combined$end,
                                   position_bp_old = combined$position_bp,
                                   name = combined$name,
                                   position_cM = combined$position_cM,
                                   stringsAsFactors = FALSE)

lifted_map_positions <- lifted_map_positions[order(lifted_map_positions$chr,
                                                   lifted_map_positions$position_bp),]

lifted_map_chr <- split(lifted_map_positions, lifted_map_positions$chr)


## Number of markers lifted

n_markers_lifted <- sum(map_dbl(lifted_map_chr, nrow))


## Filter by number of disagreements per marker

lifted_map_chr_filter1 <- iterative_filter_marker_flip(lifted_map_chr)


## Number of markers after in first filter

n_markers_filter1 <- sum(map_dbl(lifted_map_chr_filter1, nrow))



## Remove problematic chr16

lifted_map_chr_filter1 <- lifted_map_chr_filter1[names(lifted_map_chr_filter1) != "chr16"]





## Create windows from retained lifted markers

windows_lifted <- map_dfr(lifted_map_chr_filter1,
                          function(on_chr) {
                              get_windows_chr(chr = on_chr$chr,
                                              position_bp = on_chr$position_bp,
                                              position_cM = on_chr$position_cM,
                                              marker_name = on_chr$name)
                          })




plot_windows_lifted_diagnostic <- qplot(x = window_length_bp,
                                        y = window_length_cM,
                                        data = windows_lifted)

## Filtering of windows

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




## Create map file (single file format)

map_file <- tibble(marker = windows_lifted_filtered$window_end_marker,
                   chr = windows_lifted_filtered$chr,
                   position_bp = windows_lifted_filtered$end,
                   position_cM = windows_lifted_filtered$window_end_position_cM)

write.table(map_file,
            file = "outputs/elferink2010_GRCg6a.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)





## Create map file (SHAPEIT2 format)

windows_lifted_filtered_chr <- split(windows_lifted_filtered,
                                     windows_lifted_filtered$chr)


make_shapeit_map <- function(windows) {
    tibble(pos = windows$end,
           rate = windows$rec_rate,
           cM = windows$window_end_position_cM)
}

system("mkdir outputs/elferink2010_shapeit_GRCg6a")

for (chr_ix in 1:length(windows_lifted_filtered_chr)) {
    write.table(make_shapeit_map(windows_lifted_filtered_chr[[chr_ix]]),
                file = paste("outputs/elferink2010_shapeit_GRCg6a/",
                             names(windows_lifted_filtered_chr)[chr_ix],
                             ".txt",
                             sep = ""),
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE,
                sep = "\t")
}
