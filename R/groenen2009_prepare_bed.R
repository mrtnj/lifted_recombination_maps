## Prepare supplementary table from Groenen et al. (2009) for lifting

library(dplyr)
library(purrr)
library(readr)


## Bed file of single marker positions

groenen2009 <- read_tsv("data/Supplemental_Table_S1_Groenen.txt",
                        col_types = "ccnnnnc")[,-7]

colnames(groenen2009) <- c("chromosome", "name", "female_cM", "male_cM", "average_cM",
                           "position_galgal2")

bed <- data.frame(chr = paste("chr", groenen2009$chromosome, sep = ""),
                  start = groenen2009$position_galgal2 - 1,
                  end = groenen2009$position_galgal2,
                  name = groenen2009$name,
                  stringsAsFactors = FALSE)

bed <- filter(bed, !is.na(start))


write.table(bed,
            file = "outputs/groenen2009.bed",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)



## Bed file of marker intervals


groenen2009_chr <- split(groenen2009, groenen2009$chromosome)



get_marker_itervals_sorted <- function(chr) {
    ##chr <- chr[order(chr$position_galgal2),]
    window_start <- chr$position_galgal2[-nrow(chr)]
    window_end <- chr$position_galgal2[-1]
    window_length_bp <- chr$position_galgal2[-1] - chr$position_galgal2[-nrow(chr)]
    window_length_cM <- chr$average_cM[-1] - chr$average_cM[-nrow(chr)]
    
    data.frame(chr = unique(chr$chromosome),
               start = window_start,
               end = window_end,
               window_length_bp = window_length_bp,
               window_length_cM = window_length_cM,
               rec_rate = window_length_cM/window_length_bp * 1e6)
}
    

windows <- map_dfr(groenen2009_chr, get_marker_itervals_sorted)

windows_filtered <- filter(windows,
                           window_length_bp > 0 & window_length_bp < 1e6 &
                               window_length_cM >= 0 &
                               rec_rate < 200)

windows_filtered$window_id <- paste(windows_filtered$chr,
                                    windows_filtered$start,
                                    windows_filtered$end,
                                    sep = "_")

plot_recomb_windows_filtered <- qplot(x = (start + end)/2,
                                      y = rec_rate,
                                      data = windows_filtered) +
    facet_wrap(~chr, scale = "free_x")


bed_windows <- data.frame(chr = paste("chr", windows_filtered$chr, sep = ""),
                          start = windows_filtered$start - 1,
                          end = windows_filtered$end,
                          name = windows_filtered$window_id,
                          stringsAsFactors = FALSE)


write.table(windows_filtered,
            file = "outputs/groenen2009_windows_filtered.txt",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)


write.table(bed_windows,
            file = "outputs/groenen2009_windows.bed",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)
