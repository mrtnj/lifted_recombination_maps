## Combine lifted positions with cM positions from Groenen et al. (2009)

library(GenomicRanges)
library(dplyr)
library(readr)
library(ggplot2)


## Chromosome sizes

chr_length_galgal4 <- read_tsv("data/GCF_000002315.3_Gallus_gallus-4.0_assembly_report.txt",
                               comment = "#",
                               col_names = FALSE)
chr_length_galgal4 <- chr_length_galgal4[1:33, c(10, 9)]
colnames(chr_length_galgal4) <- c("chr", "length")

chr_length_galgal4_vector <- chr_length_galgal4$length
names(chr_length_galgal4_vector) <- chr_length_galgal4$chr


## Marker positions

groenen2009 <- read_tsv("data/Supplemental_Table_S1_Groenen.txt",
                        col_types = "ccnnnnc")[,-7]

colnames(groenen2009) <- c("chromosome", "name", "female_cM", "male_cM", "average_cM",
                           "position_galgal2")



galgal4 <- read_tsv("outputs/groenen2009_hglft_galgal4.bed",
                    col_names = FALSE)


colnames(galgal4) <- c("chr", "start", "end", "name")


## Checking

a <- inner_join(groenen2009, galgal4, by = "name")
qplot(x = end, y = female_cM, data = a) + facet_wrap(~chr, scale = "free")



## Lifted windows

windows <- read_delim("outputs/groenen2009_windows_filtered.txt",
                      col_types = "cnnnnnc",
                      delim = " ")

lifted_coordinates <- read_tsv("outputs/groenen2009_windows_hglft_galgal4.bed",
                               col_names = FALSE)

colnames(lifted_coordinates) <- c("chr_lifted", "start_lifted", "end_lifted", "window_id")

windows_lifted <- inner_join(windows, lifted_coordinates,
                             by = "window_id")


plot_windows_lifted <- qplot(x = (start_lifted + end_lifted)/2/1e6,
                             y = rec_rate,
                             data = windows_lifted) +
    facet_wrap(~ chr, scale = "free_x")




## Average windows

rec_ranges <- GRanges(seqnames = windows_lifted$chr_lifted,
                      ranges = IRanges(windows_lifted$start_lifted,
                                       windows_lifted$end_lifted),
                      mcols = data.frame(rec_rate = windows_lifted$rec_rate))

windows_galgal4 <- tileGenome(seqlengths = chr_length_galgal4_vector,
                              tilewidth = 500e3,
                              cut.last.tile.in.chrom = TRUE)


rec_windows <- data.frame(as.data.frame(windows_galgal4),
                          average_rec = numeric(length(windows_galgal4)),
                          total_overlap_length = numeric(length(windows_galgal4)))

for (window_ix in 1:length(windows_galgal4)) {
    rec_in_window <- subsetByOverlaps(rec_ranges, windows_galgal4[window_ix])
    
    ## Cut all overlaps down to the start and end of the window
    window_start <- start(windows_galgal4[window_ix])
    window_end <- end(windows_galgal4[window_ix])
    
    start(rec_in_window) <- ifelse(start(rec_in_window) < window_start,
                                    window_start,
                                    start(rec_in_window))
    
    end(rec_in_window) <- ifelse(end(rec_in_window) > window_end,
                                  window_end,
                                  end(rec_in_window))
    
    total_length <- sum(width(rec_in_window))
    weighted_average <-
        sum(width(rec_in_window) * rec_in_window$mcols.rec_rate) / total_length
    
    rec_windows$average_rec[window_ix] <- weighted_average
    rec_windows$total_overlap_length[window_ix] <- total_length
}



colnames(rec_windows)[1] <- "chr"

write.table(rec_windows,
            file = "outputs/groenen2009_windows_500kbp_galgal4.txt",
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
