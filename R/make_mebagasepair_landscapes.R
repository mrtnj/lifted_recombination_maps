
library(dplyr)
library(tibble)
library(readr)



cattle <- read_tsv("outputs/brekke2023_ARS-UCD1.2.txt")

assembly_report <- read_tsv("data/GCA_002263795.2_ARS-UCD1.2_assembly_report.txt",
                            comment = "#",
                            col_names = FALSE)


chr_length <- filter(assembly_report[, c(10, 9)],
                     X10 %in% paste("chr", 1:29, sep = ""))
colnames(chr_length) <- c("chr", "length")


windows <- vector(mode = "list", length = 29)
for (chr_ix in 1:29) {
    chr_name <- paste("chr", chr_ix, sep = "")
    
    len <- chr_length$length[chr_length$chr == chr_name]
    on_chr <- cattle[cattle$chr == chr_name,]
    start <- seq(from = 1, to = len, by = 1e6)
    end <- start + 1e6 - 1
    end[length(end)] <- len
    
    windows[[chr_ix]] <- tibble(chr = chr_name, start = start, end = end)
    windows[[chr_ix]]$cM_per_bp <- 0
    for (window_ix in 1:nrow(windows[[chr_ix]])) {
        in_window <- on_chr[on_chr$position_bp >= windows[[chr_ix]]$start[window_ix] &
                                on_chr$position_bp < windows[[chr_ix]]$end[window_ix],]
        
        length_cM <- max(in_window$position_cM) - min(in_window$position_cM)
        length_bp <- max(in_window$position_bp) - min(in_window$position_bp)
        
        windows[[chr_ix]]$cM_per_bp[window_ix] <- length_cM / length_bp
    }
    windows[[chr_ix]]$cM_per_bp[is.nan(windows[[chr_ix]]$cM_per_bp)] <-
        mean(windows[[chr_ix]]$cM_per_bp, na.rm = TRUE)
    
    for (window_ix in 2:(nrow(windows[[chr_ix]]) - 1)) {
        windows[[chr_ix]]$cM_per_bp[window_ix] <-
            mean(windows[[chr_ix]]$cM_per_bp[(window_ix - 1):(window_ix + 1)])
    }
}



windows_df <- Reduce(rbind, windows)



write.table(windows_df,
            file = "outputs/brekke2023_ARS-UCD1.2_Mbp_smoothed.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)
