## Prepare bed file from Elferink et al. (2010) Additional files 1 and 2 for lifting

library(dplyr)
library(purrr)
library(readr)

source("R/helper_functions.R")


## Read additional files

elferink_map <- read_tsv("data/Additional_file_2_Elferink_combined_sex_average.txt")
colnames(elferink_map) <- c("chr_linkage_map", "marker", "position_cM")

elferink_position <- read_tsv("data/Additional_file_1_Elferink_position.txt")
colnames(elferink_position) <- c("number", "marker", "chr_assembly", "new_chr",
                                 "position_bp", "status", "sequence")


## Combine genetic and physical positions

elferink <- inner_join(elferink_map[, c("chr_linkage_map", "marker", "position_cM")],
                       elferink_position[, c("chr_assembly", "marker", "position_bp")],
                       by = "marker")

elferink <- filter(elferink,
                   chr_linkage_map == chr_assembly &
                       chr_assembly %in% c(paste("GGA", c(1:28), sep = ""), "GGZ"))

elferink$chr <- sub(elferink$chr_assembly, pattern = "GGA|GG", replacement = "")


## Write out table

write.table(elferink,
            file = "outputs/elferink2010_marker_positions.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)


bed <- data.frame(chr = paste("chr", elferink$chr, sep = ""),
                  start = elferink$position_bp - 1,
                  end = elferink$position_bp,
                  name = elferink$marker)


write.table(bed,
            file = "outputs/elferink2010_markers.bed",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)
