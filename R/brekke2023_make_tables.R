## Make output tables for Brekke et al. (2023)

library(assertthat)
library(dplyr)
library(readxl)
library(ggplot2)
library(patchwork)
library(purrr)
library(readr)

source("R/helper_functions.R")



positions <- read_excel("data/all_chr_linkage_map_JDS.xls")

positions$sex_average <- (positions$male_cM + positions$female_cM) * 0.5



## Create map file (single file format)

map_file <- tibble(marker = positions$snpid,
                   chr = paste("chr", positions$chr, sep = ""),
                   position_bp = positions$Mb,
                   position_cM = positions$sex_average)

write.table(map_file,
            file = "outputs/brekke2023_ARS-UCD1.2.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)



positions1.3 <- read_tsv("blat/brekke2023_positions1.3.txt")

combined <- inner_join(map_file, positions1.3)

assert_that(all(combined$position_bp == combined$position_bp1.3))

assert_that(all(combined$chr == paste("chr", combined$chr1.3, sep = "")))


map_file1.3 <- tibble(marker = combined$marker,
                      chr = combined$chr,
                      position_bp = combined$position_bp1.3,
                      position_cM = combined$position_cM)


write.table(map_file1.3,
            file = "outputs/brekke2023_ARS-UCD1.3.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)


## Create map file (SHAPEIT2 format)

# windows_lifted_filtered_chr <- split(windows_lifted_filtered,
#                                      windows_lifted_filtered$chr)
# 
# 
# make_shapeit_map <- function(windows) {
#     tibble(pos = windows$end,
#            rate = windows$rec_rate,
#            cM = windows$window_end_position_cM)
# }
# 
# 
# system("mkdir outputs/brekke2023_shapeit_ARS-UCD1.2")
# 
# for (chr_ix in 1:length(windows_lifted_filtered_chr)) {
#     write.table(make_shapeit_map(windows_lifted_filtered_chr[[chr_ix]]),
#                 file = paste("outputs/brekke2023_shapeit_ARS-UCD1.2/",
#                              names(windows_lifted_filtered_chr)[chr_ix],
#                              ".txt",
#                              sep = ""),
#                 row.names = FALSE,
#                 col.names = FALSE,
#                 quote = FALSE,
#                 sep = "\t")
# }
# 
# 
# 
# 
