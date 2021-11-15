## Prepare Dryad table from Ma et al. (2015) for lifting

library(dplyr)
library(purrr)
library(readr)
library(tibble)

source("R/helper_functions.R")



## Read Dryad files

ma_dryad <- read_delim("data/dryad_data_cattle_rmap_Ma2015.txt",
                       delim = " ")


## Convert to cM

ma_dryad$average_r <- (ma_dryad$map_f + ma_dryad$map_m) * 0.5


ma_dryad_chr <- split(ma_dryad, ma_dryad$Chr)

ma_dryad_cM <- map_dfr(ma_dryad_chr,
                       function(chr) transform(chr,
                                               position_cM = cumsum(haldane_cM(chr$average_r))))


## Format table of genetic and physical positions

ma <- tibble(marker = ma_dryad_cM$Name,
             position_cM = ma_dryad_cM$position_cM,
             chr_assembly = ma_dryad_cM$Chr,
             position_bp = ma_dryad_cM$Location)





## Marker positions

bed <- data.frame(chr = paste("chr", ma_dryad$Chr, sep = ""),
                  start = ma_dryad$Location - 1,
                  end = ma_dryad$Location,
                  name = ma_dryad$Name,
                  stringsAsFactors = FALSE)


write.table(bed,
            file = "outputs/ma2015_markers.bed",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)




write.table(ma,
            file = "outputs/ma2015_marker_positions.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)


