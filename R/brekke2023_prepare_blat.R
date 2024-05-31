

library(Biostrings)
library(purrr)
library(readr)



genome <- readBStringSet("genomes/Bos_taurus.ARS-UCD1.3.dna_sm.toplevel.fa")

names_split <- strsplit(names(genome), split = " ")

names(genome) <- paste("chr", map_chr(names_split, 1), sep = "")



map <- read_tsv("outputs/brekke2023_ARS-UCD1.2.txt")

map$flank <- ""

n_markers <- nrow(map)


for (marker_ix in 1:n_markers) {
    
    map$flank[marker_ix] <- 
        as.character(genome[[map$chr[marker_ix]]][
            (map$position_bp[marker_ix] - 200):map$position_bp[marker_ix]])
    
    
}


dir.create("blat")

write(paste("> ", map$marker, "\n", map$flank, sep = ""),
      file = "blat/brekke2023.fasta")
