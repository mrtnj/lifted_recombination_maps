
library(assertthat)
library(readr)

psl <- read_table("blat/brekke2023.psl",
                  skip = 5,
                  col_names = FALSE)

colnames(psl) <- c("match",
                   "mismatch",
                   "rep_match",
                   "ns",
                   "q_gap_count",
                   "q_gap_bases",
                   "t_gap_count",
                   "t_gap_bases",
                   "strand",
                   "q_name",
                   "q_size",
                   "q_start",
                   "q_end",
                   "t_name",
                   "t_size",
                   "t_start",
                   "t_end",
                   "block_count",
                   "block_sizes",
                   "q_starts",
                   "t_starts")


full_length <- filter(psl, block_count == 1 & match == 201)

duplicated_markers <- full_length$q_name[duplicated(full_length$q_name)]


full_length_unique <- filter(full_length, ! q_name %in% duplicated_markers)


assert_that(all(full_length_unique$strand == "+"))

marker_positions <- tibble(marker = full_length$q_name,
                           chr1.3 = full_length$t_name,
                           position_bp1.3 = full_length$t_end)

options(scipen = 1e6)

write.table(marker_positions,
            file = "blat/brekke2023_positions1.3.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)
