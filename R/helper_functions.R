

haldane_cM <- function(r) -50 * log(1 - 2 * r)


## Read a tab-separate bed file (as expected from UCSC genome browser)

read_bed <- function(filename) {
    bed <- read_tsv(filename,
                    col_names = FALSE,
                    col_types = "cnnc")
    colnames(bed) <- c("chr", "start", "end", "name")
    bed
}


## Create window ID for identifying lifted windows in BED file

make_window_id <- function(windows) paste(windows$chr,
                                          windows$start,
                                          windows$end,
                                          sep = "_")


## Create windows on a chromosome

get_windows_chr <- function(chr,
                            position_bp,
                            position_cM,
                            marker_name) {
    
    n_windows <- length(chr)
    assert_that(length(unique(chr)) == 1)
    assert_that(length(position_bp) == n_windows)
    assert_that(length(position_cM) == n_windows)
    assert_that(length(marker_name) == n_windows)
    assert_that(!is.unsorted(position_bp))
        
    window_start <- position_bp[-n_windows]
    window_end <- position_bp[-1]
    window_end_marker <- marker_name[-1]
    window_length_bp <- window_end - window_start
    window_length_cM <- position_cM[-1] - position_cM[-n_windows]
    window_end_position_cM <- position_cM[-1]
    
    data.frame(chr = unique(chr),
               start = window_start,
               end = window_end,
               window_length_bp = window_length_bp,
               window_length_cM = window_length_cM,
               window_end_marker,
               window_end_position_cM,
               rec_rate = window_length_cM/window_length_bp * 1e6)
}


###############################################################################
## Functions for filtering by number of markers disagreeing with other markers

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
