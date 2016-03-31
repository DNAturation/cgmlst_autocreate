library(ape)
library(getopt)
library(vegan)

analyze <- function() {
    # calls clean_* function
    # calculates distance
    # calls reporter
}

clean_cgmlst <- function(calls) {
    # remove loci w/ missing data

    good <- sapply(calls, function(col) length(col) == sum(col >= 1))
    
    list('clean_cgmlst' = calls[, good], 'bad' = !good)
}

clean_accessory <- function(calls) {
    # remove loci with 100% presence
    # normalize to presence/absence

    good <- sapply(calls, function(col)
                   (length(col) == sum(col >= -1)) & # no uncalled loci
                   (length(col) != sum(col != 0)))   # not 100% present
    
    normalized <- mapply(function(x) if (x != 0) {1} else {0}, calls[, good])

    list('clean_accessory' = normalized, 'bad' = !good)

}

reporter <- function() {
    # writes report on changes made to call matrices
    # writes cleaned calls table
    # writes resulting distance matrix
}

calculate_distance <- function(calls, metric) {
    # calculate normalized distance matrix with hamming (default) or jaccard
    # wraps dist.gene or vegdist
    
    if (metric == 'jaccard') {
        f <- function(x) vegdist(x, 'jaccard')
    } else {
        f <- function(x) dist.gene(x, 'percentage')
    }

    f(calls)
}

options <- function() {
    # get command line arguments
}

main <- function() {

}

