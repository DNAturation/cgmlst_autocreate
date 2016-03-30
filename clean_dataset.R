library(ape)
library(getopt)
library(vegan)

analyze <- function() {
    # calls clean_* function
    # calculates distance
    # calls reporter
}

clean_cgmlst <- function() {
    # remove loci w/ missing data
}

clean_accessory <- function() {
    # remove loci with 100% presence
    # normalize to presence/absence
}

reporter <- function() {
    # writes report on changes made to call matrices
    # writes cleaned calls table
    # writes resulting distance matrix
}

calculate_distance <- function() {
    # calculate normalized distance matrix with hamming (default) or jaccard
    # wraps dist.gene or vegdist
}

options <- function() {
    # get command line arguments
}

main <- function() {

}

