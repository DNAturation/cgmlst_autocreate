

divide_calls <- function(fname) {

    calls <- read.csv(fname, stringsAsFactors = FALSE,
                      header = TRUE, row.names = 1)


    find_missing_genes <- function(name) {

        if ('0' %in% calls[, name]) {
            TRUE
        }
        else {
            FALSE
        }

    }

    missing <- sapply(colnames(calls), find_missing_genes)

    accessory <- calls[,  missing]
    core      <- calls[, !missing]

    write.csv(accessory, file = "accessory_calls.csv", quote = FALSE)
    write.csv(core, file = "core_calls.csv", quote = FALSE)

    list(colnames(calls), missing)
}
