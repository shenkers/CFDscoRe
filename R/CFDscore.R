package_state <- new.env()

.onLoad <- function(libname, pkgname) {
    package_state$activity_scores <- readr::read_rds(fs::path_package("extdata","cfd_activity_scores.rds",package='CFDscoRe'))
}
